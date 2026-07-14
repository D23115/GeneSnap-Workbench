"""pLKO/shRNA 候选选择和 oligo 设计纯逻辑。"""

from __future__ import annotations

from datetime import datetime

from genesnap_workbench.domain.shrna import (
    BlastScreenStatus,
    ShRNACandidate,
    ShRNABlastResolution,
    ShRNADesignInput,
    ShRNADesignVersion,
    ShRNAOligoPair,
    ShRNATargetDesign,
    ShRNATargetSelectionResult,
)

from .dna import reverse_complement, sha256_sequence
from genesnap_workbench.vector_library.models import ShRNAVectorProtocol, VectorRecord
from genesnap_workbench.vector_library.shrna import (
    simulate_shrna_plasmid,
    validate_shrna_protocol,
)


SHRNA_SELECTION_RULE_VERSION = "shrna-selection-v1"
PLKO_OLIGO_RULE_VERSION = "plko-hairpin-oligo-v1"


def _rank_candidates(candidates: tuple[ShRNACandidate, ...]) -> tuple[ShRNACandidate, ...]:
    unique: dict[str, ShRNACandidate] = {}
    for item in candidates:
        current = unique.get(item.target_sequence)
        if current is None or (item.intrinsic_score, -item.source_rank) > (
            current.intrinsic_score,
            -current.source_rank,
        ):
            unique[item.target_sequence] = item
    return tuple(
        sorted(
            unique.values(),
            key=lambda item: (
                -item.intrinsic_score,
                item.source_rank,
                item.start_position if item.start_position is not None else 10**12,
                item.candidate_id,
            ),
        ),
    )


def select_initial_candidates(
    candidates: tuple[ShRNACandidate, ...],
    *,
    target_count: int,
    minimum_spacing_bp: int = 100,
) -> ShRNATargetSelectionResult:
    """按得分选择首轮 target，并尽量保持严格大于 100 bp 的间距。"""
    if target_count <= 0:
        raise ValueError("target_count must be positive")
    if minimum_spacing_bp < 0:
        raise ValueError("minimum_spacing_bp must not be negative")
    ranked = _rank_candidates(candidates)
    eligible = tuple(
        item for item in ranked if item.blast_status is not BlastScreenStatus.FAIL
    )
    if len(eligible) < target_count:
        raise ValueError("可用 shRNA 候选数量少于需求数量")

    selected: list[ShRNACandidate] = []
    for item in eligible:
        if item.start_position is not None and all(
            existing.start_position is not None
            and abs(item.start_position - existing.start_position) > minimum_spacing_bp
            for existing in selected
        ):
            selected.append(item)
            if len(selected) == target_count:
                break

    spacing_relaxed = len(selected) < target_count
    if spacing_relaxed:
        selected_sequences = {item.target_sequence for item in selected}
        for item in eligible:
            if item.target_sequence in selected_sequences:
                continue
            selected.append(item)
            selected_sequences.add(item.target_sequence)
            if len(selected) == target_count:
                break

    selected_sequences = {item.target_sequence for item in selected}
    return ShRNATargetSelectionResult(
        selected=tuple(selected),
        remaining=tuple(
            item for item in eligible if item.target_sequence not in selected_sequences
        ),
        spacing_relaxed=spacing_relaxed,
        minimum_spacing_bp=minimum_spacing_bp,
    )


def advance_blast_selection(
    current_selected: tuple[ShRNACandidate, ...],
    candidate_pool: tuple[ShRNACandidate, ...],
    *,
    target_count: int,
) -> ShRNABlastResolution:
    """处理一轮 BLAST 结果；失败后按得分补选且不再强制间距。"""
    if len(current_selected) != target_count:
        raise ValueError("current_selected 数量与 target_count 不一致")
    if any(
        item.blast_status is BlastScreenStatus.UNAVAILABLE
        for item in current_selected
    ):
        return ShRNABlastResolution(
            selected=current_selected,
            needs_screening=(),
            completed=False,
            requires_confirmation=True,
            note="自动 BLAST 不可用，保留首轮 target 并等待人工确认",
        )

    retained = [
        item
        for item in current_selected
        if item.blast_status is not BlastScreenStatus.FAIL
    ]
    selected_ids = {item.candidate_id for item in current_selected}
    if len(retained) < target_count:
        replacements = (
            item
            for item in _rank_candidates(candidate_pool)
            if item.candidate_id not in selected_ids
            and item.blast_status is not BlastScreenStatus.FAIL
        )
        for item in replacements:
            retained.append(item)
            if len(retained) == target_count:
                break
    if len(retained) < target_count:
        raise ValueError("BLAST 淘汰后没有足够候选可补选")

    selected = tuple(retained)
    needs_screening = tuple(
        item
        for item in selected
        if item.blast_status is BlastScreenStatus.PENDING
    )
    completed = all(
        item.blast_status
        in {BlastScreenStatus.PASS, BlastScreenStatus.MANUALLY_ACCEPTED}
        for item in selected
    )
    return ShRNABlastResolution(
        selected=selected,
        needs_screening=needs_screening,
        completed=completed,
        requires_confirmation=False,
        note=(
            "全部 target 已通过 BLAST"
            if completed
            else "已按得分补选 target，需要继续 BLAST"
        ),
    )


def build_shrna_oligo_pair(
    *,
    gene_symbol: str,
    target_no: int,
    target_id: str,
    target_sequence: str,
) -> ShRNAOligoPair:
    target = target_sequence.strip().upper()
    reverse_target = reverse_complement(target)
    return ShRNAOligoPair(
        target_id=target_id,
        forward_name=f"{gene_symbol}-{target_no}-F",
        forward_sequence=f"CCGG{target}CTCGAG{reverse_target}TTTTTG",
        reverse_name=f"{gene_symbol}-{target_no}-R",
        reverse_sequence=f"AATTCAAAAA{target}CTCGAG{reverse_target}",
    )


def create_shrna_design(
    design_input: ShRNADesignInput,
    selected_candidates: tuple[ShRNACandidate, ...],
    vector: VectorRecord,
    protocol: ShRNAVectorProtocol,
    *,
    design_version_id: str,
    created_at: datetime,
) -> ShRNADesignVersion:
    if len(selected_candidates) != design_input.target_count:
        raise ValueError("selected target 数量与需求数量不一致")
    if any(item.blast_status is BlastScreenStatus.FAIL for item in selected_candidates):
        raise ValueError("BLAST 明确失败的 target 不能进入正式设计")
    if design_input.vector_protocol_version_id != protocol.protocol_version_id:
        raise ValueError("设计输入与载体 protocol 版本不一致")
    validation = validate_shrna_protocol(vector, protocol)
    if not validation.is_valid or protocol.status != "enabled":
        reasons = "; ".join(issue.message for issue in validation.errors)
        raise ValueError(reasons or "shRNA 载体 protocol 尚未启用")

    pending_blast = tuple(
        item
        for item in selected_candidates
        if item.blast_status in {
            BlastScreenStatus.PENDING,
            BlastScreenStatus.UNAVAILABLE,
        }
    )
    blast_warnings = (
        (f"{len(pending_blast)} 条 target 未完成自动 BLAST，需要人工确认",)
        if pending_blast
        else ()
    )
    oligo_warnings = tuple(
        item.oligo_comparison_note
        for item in selected_candidates
        if item.oligo_comparison_note
        and item.oligo_comparison_note.startswith("WARNING")
    )
    warnings = blast_warnings + oligo_warnings

    targets: list[ShRNATargetDesign] = []
    for target_no, candidate in enumerate(selected_candidates, start=1):
        target_id = f"{design_input.project_id}-target-{target_no}"
        local_oligos = build_shrna_oligo_pair(
            gene_symbol=design_input.gene_symbol,
            target_no=target_no,
            target_id=target_id,
            target_sequence=candidate.target_sequence,
        )
        oligos = (
            ShRNAOligoPair(
                target_id=target_id,
                forward_name=local_oligos.forward_name,
                forward_sequence=candidate.forward_oligo_sequence,
                reverse_name=local_oligos.reverse_name,
                reverse_sequence=candidate.reverse_oligo_sequence,
            )
            if candidate.forward_oligo_sequence is not None
            and candidate.reverse_oligo_sequence is not None
            else local_oligos
        )
        targets.append(
            ShRNATargetDesign(
                target_id=target_id,
                target_no=target_no,
                candidate=candidate,
                oligos=oligos,
                clone_names=tuple(
                    f"{design_input.gene_symbol}-{target_no}-{clone_no}"
                    for clone_no in range(1, design_input.clones_per_target + 1)
                ),
            ),
        )

    return ShRNADesignVersion(
        design_version_id=design_version_id,
        project_id=design_input.project_id,
        gene_symbol=design_input.gene_symbol,
        species=design_input.species,
        vector_record_id=vector.vector_record_id,
        vector_checksum=vector.normalized_circular_sha256,
        vector_protocol_version_id=design_input.vector_protocol_version_id,
        created_at=created_at,
        cds_checksum=sha256_sequence(design_input.cds_sequence),
        transcript_accession=design_input.transcript_accession,
        gene_id=design_input.gene_id,
        ccds_id=design_input.ccds_id,
        target_count=design_input.target_count,
        clones_per_target=design_input.clones_per_target,
        targets=tuple(targets),
        plasmid_simulations=tuple(
            simulate_shrna_plasmid(
                vector,
                protocol,
                target_id=target.target_id,
                forward_oligo=target.oligos.forward_sequence,
            )
            for target in targets
        ),
        rule_versions=(SHRNA_SELECTION_RULE_VERSION, PLKO_OLIGO_RULE_VERSION),
        design_warnings=warnings,
        requires_confirmation=bool(pending_blast),
    )
