"""End-to-end online shRNA target discovery and specificity screening."""

from __future__ import annotations

from collections.abc import Callable
from dataclasses import dataclass, replace

from genesnap_workbench.domain.shrna import BlastScreenStatus, ShRNACandidate
from genesnap_workbench.sequence_core.shrna import (
    build_shrna_oligo_pair,
    select_initial_candidates,
)

from .broad_gpp import BroadGPPClient, BroadGPPError, BroadHairpinCandidate
from .ncbi_blast import BlastClassification, NCBIBlastClient, NCBIBlastError


ProgressCallback = Callable[[str], None]


@dataclass(frozen=True, slots=True)
class ShRNAOnlineDesignResult:
    selected_candidates: tuple[ShRNACandidate, ...]
    candidate_pool: tuple[ShRNACandidate, ...]
    requires_manual_confirmation: bool
    notes: tuple[str, ...]


class ShRNAOnlineDesigner:
    def __init__(
        self,
        *,
        broad_client: BroadGPPClient | None = None,
        blast_client: NCBIBlastClient | None = None,
        blast_batch_size: int = 20,
    ) -> None:
        if blast_batch_size <= 0:
            raise ValueError("blast_batch_size must be positive")
        self.broad_client = broad_client or BroadGPPClient()
        self.blast_client = blast_client or NCBIBlastClient()
        self.blast_batch_size = blast_batch_size

    @staticmethod
    def _to_domain(item: BroadHairpinCandidate) -> ShRNACandidate:
        return ShRNACandidate(
            candidate_id=f"broad-{item.source_rank}",
            target_sequence=item.target_sequence,
            start_position=item.start_position,
            intrinsic_score=item.intrinsic_score,
            source_rank=item.source_rank,
            blast_status=BlastScreenStatus.PENDING,
            blast_note="Broad GPP 候选，等待 NCBI BLAST",
        )

    def design(
        self,
        *,
        cds_sequence: str,
        gene_symbol: str,
        species: str,
        target_count: int,
        progress: ProgressCallback | None = None,
    ) -> ShRNAOnlineDesignResult:
        emit = progress or (lambda message: None)
        emit("正在向 Broad GPP 提交 CDS 并读取候选 target…")
        broad_candidates = self.broad_client.design_hairpins(cds_sequence)
        pool = tuple(self._to_domain(item) for item in broad_candidates)
        broad_by_sequence = {item.target_sequence: item for item in broad_candidates}
        initial = select_initial_candidates(pool, target_count=target_count).selected
        ordered_for_screening = initial + tuple(
            item for item in pool if item.candidate_id not in {selected.candidate_id for selected in initial}
        )
        updated = {item.candidate_id: item for item in pool}
        notes: list[str] = [f"Broad GPP 返回 {len(pool)} 条候选"]
        blast_unavailable = False

        for start in range(0, len(ordered_for_screening), self.blast_batch_size):
            batch = tuple(
                updated[item.candidate_id]
                for item in ordered_for_screening[start : start + self.blast_batch_size]
            )
            emit(
                f"正在通过 NCBI BLAST 检查候选 {start + 1}-"
                f"{min(start + len(batch), len(ordered_for_screening))}…",
            )
            try:
                classifications = self.blast_client.screen_targets(
                    batch,
                    expected_gene_symbol=gene_symbol,
                    species=species,
                )
            except NCBIBlastError as error:
                blast_unavailable = True
                notes.append(str(error))
                for item in initial:
                    updated[item.candidate_id] = replace(
                        item,
                        blast_status=BlastScreenStatus.UNAVAILABLE,
                        blast_note=str(error),
                    )
                break
            for item in batch:
                result = classifications[item.target_sequence]
                updated[item.candidate_id] = replace(
                    item,
                    blast_status=result.status,
                    first_offtarget_gene=result.first_offtarget_gene,
                    first_offtarget_mismatches=result.first_offtarget_mismatches,
                    blast_note=result.note,
                )
            selected = self._resolve_selected(initial, tuple(updated.values()), target_count)
            if len(selected) == target_count:
                break

        selected = (
            tuple(updated[item.candidate_id] for item in initial)
            if blast_unavailable
            else self._resolve_selected(initial, tuple(updated.values()), target_count)
        )
        if len(selected) < target_count:
            notes.append(f"自动 BLAST 后仅有 {len(selected)} 条候选通过，未达到 {target_count} 条")
        selected_with_oligos: list[ShRNACandidate] = []
        oligo_mismatch = False
        for index, item in enumerate(selected, start=1):
            emit(f"正在读取 Broad oligo 详情 {index}/{len(selected)}…")
            local = build_shrna_oligo_pair(
                gene_symbol=gene_symbol,
                target_no=index,
                target_id=item.candidate_id,
                target_sequence=item.target_sequence,
            )
            broad_item = broad_by_sequence[item.target_sequence]
            try:
                broad_oligos = self.broad_client.fetch_oligos(broad_item)
            except BroadGPPError as error:
                selected_with_oligos.append(
                    replace(
                        item,
                        forward_oligo_sequence=local.forward_sequence,
                        reverse_oligo_sequence=local.reverse_sequence,
                        oligo_source="local_plko_fallback",
                        oligo_comparison_note=f"Broad oligo 详情不可用：{error}",
                    ),
                )
                continue
            matches_local = (
                broad_oligos.forward_sequence == local.forward_sequence
                and broad_oligos.reverse_sequence == local.reverse_sequence
            )
            if not matches_local:
                oligo_mismatch = True
            selected_with_oligos.append(
                replace(
                    item,
                    forward_oligo_sequence=broad_oligos.forward_sequence,
                    reverse_oligo_sequence=broad_oligos.reverse_sequence,
                    oligo_source="broad_gpp",
                    oligo_comparison_note=(
                        "Broad Full sequence 与本地 pLKO 规则一致"
                        if matches_local
                        else "WARNING：Broad Full sequence 与本地 pLKO 规则不一致，已采用 Broad"
                    ),
                ),
            )
        final_pool = tuple(updated[item.candidate_id] for item in pool)
        requires_confirmation = (
            blast_unavailable
            or len(selected_with_oligos) < target_count
            or oligo_mismatch
        )
        return ShRNAOnlineDesignResult(
            selected_candidates=tuple(selected_with_oligos),
            candidate_pool=final_pool,
            requires_manual_confirmation=requires_confirmation,
            notes=tuple(notes),
        )

    @staticmethod
    def _resolve_selected(
        initial: tuple[ShRNACandidate, ...],
        pool: tuple[ShRNACandidate, ...],
        target_count: int,
    ) -> tuple[ShRNACandidate, ...]:
        by_id = {item.candidate_id: item for item in pool}
        selected = [
            by_id[item.candidate_id]
            for item in initial
            if by_id[item.candidate_id].blast_status is BlastScreenStatus.PASS
        ]
        initial_ids = {item.candidate_id for item in initial}
        ranked = sorted(pool, key=lambda item: (-item.intrinsic_score, item.source_rank))
        for item in ranked:
            if len(selected) == target_count:
                break
            if item.candidate_id in initial_ids or item.blast_status is not BlastScreenStatus.PASS:
                continue
            selected.append(item)
        return tuple(selected)
