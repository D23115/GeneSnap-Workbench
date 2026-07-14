"""表达类全长、截短和缺失构建的通用序列变换。"""

from __future__ import annotations

import re
from dataclasses import replace
from datetime import datetime

from genesnap_workbench.domain.expression import (
    ExpressionConstruct,
    ExpressionConstructKind,
    ExpressionDesignInput,
    ExpressionDesignConfirmation,
    ExpressionDesignVersion,
    ExpressionPCRFragment,
    ExpressionMutation,
    ExpressionVectorRules,
)
from Bio.Data import CodonTable
from Bio.Seq import Seq

from .dna import sha256_sequence


STOP_CODONS = {"TAA", "TAG", "TGA"}


def confirm_expression_design(
    design: ExpressionDesignVersion,
    *,
    confirmation_id: str,
    reason: str,
    actor: str,
    occurred_at: datetime,
) -> ExpressionDesignVersion:
    clean_reason = reason.strip()
    if not clean_reason:
        raise ValueError("表达设计人工确认必须填写原因")
    if not design.requires_confirmation:
        return design
    record = ExpressionDesignConfirmation(
        confirmation_id=confirmation_id,
        reason=clean_reason,
        actor=actor.strip() or "当前用户",
        occurred_at=occurred_at,
    )
    return replace(
        design,
        requires_confirmation=False,
        confirmation_history=design.confirmation_history + (record,),
    )


def _split_source_cds(source_cds: str) -> tuple[str, str | None]:
    if not source_cds.startswith("ATG"):
        raise ValueError("source CDS must start with ATG")
    if len(source_cds) % 3:
        raise ValueError("source CDS length must be divisible by 3")
    terminal = source_cds[-3:] if source_cds[-3:] in STOP_CODONS else None
    core = source_cds[:-3] if terminal else source_cds
    for offset in range(0, len(core), 3):
        if core[offset : offset + 3] in STOP_CODONS:
            raise ValueError(f"source CDS contains internal stop codon at nt {offset + 1}")
    return core, terminal


def _normalized_request(line: str, gene_symbol: str) -> str:
    request = line.strip()
    prefix = f"{gene_symbol}-"
    if request.casefold().startswith(prefix.casefold()):
        request = request[len(prefix) :]
    return request


def _parse_request(
    line: str,
    gene_symbol: str,
) -> tuple[
    ExpressionConstructKind,
    int | None,
    int | None,
    tuple[str, ...],
] | None:
    request = _normalized_request(line, gene_symbol)
    if request.upper() in {"FL", "WT", "FULL", "FULL-LENGTH"}:
        return ExpressionConstructKind.FULL_LENGTH, None, None, ()
    deletion = re.fullmatch(r"(?:Δ|DEL)(\d+)-(\d+)(?:AA|A)?", request, re.IGNORECASE)
    if deletion:
        return ExpressionConstructKind.DELETION, int(deletion.group(1)), int(deletion.group(2)), ()
    truncation = re.fullmatch(r"(\d+)-(\d+)(?:AA|A)?", request, re.IGNORECASE)
    if truncation:
        return ExpressionConstructKind.TRUNCATION, int(truncation.group(1)), int(truncation.group(2)), ()
    mutation_items = tuple(item.strip().upper() for item in request.split("+") if item.strip())
    if mutation_items and len(mutation_items) <= 2 and all(
        re.fullmatch(r"[A-Z*]\d+[A-Z*]", item) for item in mutation_items
    ):
        return ExpressionConstructKind.MUTATION, None, None, mutation_items
    return None


def _construct_name(
    gene_symbol: str,
    kind: ExpressionConstructKind,
    start: int | None,
    end: int | None,
    mutation_notations: tuple[str, ...] = (),
) -> str:
    if kind is ExpressionConstructKind.FULL_LENGTH:
        return f"{gene_symbol}-FL"
    if kind is ExpressionConstructKind.TRUNCATION:
        return f"{gene_symbol}-{start}-{end}aa"
    if kind is ExpressionConstructKind.MUTATION:
        return f"{gene_symbol}-{'-'.join(mutation_notations)}"
    return f"{gene_symbol}-Δ{start}-{end}"


def _transform_core(
    core: str,
    kind: ExpressionConstructKind,
    start: int | None,
    end: int | None,
) -> tuple[str, bool]:
    amino_acids = len(core) // 3
    if kind is ExpressionConstructKind.FULL_LENGTH:
        return core, False
    assert start is not None and end is not None
    if start < 1 or end < start or end > amino_acids:
        raise ValueError(
            f"construct amino-acid interval {start}-{end} exceeds 1-{amino_acids}",
        )
    if kind is ExpressionConstructKind.TRUNCATION:
        transformed = core[(start - 1) * 3 : end * 3]
    else:
        transformed = core[: (start - 1) * 3] + core[end * 3 :]
    reintroduced = not transformed.startswith("ATG")
    if reintroduced:
        transformed = "ATG" + transformed
    return transformed, reintroduced


def _codons_by_amino_acid() -> dict[str, tuple[str, ...]]:
    table = CodonTable.unambiguous_dna_by_name["Standard"]
    grouped: dict[str, list[str]] = {}
    for codon, amino_acid in table.forward_table.items():
        grouped.setdefault(amino_acid, []).append(codon)
    grouped["*"] = list(table.stop_codons)
    return {key: tuple(sorted(value)) for key, value in grouped.items()}


CODONS_BY_AMINO_ACID = _codons_by_amino_acid()


def _apply_mutations(
    core: str,
    notations: tuple[str, ...],
) -> tuple[str, tuple[ExpressionMutation, ...]]:
    translated = str(Seq(core).translate())
    sequence = list(core)
    records: list[ExpressionMutation] = []
    used_positions: set[int] = set()
    for notation in notations:
        match = re.fullmatch(r"([A-Z*])(\d+)([A-Z*])", notation)
        assert match is not None
        original_aa, position_text, new_aa = match.groups()
        position = int(position_text)
        if position < 1 or position > len(translated):
            raise ValueError(f"mutation {notation} exceeds protein length {len(translated)}")
        if position in used_positions:
            raise ValueError(f"multiple mutations target amino acid {position}")
        used_positions.add(position)
        observed = translated[position - 1]
        if observed != original_aa:
            raise ValueError(
                f"mutation {notation} expects {original_aa} at aa {position}, observed {observed}",
            )
        candidates = CODONS_BY_AMINO_ACID.get(new_aa)
        if not candidates:
            raise ValueError(f"unsupported target amino acid: {new_aa}")
        start = (position - 1) * 3
        original_codon = core[start : start + 3]
        new_codon = min(
            candidates,
            key=lambda codon: (
                sum(left != right for left, right in zip(original_codon, codon, strict=True)),
                abs((codon.count("G") + codon.count("C")) - 2),
                codon,
            ),
        )
        sequence[start : start + 3] = new_codon
        records.append(
            ExpressionMutation(
                notation=notation,
                amino_acid_position=position,
                original_amino_acid=original_aa,
                new_amino_acid=new_aa,
                original_codon=original_codon,
                new_codon=new_codon,
            ),
        )
    return "".join(sequence), tuple(records)


def _plan_fragments(
    insert_sequence: str,
    max_bp: int,
) -> tuple[ExpressionPCRFragment, ...]:
    if len(insert_sequence) <= max_bp:
        intervals = ((0, len(insert_sequence)),)
    else:
        midpoint = len(insert_sequence) // 2
        intervals = ((0, midpoint), (midpoint, len(insert_sequence)))
    return tuple(
        ExpressionPCRFragment(
            fragment_no=index,
            start=start,
            end=end,
            sequence=insert_sequence[start:end],
        )
        for index, (start, end) in enumerate(intervals, start=1)
    )


def create_expression_design(
    design_input: ExpressionDesignInput,
    rules: ExpressionVectorRules,
    *,
    design_version_id: str,
    created_at: datetime,
) -> ExpressionDesignVersion:
    core, original_stop = _split_source_cds(design_input.source_cds)
    constructs: list[ExpressionConstruct] = []
    unparsed: list[str] = []
    warnings: list[str] = []
    for line in design_input.construct_lines:
        parsed = _parse_request(line, design_input.gene_symbol)
        if parsed is None:
            unparsed.append(line)
            continue
        kind, start, end, mutation_notations = parsed
        if kind is ExpressionConstructKind.MUTATION:
            coding_core, mutations = _apply_mutations(core, mutation_notations)
            reintroduced = False
        else:
            coding_core, reintroduced = _transform_core(core, kind, start, end)
            mutations = ()
        preserve_stop = rules.stop_codon_rule == "preserve"
        terminal_stop = original_stop or "TAA"
        coding_sequence = coding_core + (terminal_stop if preserve_stop else "")
        insert_sequence = rules.kozak_sequence + coding_sequence
        name = _construct_name(
            design_input.gene_symbol,
            kind,
            start,
            end,
            mutation_notations,
        )
        constructs.append(
            ExpressionConstruct(
                construct_id=f"{design_input.project_id}-construct-{len(constructs) + 1}",
                construct_name=name,
                request_line=line,
                kind=kind,
                coding_sequence=coding_sequence,
                insert_sequence=insert_sequence,
                start_codon_reintroduced=reintroduced,
                terminal_stop_present=coding_sequence[-3:] in STOP_CODONS,
                c_terminal_fusion_name=rules.c_terminal_fusion_name,
                mutations=mutations,
                fragments=_plan_fragments(
                    insert_sequence,
                    rules.single_fragment_max_bp,
                ),
            ),
        )
    if unparsed:
        warnings.append(f"{len(unparsed)} 行构建需求无法自动解析，需要人工确认")
    mutation_count = sum(bool(item.mutations) for item in constructs)
    if mutation_count:
        warnings.append(f"{mutation_count} 个点突变构建需要人工确认突变序列")
    return ExpressionDesignVersion(
        design_version_id=design_version_id,
        project_id=design_input.project_id,
        gene_symbol=design_input.gene_symbol,
        species=design_input.species,
        transcript_accession=design_input.transcript_accession,
        source_cds_checksum=sha256_sequence(design_input.source_cds),
        protocol_version_id=rules.protocol_version_id,
        created_at=created_at,
        constructs=tuple(constructs),
        unparsed_lines=tuple(unparsed),
        design_warnings=tuple(warnings),
        requires_confirmation=bool(unparsed or mutation_count),
    )
