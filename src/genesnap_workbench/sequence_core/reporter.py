"""GL002 promoter WT、逐级删除和区域突变的序列变换。"""

from __future__ import annotations

import re
from dataclasses import replace
from datetime import datetime

from genesnap_workbench.domain.reporter import (
    PromoterMutationDefinition,
    ReporterConstruct,
    ReporterDesignConfirmation,
    ReporterDesignInput,
    ReporterDesignVersion,
)

from .dna import sha256_sequence


def confirm_reporter_design(
    design: ReporterDesignVersion,
    *,
    confirmation_id: str,
    reason: str,
    actor: str,
    occurred_at: datetime,
) -> ReporterDesignVersion:
    clean_reason = reason.strip()
    if not clean_reason:
        raise ValueError("reporter 设计人工确认必须填写原因")
    if not design.requires_confirmation:
        return design
    confirmation = ReporterDesignConfirmation(
        confirmation_id=confirmation_id,
        reason=clean_reason,
        actor=actor.strip() or "当前用户",
        occurred_at=occurred_at,
    )
    return replace(
        design,
        requires_confirmation=False,
        confirmation_history=design.confirmation_history + (confirmation,),
    )


def _parse_mutations(
    design_input: ReporterDesignInput,
) -> dict[str, PromoterMutationDefinition]:
    definitions = {}
    for line in design_input.mutation_definitions:
        match = re.fullmatch(
            r"([A-Za-z][A-Za-z0-9_-]*):\s*(\d+)\s*-\s*(\d+)\s*=\s*([ACGTacgt]+)",
            line,
        )
        if match is None:
            raise ValueError(f"无法解析突变定义：{line}")
        name = match.group(1)
        start = int(match.group(2))
        end = int(match.group(3))
        replacement = match.group(4).upper()
        if start < 1 or end < start or end > len(design_input.promoter_sequence):
            raise ValueError(f"突变 {name} 的坐标超出 promoter 范围")
        if name in definitions:
            raise ValueError(f"突变名称重复：{name}")
        definitions[name] = PromoterMutationDefinition(
            name=name,
            start=start,
            end=end,
            replacement_sequence=replacement,
            original_sequence=design_input.promoter_sequence[start - 1 : end],
        )
    return definitions


def _parse_construct_line(
    line: str,
    definitions: dict[str, PromoterMutationDefinition],
) -> tuple[int | None, tuple[str, ...]]:
    tokens = tuple(item.strip() for item in line.split("+") if item.strip())
    if not tokens:
        raise ValueError("构建需求不能为空")
    retained_length = None
    mutation_names = []
    for token in tokens:
        if token.upper() == "WT":
            if retained_length is not None:
                raise ValueError(f"构建需求重复定义 promoter 长度：{line}")
            retained_length = None
        elif re.fullmatch(r"P\d+", token, re.IGNORECASE):
            if retained_length is not None:
                raise ValueError(f"构建需求重复定义 promoter 长度：{line}")
            retained_length = int(token[1:])
        elif token in definitions:
            mutation_names.append(token)
        else:
            raise ValueError(f"无法解析构建需求：{line}")
    return retained_length, tuple(mutation_names)


def create_reporter_design(
    design_input: ReporterDesignInput,
    *,
    protocol_version_id: str,
    design_version_id: str,
    created_at: datetime,
) -> ReporterDesignVersion:
    definitions = _parse_mutations(design_input)
    constructs = []
    mutation_construct_count = 0
    full_length = len(design_input.promoter_sequence)
    for line in design_input.construct_lines:
        requested_length, mutation_names = _parse_construct_line(line, definitions)
        retained_length = requested_length or full_length
        if retained_length < 1 or retained_length > full_length:
            raise ValueError(f"P{retained_length} 超出输入 promoter 长度")
        retained_start = full_length - retained_length
        selected_definitions = tuple(definitions[name] for name in mutation_names)
        for mutation in selected_definitions:
            if mutation.start - 1 < retained_start:
                raise ValueError(
                    f"突变 {mutation.name} 不在 P{retained_length} 保留范围内",
                )
        sequence = design_input.promoter_sequence[retained_start:]
        for mutation in sorted(selected_definitions, key=lambda item: item.start, reverse=True):
            local_start = mutation.start - 1 - retained_start
            local_end = mutation.end - retained_start
            sequence = (
                sequence[:local_start]
                + mutation.replacement_sequence
                + sequence[local_end:]
            )
        if mutation_names:
            mutation_construct_count += 1
        length_label = "WT" if requested_length is None else str(retained_length)
        mutation_label = f"-{'-'.join(mutation_names)}" if mutation_names else ""
        constructs.append(
            ReporterConstruct(
                construct_id=f"{design_input.project_id}-construct-{len(constructs) + 1}",
                construct_name=(
                    f"{design_input.gene_symbol}-promoter-{length_label}{mutation_label}"
                ),
                request_line=line,
                retained_promoter_length=retained_length,
                retained_source_start=retained_start,
                mutation_names=mutation_names,
                insert_sequence=sequence,
            ),
        )
    warnings = ()
    if mutation_construct_count:
        warnings = (
            f"{mutation_construct_count} 个 promoter 突变构建需要人工确认最终替换序列",
        )
    return ReporterDesignVersion(
        design_version_id=design_version_id,
        project_id=design_input.project_id,
        gene_symbol=design_input.gene_symbol,
        species=design_input.species,
        transcript_accession=design_input.transcript_accession,
        promoter_source_checksum=sha256_sequence(design_input.promoter_sequence),
        protocol_version_id=protocol_version_id,
        created_at=created_at,
        mutation_definitions=tuple(definitions.values()),
        constructs=tuple(constructs),
        design_warnings=warnings,
        requires_confirmation=bool(mutation_construct_count),
    )
