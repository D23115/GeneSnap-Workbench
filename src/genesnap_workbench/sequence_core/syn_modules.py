"""Route SYN targets into a single pool or contiguous synthesis modules."""

from __future__ import annotations

from dataclasses import dataclass
from math import ceil

from genesnap_workbench.domain.syn import (
    SYNModule,
    SYNModulePlan,
    SYNRoute,
    SYNSequenceQCResult,
)

from .dna import normalize_dna, sha256_sequence


@dataclass(frozen=True, slots=True)
class SYNModuleRules:
    single_pool_default_max_bp: int = 1000
    single_pool_warning_max_bp: int = 1200
    module_min_bp: int = 500
    module_target_max_bp: int = 900
    module_hard_max_bp: int = 1000
    direct_assembly_module_max: int = 4
    warning_module_count: int = 5
    module_overlap_min_bp: int = 20
    module_overlap_target_bp: int = 25
    module_overlap_max_bp: int = 30


def _module_sizes(sequence_length: int, rules: SYNModuleRules) -> tuple[int, ...]:
    module_count = ceil(sequence_length / rules.module_target_max_bp)
    base_size, remainder = divmod(sequence_length, module_count)
    sizes = tuple(
        base_size + (1 if index < remainder else 0)
        for index in range(module_count)
    )
    if any(
        size < rules.module_min_bp or size > rules.module_hard_max_bp
        for size in sizes
    ):
        raise ValueError("Unable to partition target within module size constraints")
    return sizes


def _boundary_risk_score(
    coordinate: int,
    qc_result: SYNSequenceQCResult | None,
) -> tuple[int, int]:
    if qc_result is None:
        return (0, 0)
    covering = [
        risk
        for risk in qc_result.risks
        if risk.start <= coordinate < risk.end
    ]
    return (
        sum(risk.severity == "high_risk" for risk in covering),
        sum(risk.severity == "warning" for risk in covering),
    )


def _is_unique(sequence: str, motif: str) -> bool:
    first = sequence.find(motif)
    return first >= 0 and sequence.find(motif, first + 1) < 0


def _find_unique_module_overlap(
    sequence: str,
    boundary: int,
    rules: SYNModuleRules,
) -> tuple[int, int] | None:
    lengths = sorted(
        range(rules.module_overlap_min_bp, rules.module_overlap_max_bp + 1),
        key=lambda length: (
            abs(length - rules.module_overlap_target_bp),
            length,
        ),
    )
    for length in lengths:
        centered_start = boundary - length // 2
        for offset in (0, -1, 1, -2, 2, -3, 3):
            start = centered_start + offset
            end = start + length
            if start < 0 or end > len(sequence) or not (start < boundary < end):
                continue
            if _is_unique(sequence, sequence[start:end]):
                return (start, end)
    return None


def _adjust_boundaries_for_qc(
    sequence: str,
    sizes: tuple[int, ...],
    rules: SYNModuleRules,
    qc_result: SYNSequenceQCResult | None,
) -> tuple[tuple[int, ...], tuple[str, ...], tuple[tuple[int, int], ...]]:
    sequence_length = len(sequence)
    desired_boundaries: list[int] = []
    running = 0
    for size in sizes[:-1]:
        running += size
        desired_boundaries.append(running)

    selected: list[int] = []
    reasons: list[str] = []
    overlaps: list[tuple[int, int]] = []
    overlap_cache: dict[int, tuple[int, int] | None] = {}
    previous = 0
    module_count = len(sizes)
    for index, desired in enumerate(desired_boundaries, start=1):
        remaining_modules = module_count - index
        lower = max(
            previous + rules.module_min_bp,
            sequence_length - remaining_modules * rules.module_hard_max_bp,
        )
        upper = min(
            previous + rules.module_target_max_bp,
            sequence_length - remaining_modules * rules.module_min_bp,
        )
        boundary: int | None = None
        overlap: tuple[int, int] | None = None
        best_fallback: tuple[tuple[int, int, int, int], int, tuple[int, int]] | None = None
        maximum_distance = max(desired - lower, upper - desired)
        for distance in range(maximum_distance + 1):
            for coordinate in (desired - distance, desired + distance):
                if coordinate < lower or coordinate > upper:
                    continue
                if coordinate not in overlap_cache:
                    overlap_cache[coordinate] = _find_unique_module_overlap(
                        sequence,
                        coordinate,
                        rules,
                    )
                candidate_overlap = overlap_cache[coordinate]
                if candidate_overlap is None:
                    continue
                high_risk, warning = _boundary_risk_score(
                    coordinate,
                    qc_result,
                )
                if high_risk == 0 and warning == 0:
                    boundary = coordinate
                    overlap = candidate_overlap
                    break
                fallback_score = (high_risk, warning, distance, coordinate)
                if best_fallback is None or fallback_score < best_fallback[0]:
                    best_fallback = (
                        fallback_score,
                        coordinate,
                        candidate_overlap,
                    )
            if boundary is not None:
                break
        if boundary is None and best_fallback is not None:
            _, boundary, overlap = best_fallback
        if boundary is None or overlap is None:
            raise ValueError("20-30 bp 范围内找不到唯一模块同源区")
        selected.append(boundary)
        overlaps.append(overlap)
        reasons.append(
            "low_risk_adjusted_boundary"
            if boundary != desired
            else "length_balanced_boundary"
        )
        previous = boundary
    return tuple(selected), tuple(reasons), tuple(overlaps)


def plan_syn_modules(
    sequence: str,
    rules: SYNModuleRules,
    *,
    design_version_id: str,
    qc_result: SYNSequenceQCResult | None = None,
) -> SYNModulePlan:
    """Choose the default length-based SYN route and module boundaries."""
    normalized = normalize_dna(sequence)
    sequence_length = len(normalized)

    if sequence_length <= rules.single_pool_default_max_bp:
        route = SYNRoute.SINGLE_POOL
        sizes = (sequence_length,)
        requires_confirmation = False
        routing_reason = "length_at_or_below_1000_bp"
    elif sequence_length <= rules.single_pool_warning_max_bp:
        route = SYNRoute.SINGLE_POOL
        sizes = (sequence_length,)
        requires_confirmation = True
        routing_reason = "length_1001_to_1200_bp_single_pool_warning"
    else:
        route = SYNRoute.MODULAR
        sizes = _module_sizes(sequence_length, rules)
        requires_confirmation = len(sizes) == rules.warning_module_count
        routing_reason = (
            "more_than_5_modules_staged_assembly"
            if len(sizes) > rules.warning_module_count
            else "length_over_1200_bp_modular"
        )

    if qc_result is not None and qc_result.design_version_id != design_version_id:
        raise ValueError("QC design_version_id does not match module plan")
    boundaries, boundary_reasons, module_overlaps = _adjust_boundaries_for_qc(
        normalized,
        sizes,
        rules,
        qc_result,
    )

    modules: list[SYNModule] = []
    start = 0
    ends = (*boundaries, sequence_length)
    for ordinal, end in enumerate(ends, start=1):
        modules.append(
            SYNModule(
                design_version_id=design_version_id,
                module_id=f"module-{ordinal:02d}",
                ordinal=ordinal,
                start=start,
                end=end,
                sequence_checksum=sha256_sequence(normalized[start:end]),
                oligo_ids=(),
                boundary_reason=(
                    boundary_reasons[ordinal - 1]
                    if ordinal <= len(boundary_reasons)
                    else "length_balanced_boundary"
                ),
                left_overlap=(
                    module_overlaps[ordinal - 2]
                    if ordinal > 1
                    else None
                ),
                right_overlap=(
                    module_overlaps[ordinal - 1]
                    if ordinal <= len(module_overlaps)
                    else None
                ),
            ),
        )
        start = end

    return SYNModulePlan(
        design_version_id=design_version_id,
        route=route,
        modules=tuple(modules),
        requires_confirmation=requires_confirmation,
        routing_reason=routing_reason,
    )
