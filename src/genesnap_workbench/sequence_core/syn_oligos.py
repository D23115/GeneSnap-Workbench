"""Deterministic overlapping-oligo design for SYN modules."""

from __future__ import annotations

from dataclasses import dataclass, replace
from decimal import Decimal
from typing import Protocol

from Bio import __version__ as biopython_version
from Bio.SeqUtils import MeltingTemp

from genesnap_workbench.domain.syn import (
    SYNAssemblyOligo,
    SYNModule,
    SYNModulePlan,
    SYNThermodynamicMetadata,
)

from .dna import normalize_dna, reverse_complement


class ThermodynamicAnalyzer(Protocol):
    def analyze(self, sequence: str) -> SYNThermodynamicMetadata: ...


@dataclass(frozen=True, slots=True)
class BiopythonThermodynamicAnalyzer:
    sodium_mM: float = 50.0
    magnesium_mM: float = 1.5
    strand_concentration_nM: float = 250.0

    def analyze(self, sequence: str) -> SYNThermodynamicMetadata:
        normalized = normalize_dna(sequence)
        tm = MeltingTemp.Tm_NN(
            normalized,
            nn_table=MeltingTemp.DNA_NN4,
            Na=self.sodium_mM,
            Mg=self.magnesium_mM,
            dnac1=self.strand_concentration_nM,
            dnac2=self.strand_concentration_nM,
        )
        return SYNThermodynamicMetadata(
            analyzer_name="Bio.SeqUtils.MeltingTemp.Tm_NN",
            analyzer_version=biopython_version,
            tm_celsius=Decimal(str(round(tm, 4))),
            parameters=(
                ("nn_table", "DNA_NN4"),
                ("sodium_mM", str(self.sodium_mM)),
                ("magnesium_mM", str(self.magnesium_mM)),
                (
                    "strand_concentration_nM",
                    str(self.strand_concentration_nM),
                ),
            ),
        )


@dataclass(frozen=True, slots=True)
class SYNOligoRules:
    preferred_oligo_min: int = 58
    preferred_oligo_max: int = 60
    max_oligo_length: int = 65
    target_overlap_length: int = 20
    min_overlap_length: int = 18
    max_overlap_length: int = 25
    target_overlap_tm_celsius: float = 60.0
    normal_tm_min_celsius: float = 58.0
    normal_tm_max_celsius: float = 65.0
    high_tm_max_celsius: float = 70.0
    preferred_tm_spread_celsius: float = 4.0
    high_risk_tm_spread_celsius: float = 6.0

    def __post_init__(self) -> None:
        if self.max_oligo_length > 65:
            raise ValueError("Assembly oligo hard limit is 65 nt")
        if self.min_overlap_length < 18 or self.max_overlap_length > 25:
            raise ValueError("SYN overlap length must stay within 18-25 bp")
        if not (
            self.min_overlap_length
            <= self.target_overlap_length
            <= self.max_overlap_length
        ):
            raise ValueError("Target overlap length must be within configured range")


@dataclass(frozen=True, slots=True)
class SYNOverlap:
    overlap_id: str
    module_id: str
    left_oligo_id: str
    right_oligo_id: str
    start: int
    end: int
    sequence: str
    occurrence_count: int
    thermodynamic_metadata: SYNThermodynamicMetadata


@dataclass(frozen=True, slots=True)
class SYNOligoDesignResult:
    design_version_id: str
    module_plan: SYNModulePlan
    oligos: tuple[SYNAssemblyOligo, ...]
    overlaps: tuple[SYNOverlap, ...]
    requires_confirmation: bool
    warnings: tuple[str, ...]


class SYNOligoDesignFailure(ValueError):
    def __init__(
        self,
        reasons: tuple[str, ...],
        *,
        failure_coordinate: int | None = None,
    ) -> None:
        self.reasons = reasons
        self.failure_coordinate = failure_coordinate
        super().__init__("; ".join(reasons))


@dataclass(frozen=True, slots=True)
class _Candidate:
    intervals: tuple[tuple[int, int], ...]
    overlap_intervals: tuple[tuple[int, int], ...]
    overlap_metadata: tuple[SYNThermodynamicMetadata, ...]
    score: tuple[float | int, ...]
    warnings: tuple[str, ...]
    requires_confirmation: bool


def _overlap_order(rules: SYNOligoRules) -> tuple[int, ...]:
    return tuple(
        sorted(
            range(rules.min_overlap_length, rules.max_overlap_length + 1),
            key=lambda value: (
                abs(value - rules.target_overlap_length),
                value,
            ),
        ),
    )


def _count_occurrences(sequence: str, motif: str) -> int:
    count = 0
    start = sequence.find(motif)
    while start >= 0:
        count += 1
        start = sequence.find(motif, start + 1)
    return count


def _distributed_lengths(total_bases: int, count: int) -> tuple[int, ...]:
    base, remainder = divmod(total_bases, count)
    return tuple(base + (1 if index < remainder else 0) for index in range(count))


def _candidate_intervals(
    module_length: int,
    oligo_count: int,
    overlap_length: int,
    rules: SYNOligoRules,
) -> tuple[tuple[int, int], ...] | None:
    total_oligo_bases = module_length + (oligo_count - 1) * overlap_length
    lengths = _distributed_lengths(total_oligo_bases, oligo_count)
    if (
        min(lengths) < rules.preferred_oligo_min
        or max(lengths) > rules.max_oligo_length
    ):
        return None

    intervals: list[tuple[int, int]] = []
    start = 0
    for length in lengths:
        end = start + length
        intervals.append((start, end))
        start = end - overlap_length
    if intervals[-1][1] != module_length:
        return None
    return tuple(intervals)


def _evaluate_candidate(
    full_sequence: str,
    module_sequence: str,
    intervals: tuple[tuple[int, int], ...],
    rules: SYNOligoRules,
    thermodynamics: ThermodynamicAnalyzer,
) -> _Candidate | None:
    overlap_intervals = tuple(
        (right_start, left_end)
        for (_, left_end), (right_start, _) in zip(intervals, intervals[1:])
    )
    overlap_sequences = tuple(
        module_sequence[start:end] for start, end in overlap_intervals
    )
    if any(_count_occurrences(full_sequence, item) != 1 for item in overlap_sequences):
        return None

    metadata = tuple(thermodynamics.analyze(item) for item in overlap_sequences)
    tms = tuple(float(item.tm_celsius) for item in metadata)
    spread = max(tms) - min(tms) if tms else 0.0
    high_risk_count = sum(
        tm < rules.normal_tm_min_celsius or tm > rules.high_tm_max_celsius
        for tm in tms
    )
    warning_count = sum(
        rules.normal_tm_max_celsius < tm <= rules.high_tm_max_celsius
        for tm in tms
    )
    spread_high_risk = spread > rules.high_risk_tm_spread_celsius
    spread_warning = spread > rules.preferred_tm_spread_celsius
    lengths = tuple(end - start for start, end in intervals)
    long_oligo_count = sum(length > rules.preferred_oligo_max for length in lengths)

    warnings: list[str] = []
    if high_risk_count:
        warnings.append(f"{high_risk_count} 个 overlap Tm 超出 58-70°C")
    if warning_count:
        warnings.append(f"{warning_count} 个 overlap Tm 位于 65-70°C")
    if spread_warning:
        warnings.append(f"overlap Tm 极差为 {spread:.2f}°C")
    if long_oligo_count:
        warnings.append(f"{long_oligo_count} 条 oligo 使用 61-65 nt 软范围")

    score = (
        high_risk_count + int(spread_high_risk),
        long_oligo_count,
        warning_count + int(spread_warning),
        round(spread, 4),
        round(
            sum(abs(tm - rules.target_overlap_tm_celsius) for tm in tms),
            4,
        ),
        abs(len(overlap_sequences[0]) - rules.target_overlap_length)
        if overlap_sequences
        else 0,
        sum(abs(length - 59) for length in lengths),
        len(intervals),
    )
    return _Candidate(
        intervals=intervals,
        overlap_intervals=overlap_intervals,
        overlap_metadata=metadata,
        score=score,
        warnings=tuple(warnings),
        requires_confirmation=(
            high_risk_count > 0 or spread_high_risk or long_oligo_count > 0
        ),
    )


def _design_module_candidate(
    full_sequence: str,
    module: SYNModule,
    rules: SYNOligoRules,
    thermodynamics: ThermodynamicAnalyzer,
) -> _Candidate:
    module_sequence = full_sequence[module.start:module.end]
    module_length = len(module_sequence)
    candidates: list[_Candidate] = []
    saw_length_feasible = False

    maximum_oligo_count = max(
        1,
        module_length // max(1, rules.preferred_oligo_min - rules.max_overlap_length)
        + 3,
    )
    for overlap_length in _overlap_order(rules):
        for oligo_count in range(1, maximum_oligo_count + 1):
            intervals = _candidate_intervals(
                module_length,
                oligo_count,
                overlap_length,
                rules,
            )
            if intervals is None:
                continue
            saw_length_feasible = True
            candidate = _evaluate_candidate(
                full_sequence,
                module_sequence,
                intervals,
                rules,
                thermodynamics,
            )
            if candidate is not None:
                candidates.append(candidate)

    if not candidates:
        reasons = (
            ("18-25 bp 范围内的候选 overlap 不唯一",)
            if saw_length_feasible
            else ("65 nt assembly oligo 硬上限内不存在可行长度方案",)
        )
        raise SYNOligoDesignFailure(reasons, failure_coordinate=module.start)
    return min(candidates, key=lambda item: item.score)


def design_assembly_oligos(
    sequence: str,
    module_plan: SYNModulePlan,
    rules: SYNOligoRules,
    thermodynamics: ThermodynamicAnalyzer,
    *,
    design_version_id: str,
    project_id: str,
    target_name: str,
) -> SYNOligoDesignResult:
    """Design all module oligos and return an updated module plan."""
    normalized = normalize_dna(sequence)
    if module_plan.design_version_id != design_version_id:
        raise ValueError("Module plan design_version_id does not match design")

    oligos: list[SYNAssemblyOligo] = []
    overlaps: list[SYNOverlap] = []
    updated_modules: list[SYNModule] = []
    warnings: list[str] = []
    requires_confirmation = module_plan.requires_confirmation

    for module in module_plan.modules:
        candidate = _design_module_candidate(
            normalized,
            module,
            rules,
            thermodynamics,
        )
        module_oligo_ids: list[str] = []
        local_oligos: list[SYNAssemblyOligo] = []
        for local_index, (local_start, local_end) in enumerate(
            candidate.intervals,
            start=1,
        ):
            global_index = len(oligos) + 1
            oligo_id = f"{project_id}-{target_name}-ASM-{global_index:02d}"
            module_oligo_ids.append(oligo_id)
            global_start = module.start + local_start
            global_end = module.start + local_end
            target_fragment = normalized[global_start:global_end]
            strand = "forward" if local_index % 2 else "reverse"
            oligo_sequence = (
                target_fragment
                if strand == "forward"
                else reverse_complement(target_fragment)
            )
            left_overlap = (
                (
                    module.start + candidate.overlap_intervals[local_index - 2][0],
                    module.start + candidate.overlap_intervals[local_index - 2][1],
                )
                if local_index > 1
                else None
            )
            right_overlap = (
                (
                    module.start + candidate.overlap_intervals[local_index - 1][0],
                    module.start + candidate.overlap_intervals[local_index - 1][1],
                )
                if local_index <= len(candidate.overlap_intervals)
                else None
            )
            metadata_index = min(
                max(local_index - 2, 0),
                max(len(candidate.overlap_metadata) - 1, 0),
            )
            metadata = (
                candidate.overlap_metadata[metadata_index]
                if candidate.overlap_metadata
                else thermodynamics.analyze(target_fragment)
            )
            local_oligos.append(
                SYNAssemblyOligo(
                    oligo_id=oligo_id,
                    design_version_id=design_version_id,
                    name=f"{target_name}-A{global_index:02d}",
                    sequence=oligo_sequence,
                    strand=strand,
                    start=global_start,
                    end=global_end,
                    pool_id=f"pool-{module.ordinal:02d}",
                    module_id=module.module_id,
                    overlap_left=left_overlap,
                    overlap_right=right_overlap,
                    tm_metadata=metadata,
                ),
            )

        for overlap_index, ((local_start, local_end), metadata) in enumerate(
            zip(candidate.overlap_intervals, candidate.overlap_metadata),
            start=1,
        ):
            global_start = module.start + local_start
            global_end = module.start + local_end
            overlap_sequence = normalized[global_start:global_end]
            overlaps.append(
                SYNOverlap(
                    overlap_id=f"{module.module_id}-overlap-{overlap_index:02d}",
                    module_id=module.module_id,
                    left_oligo_id=module_oligo_ids[overlap_index - 1],
                    right_oligo_id=module_oligo_ids[overlap_index],
                    start=global_start,
                    end=global_end,
                    sequence=overlap_sequence,
                    occurrence_count=_count_occurrences(normalized, overlap_sequence),
                    thermodynamic_metadata=metadata,
                ),
            )
        oligos.extend(local_oligos)
        updated_modules.append(replace(module, oligo_ids=tuple(module_oligo_ids)))
        warnings.extend(candidate.warnings)
        requires_confirmation = requires_confirmation or candidate.requires_confirmation

    updated_plan = replace(module_plan, modules=tuple(updated_modules))
    reconstructed = reconstruct_from_assembly_oligos(tuple(oligos))
    if reconstructed != normalized:
        raise SYNOligoDesignFailure(
            ("设计的 assembly oligo 无法完整重建目标序列",),
        )
    return SYNOligoDesignResult(
        design_version_id=design_version_id,
        module_plan=updated_plan,
        oligos=tuple(oligos),
        overlaps=tuple(overlaps),
        requires_confirmation=requires_confirmation,
        warnings=tuple(dict.fromkeys(warnings)),
    )


def reconstruct_from_assembly_oligos(
    oligos: tuple[SYNAssemblyOligo, ...],
) -> str:
    """Reconstruct target-oriented DNA from ordered alternating oligos."""
    if not oligos:
        raise ValueError("At least one assembly oligo is required")
    ordered = sorted(oligos, key=lambda item: (item.start, item.end))
    first = ordered[0]
    if first.start != 0:
        raise ValueError("Assembly oligos must start at coordinate 0")

    assembled = ""
    for oligo in ordered:
        fragment = (
            normalize_dna(oligo.sequence)
            if oligo.strand == "forward"
            else reverse_complement(oligo.sequence)
        )
        if len(fragment) != oligo.end - oligo.start:
            raise ValueError(f"Oligo coordinates do not match sequence: {oligo.oligo_id}")
        if oligo.start > len(assembled):
            raise ValueError(f"Gap before assembly oligo: {oligo.oligo_id}")
        overlap_length = len(assembled) - oligo.start
        if overlap_length and assembled[oligo.start:] != fragment[:overlap_length]:
            raise ValueError(f"Overlap mismatch at assembly oligo: {oligo.oligo_id}")
        assembled += fragment[overlap_length:]
    return assembled
