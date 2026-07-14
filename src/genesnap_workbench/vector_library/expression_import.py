"""Resolve expression-vector insertion boundaries from sites or homology arms."""

from __future__ import annotations

from dataclasses import dataclass

from genesnap_workbench.sequence_core.dna import normalize_dna, reverse_complement

from .models import RestrictionSite


COMMON_EXPRESSION_RESTRICTION_SITES: tuple[RestrictionSite, ...] = (
    RestrictionSite("NheI", "GCTAGC", 1),
    RestrictionSite("BamHI", "GGATCC", 1),
    RestrictionSite("EcoRI", "GAATTC", 1),
    RestrictionSite("HindIII", "AAGCTT", 1),
    RestrictionSite("XhoI", "CTCGAG", 1),
    RestrictionSite("AgeI", "ACCGGT", 1),
    RestrictionSite("NotI", "GCGGCCGC", 2),
    RestrictionSite("KpnI", "GGTACC", 5),
    RestrictionSite("BglII", "AGATCT", 1),
    RestrictionSite("MluI", "ACGCGT", 1),
    RestrictionSite("SmaI", "CCCGGG", 3),
    RestrictionSite("EcoRV", "GATATC", 3),
)


@dataclass(frozen=True, slots=True)
class RestrictionSiteOccurrence:
    """One recognition-site match in a normalized, zero-indexed vector sequence."""

    enzyme_name: str
    recognition_sequence: str
    cut_offset: int
    start: int
    end: int
    cut_position: int


@dataclass(frozen=True, slots=True)
class ExpressionInsertionResolution:
    """Resolved expression insertion settings and UI-facing provenance metadata."""

    left_boundary: int
    right_boundary: int
    left_primer_homology: str
    right_primer_homology: str
    source: str
    left_site: RestrictionSiteOccurrence | None
    right_site: RestrictionSiteOccurrence | None
    left_homology_length: int
    right_homology_length: int
    left_homology_gc_percent: float
    right_homology_gc_percent: float
    left_homology_tm: int
    right_homology_tm: int


def scan_restriction_sites(
    vector_sequence: str,
) -> tuple[RestrictionSiteOccurrence, ...]:
    """Return every common expression-site occurrence in circular-vector order."""
    sequence = normalize_dna(vector_sequence)
    sequence_length = len(sequence)
    max_site_length = max(len(site.sequence) for site in COMMON_EXPRESSION_RESTRICTION_SITES)
    wrapped_sequence = sequence + (
        sequence * ((max_site_length - 1 + sequence_length - 1) // sequence_length)
    )[: max_site_length - 1]
    occurrences: list[RestrictionSiteOccurrence] = []
    for site in COMMON_EXPRESSION_RESTRICTION_SITES:
        start = 0
        while True:
            position = wrapped_sequence.find(site.sequence, start)
            if position < 0 or position >= sequence_length:
                break
            occurrences.append(
                RestrictionSiteOccurrence(
                    enzyme_name=site.name,
                    recognition_sequence=site.sequence,
                    cut_offset=site.cut_offset,
                    start=position,
                    end=(position + len(site.sequence)) % sequence_length,
                    cut_position=(position + site.cut_offset) % sequence_length,
                ),
            )
            start = position + 1
    return tuple(
        sorted(
            occurrences,
            key=lambda occurrence: (
                occurrence.start,
                occurrence.enzyme_name,
                occurrence.recognition_sequence,
            ),
        ),
    )


def scan_expression_restriction_sites(
    vector_sequence: str,
) -> tuple[RestrictionSiteOccurrence, ...]:
    """Compatibility alias for :func:`scan_restriction_sites`."""
    return scan_restriction_sites(vector_sequence)


def resolve_restriction_insertion(
    vector_sequence: str,
    left_occurrence: RestrictionSiteOccurrence,
    right_occurrence: RestrictionSiteOccurrence,
) -> ExpressionInsertionResolution:
    """Resolve explicitly selected restriction-site occurrences into insertion flanks."""
    sequence = normalize_dna(vector_sequence)
    _validate_occurrence(sequence, left_occurrence)
    _validate_occurrence(sequence, right_occurrence)
    if _spans_circular_origin(sequence, left_occurrence) or _spans_circular_origin(
        sequence,
        right_occurrence,
    ):
        raise ValueError("所选酶切位点跨越载体首尾，当前线性插入区间模型无法解析")
    left_boundary = left_occurrence.cut_position
    right_boundary = right_occurrence.cut_position
    _validate_boundaries(sequence, left_boundary, right_boundary)
    left_arm = _select_left_homology(sequence, left_boundary)
    right_top_arm = _select_right_homology(sequence, right_boundary)
    return _resolution(
        left_boundary=left_boundary,
        right_boundary=right_boundary,
        left_primer_homology=left_arm,
        right_primer_homology=reverse_complement(right_top_arm),
        source="restriction_sites",
        left_site=left_occurrence,
        right_site=right_occurrence,
    )


def resolve_manual_homology(
    vector_sequence: str,
    forward_primer_homology: str,
    reverse_primer_homology: str,
) -> ExpressionInsertionResolution:
    """Infer a unique, ordered insertion interval from ordered-primer homology arms."""
    sequence = normalize_dna(vector_sequence)
    forward_homology = normalize_dna(forward_primer_homology)
    reverse_homology = normalize_dna(reverse_primer_homology)
    right_top_homology = reverse_complement(reverse_homology)
    left_boundaries = tuple(
        position + len(forward_homology)
        for position in _match_positions(sequence, forward_homology)
    )
    right_boundaries = _match_positions(sequence, right_top_homology)
    if not left_boundaries:
        raise ValueError("F 引物 5' 同源臂未在载体序列中找到精确匹配")
    if not right_boundaries:
        raise ValueError("R 引物 5' 同源臂未在载体序列中找到精确匹配")
    if len(left_boundaries) > 1:
        raise ValueError("F 引物 5' 同源臂存在多个匹配位置，无法唯一确定插入边界")
    if len(right_boundaries) > 1:
        raise ValueError("R 引物 5' 同源臂存在多个匹配位置，无法唯一确定插入边界")

    combinations = tuple(
        (left_boundary, right_boundary)
        for left_boundary in left_boundaries
        for right_boundary in right_boundaries
        if left_boundary <= right_boundary
    )
    if not combinations:
        raise ValueError("F/R 引物 5' 同源臂无法形成左边界不大于右边界的插入区间")

    left_boundary, right_boundary = combinations[0]
    return _resolution(
        left_boundary=left_boundary,
        right_boundary=right_boundary,
        left_primer_homology=forward_homology,
        right_primer_homology=reverse_homology,
        source="manual_homology",
        left_site=None,
        right_site=None,
    )


def _match_positions(sequence: str, motif: str) -> tuple[int, ...]:
    positions: list[int] = []
    start = 0
    while True:
        position = sequence.find(motif, start)
        if position < 0:
            return tuple(positions)
        positions.append(position)
        start = position + 1


def _validate_boundaries(sequence: str, left_boundary: int, right_boundary: int) -> None:
    if left_boundary < 0 or right_boundary > len(sequence):
        raise ValueError("酶切位点切点超出载体序列范围")
    if left_boundary > right_boundary:
        raise ValueError("左侧酶切位点切点不能位于右侧酶切位点之后")


def _validate_occurrence(
    sequence: str,
    occurrence: RestrictionSiteOccurrence,
) -> None:
    if not 0 <= occurrence.start < len(sequence):
        raise ValueError("所选酶切位点起始位置超出载体序列范围")
    if occurrence.recognition_sequence != _circular_slice(
        sequence,
        occurrence.start,
        len(occurrence.recognition_sequence),
    ):
        raise ValueError("所选酶切位点与当前载体序列不一致")
    expected_cut = (occurrence.start + occurrence.cut_offset) % len(sequence)
    if occurrence.cut_position != expected_cut:
        raise ValueError("所选酶切位点切点信息无效")


def _spans_circular_origin(
    sequence: str,
    occurrence: RestrictionSiteOccurrence,
) -> bool:
    return occurrence.start + len(occurrence.recognition_sequence) > len(sequence)


def _circular_slice(sequence: str, start: int, length: int) -> str:
    return "".join(sequence[(start + offset) % len(sequence)] for offset in range(length))


def _select_left_homology(sequence: str, boundary: int) -> str:
    if boundary < 18:
        raise ValueError("左侧插入边界前不足 18 bp，无法生成同源臂")
    candidates = (
        sequence[boundary - length : boundary]
        for length in range(18, min(25, boundary) + 1)
    )
    return min(candidates, key=_homology_score)


def _select_right_homology(sequence: str, boundary: int) -> str:
    available = len(sequence) - boundary
    if available < 18:
        raise ValueError("右侧插入边界后不足 18 bp，无法生成同源臂")
    candidates = (
        sequence[boundary : boundary + length]
        for length in range(18, min(25, available) + 1)
    )
    return min(candidates, key=_homology_score)


def _homology_score(sequence: str) -> tuple[int, int, int, float, str]:
    gc_percent = _gc_percent(sequence)
    return (
        0 if 40.0 <= gc_percent <= 60.0 else 1,
        abs(_wallace_tm(sequence) - 60),
        abs(len(sequence) - 20),
        min(abs(gc_percent - 40.0), abs(gc_percent - 60.0)),
        sequence,
    )


def _resolution(
    *,
    left_boundary: int,
    right_boundary: int,
    left_primer_homology: str,
    right_primer_homology: str,
    source: str,
    left_site: RestrictionSiteOccurrence | None,
    right_site: RestrictionSiteOccurrence | None,
) -> ExpressionInsertionResolution:
    return ExpressionInsertionResolution(
        left_boundary=left_boundary,
        right_boundary=right_boundary,
        left_primer_homology=left_primer_homology,
        right_primer_homology=right_primer_homology,
        source=source,
        left_site=left_site,
        right_site=right_site,
        left_homology_length=len(left_primer_homology),
        right_homology_length=len(right_primer_homology),
        left_homology_gc_percent=_gc_percent(left_primer_homology),
        right_homology_gc_percent=_gc_percent(right_primer_homology),
        left_homology_tm=_wallace_tm(left_primer_homology),
        right_homology_tm=_wallace_tm(right_primer_homology),
    )


def _gc_percent(sequence: str) -> float:
    return (sequence.count("G") + sequence.count("C")) / len(sequence) * 100


def _wallace_tm(sequence: str) -> int:
    return 2 * (sequence.count("A") + sequence.count("T")) + 4 * (
        sequence.count("G") + sequence.count("C")
    )
