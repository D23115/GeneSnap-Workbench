"""Sequence-level QC rules for the SYN v0 workflow."""

from __future__ import annotations

from collections import defaultdict
from dataclasses import dataclass
from decimal import Decimal
import re

from genesnap_workbench.domain.syn import SYNQCRisk, SYNSequenceQCResult

from .dna import normalize_dna, reverse_complement, sha256_sequence


PASS = "pass"
WARNING = "warning"
HIGH_RISK = "high_risk"
INFO = "info"


@dataclass(frozen=True, slots=True)
class SYNQCRules:
    rules_version: str = "syn-qc-v1"
    gc_window_size: int = 41
    minimum_repeat_length: int = 15
    restriction_sites: tuple[tuple[str, str], ...] = (
        ("EcoRV", "GATATC"),
        ("SmaI", "CCCGGG"),
        ("EcoRI", "GAATTC"),
        ("HindIII", "AAGCTT"),
    )


def classify_local_gc_percent(gc_percent: float) -> str:
    """Classify local GC percent using the SYN protocol boundaries."""
    if gc_percent < 20.0 or gc_percent > 80.0:
        return HIGH_RISK
    if gc_percent < 25.0 or gc_percent > 75.0:
        return WARNING
    return PASS


def classify_repeat_length(length: int) -> str:
    """Classify exact repeat length using the SYN protocol boundaries."""
    if length >= 20:
        return HIGH_RISK
    if length >= 15:
        return WARNING
    return PASS


def _gc_window_risks(sequence: str, rules: SYNQCRules) -> list[SYNQCRisk]:
    window_size = min(len(sequence), rules.gc_window_size)
    starts = range(0, len(sequence) - window_size + 1)
    groups: list[list[tuple[int, int, float, str]]] = []
    current_group: list[tuple[int, int, float, str]] = []
    for start in starts:
        window = sequence[start:start + window_size]
        gc_percent = (window.count("G") + window.count("C")) / window_size * 100
        severity = classify_local_gc_percent(gc_percent)
        if severity == PASS:
            if current_group:
                groups.append(current_group)
                current_group = []
            continue
        if current_group and current_group[-1][3] != severity:
            groups.append(current_group)
            current_group = []
        current_group.append((start, start + window_size, gc_percent, severity))
    if current_group:
        groups.append(current_group)

    return [_make_gc_region_risk(group) for group in groups]


def _make_gc_region_risk(
    group: list[tuple[int, int, float, str]],
) -> SYNQCRisk:
    start = group[0][0]
    end = group[-1][1]
    percentages = [item[2] for item in group]
    severity = group[0][3]
    minimum = min(percentages)
    maximum = max(percentages)
    observed = (
        f"{minimum:.2f}%"
        if minimum == maximum
        else f"{minimum:.2f}%-{maximum:.2f}%"
    )
    return SYNQCRisk(
        rule_key="local_gc",
        severity=severity,
        start=start,
        end=end,
        observed_value=observed,
        message=f"局部 GC {observed}（{start + 1}-{end} bp）",
        requires_confirmation=severity == HIGH_RISK,
    )


def _homopolymer_risks(sequence: str) -> list[SYNQCRisk]:
    risks: list[SYNQCRisk] = []
    for match in re.finditer(r"(A+)|(C+)|(G+)|(T+)", sequence):
        motif = match.group(0)
        base = motif[0]
        length = len(motif)
        if base in "AT" and length >= 10:
            severity = HIGH_RISK
        elif base in "AT" and length >= 8:
            severity = WARNING
        elif base in "GC" and length >= 6:
            severity = HIGH_RISK
        else:
            continue
        risks.append(
            SYNQCRisk(
                rule_key="homopolymer",
                severity=severity,
                start=match.start(),
                end=match.end(),
                observed_value=f"{base} x {length}",
                message=(
                    f"连续 {base} {length} bp（{match.start() + 1}-{match.end()} bp）"
                ),
                requires_confirmation=severity == HIGH_RISK,
            ),
        )
    return risks


def _select_non_overlapping_positions(
    positions: list[int],
    length: int,
) -> tuple[int, ...]:
    selected: list[int] = []
    next_available = 0
    for position in sorted(positions):
        if position < next_available:
            continue
        selected.append(position)
        next_available = position + length
    return tuple(selected)


def _all_repeat_positions_are_covered(
    positions: tuple[int, ...],
    length: int,
    covered_intervals: list[tuple[int, int]],
) -> bool:
    return all(
        any(
            position >= covered_start
            and position + length <= covered_end
            for covered_start, covered_end in covered_intervals
        )
        for position in positions
    )


def _repeat_risks(sequence: str, rules: SYNQCRules) -> list[SYNQCRisk]:
    risks: list[SYNQCRisk] = []
    covered_direct: list[tuple[int, int]] = []
    covered_inverted: list[tuple[int, int]] = []
    maximum_length = min(len(sequence) // 2, 100)

    for length in range(maximum_length, rules.minimum_repeat_length - 1, -1):
        positions_by_motif: dict[str, list[int]] = defaultdict(list)
        for start in range(0, len(sequence) - length + 1):
            positions_by_motif[sequence[start:start + length]].append(start)

        for motif, positions in positions_by_motif.items():
            selected = _select_non_overlapping_positions(positions, length)
            if len(selected) < 2 or _all_repeat_positions_are_covered(
                selected,
                length,
                covered_direct,
            ):
                continue
            risks.append(_make_repeat_risk("direct", motif, selected, length))
            covered_direct.extend(
                (position, position + length) for position in selected
            )

        for motif, positions in positions_by_motif.items():
            reverse = reverse_complement(motif)
            if motif >= reverse or reverse not in positions_by_motif:
                continue
            tagged_positions = sorted(
                [(position, "forward") for position in positions]
                + [
                    (position, "reverse")
                    for position in positions_by_motif[reverse]
                ],
            )
            selected_tagged: list[tuple[int, str]] = []
            next_available = 0
            for position, orientation in tagged_positions:
                if position < next_available:
                    continue
                selected_tagged.append((position, orientation))
                next_available = position + length
            orientations = {orientation for _, orientation in selected_tagged}
            selected = tuple(position for position, _ in selected_tagged)
            if (
                len(selected) < 2
                or orientations != {"forward", "reverse"}
                or _all_repeat_positions_are_covered(
                    selected,
                    length,
                    covered_inverted,
                )
            ):
                continue
            risks.append(_make_repeat_risk("inverted", motif, selected, length))
            covered_inverted.extend(
                (position, position + length) for position in selected
            )
    return risks


def _make_repeat_risk(
    orientation: str,
    motif: str,
    positions: tuple[int, ...],
    length: int,
) -> SYNQCRisk:
    severity = classify_repeat_length(length)
    position_text = ",".join(str(position + 1) for position in positions)
    range_text = "、".join(
        f"{position + 1}-{position + length}" for position in positions
    )
    return SYNQCRisk(
        rule_key="repeat",
        severity=severity,
        start=min(positions),
        end=max(positions) + length,
        observed_value=(
            f"{orientation}:{motif}; length={length}; "
            f"count={len(positions)}; positions={position_text}"
        ),
        message=f"{orientation} repeat {length} bp，位置 {range_text}",
        requires_confirmation=severity == HIGH_RISK,
    )


def _restriction_site_risks(
    sequence: str,
    rules: SYNQCRules,
) -> list[SYNQCRisk]:
    risks: list[SYNQCRisk] = []
    for name, site in rules.restriction_sites:
        normalized_site = normalize_dna(site)
        search_sites = {normalized_site, reverse_complement(normalized_site)}
        for search_site in search_sites:
            start = sequence.find(search_site)
            while start >= 0:
                risks.append(
                    SYNQCRisk(
                        rule_key="restriction_site",
                        severity=INFO,
                        start=start,
                        end=start + len(search_site),
                        observed_value=f"{name}:{search_site}",
                        message=(
                            f"检测到内部 {name} 位点（{start + 1}-{start + len(search_site)} bp）"
                        ),
                        requires_confirmation=False,
                    ),
                )
                start = sequence.find(search_site, start + 1)
    return risks


def evaluate_syn_sequence(
    sequence: str,
    rules: SYNQCRules,
    *,
    design_version_id: str,
) -> SYNSequenceQCResult:
    """Evaluate a strict DNA sequence without modifying its bases."""
    normalized = normalize_dna(sequence)
    risks = [
        *_gc_window_risks(normalized, rules),
        *_homopolymer_risks(normalized),
        *_repeat_risks(normalized, rules),
        *_restriction_site_risks(normalized, rules),
    ]
    risks.sort(
        key=lambda risk: (
            risk.start,
            risk.end,
            risk.rule_key,
            risk.observed_value,
        ),
    )
    return SYNSequenceQCResult(
        design_version_id=design_version_id,
        rules_version=rules.rules_version,
        sequence_checksum=sha256_sequence(normalized),
        sequence_length=len(normalized),
        overall_gc_percent=(
            Decimal((normalized.count("G") + normalized.count("C")) * 100)
            / Decimal(len(normalized))
        ),
        risks=tuple(risks),
        blocked_reasons=(),
        confirmable_warnings=tuple(
            risk.message for risk in risks if risk.requires_confirmation
        ),
    )
