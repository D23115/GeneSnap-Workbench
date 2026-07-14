"""Generate deterministic shRNA candidates without claiming off-target screening."""

from __future__ import annotations

from decimal import Decimal

from genesnap_workbench.domain.shrna import BlastScreenStatus, ShRNACandidate

from .dna import normalize_dna


_DNA_ALPHABET = frozenset("ACGT")
_TARGET_LENGTH = 21


def _normalize_cds(sequence: str) -> str:
    return normalize_dna(sequence)


def _contains_homopolymer(sequence: str, run_length: int = 4) -> bool:
    return any(base * run_length in sequence for base in _DNA_ALPHABET)


def _score_candidate(sequence: str) -> Decimal:
    gc_percent = Decimal(sequence.count("G") + sequence.count("C")) * Decimal("100") / Decimal(
        len(sequence),
    )
    score = Decimal("100") - abs(gc_percent - Decimal("45")) * Decimal("1.5")
    if sequence.startswith(("G", "C")):
        score += Decimal("1.5")
    if sequence.endswith(("A", "T")):
        score += Decimal("1.0")
    if "TTT" in sequence:
        score -= Decimal("4.0")
    return score.quantize(Decimal("0.01"))


def generate_shrna_candidates(
    cds_sequence: str,
    *,
    max_candidates: int = 60,
    skip_start_nt: int = 25,
) -> tuple[ShRNACandidate, ...]:
    """Return ranked 21-nt candidates for later BLAST/off-target screening.

    Positions are one-based. Candidate generation deliberately marks BLAST as
    unavailable so callers cannot confuse intrinsic ranking with off-target QC.
    """

    if max_candidates <= 0:
        raise ValueError("max_candidates must be positive")
    if skip_start_nt < 0:
        raise ValueError("skip_start_nt must not be negative")

    cds = _normalize_cds(cds_sequence)
    ranked: list[tuple[Decimal, int, str]] = []
    seen: set[str] = set()
    final_start = len(cds) - _TARGET_LENGTH
    for start_index in range(skip_start_nt, final_start + 1):
        target = cds[start_index : start_index + _TARGET_LENGTH]
        if target in seen or _contains_homopolymer(target):
            continue
        gc_percent = (target.count("G") + target.count("C")) * 100 / len(target)
        if not 30 <= gc_percent <= 60:
            continue
        seen.add(target)
        ranked.append((_score_candidate(target), start_index + 1, target))

    if not ranked:
        raise ValueError("没有找到符合基础 QC 条件的 shRNA 候选序列")

    ranked.sort(key=lambda item: (-item[0], item[1], item[2]))
    return tuple(
        ShRNACandidate(
            candidate_id=f"local-{rank:03d}",
            target_sequence=target,
            start_position=start_position,
            intrinsic_score=score,
            source_rank=rank,
            blast_status=BlastScreenStatus.UNAVAILABLE,
            blast_note="本地规则生成，尚未进行 BLAST 脱靶筛查",
        )
        for rank, (score, start_position, target) in enumerate(
            ranked[:max_candidates],
            start=1,
        )
    )
