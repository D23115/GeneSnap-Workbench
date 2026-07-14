"""表达类 insert 的整质粒/长读长测序初筛。"""

from __future__ import annotations

from dataclasses import dataclass
from enum import Enum

from Bio import Align
from Bio.Seq import Seq


class ExpressionCloneJudgmentStatus(str, Enum):
    PASS = "pass"
    FAIL = "fail"
    WARNING = "warning"


@dataclass(frozen=True, slots=True)
class ExpressionReadJudgment:
    clone_name: str
    construct_id: str
    status: ExpressionCloneJudgmentStatus
    reason: str
    orientation: str
    coverage: float
    identity: float
    substitution_count: int
    insertion_count: int
    deletion_count: int
    frameshift: bool
    premature_stop: bool


def _normalize_read(sequence: str) -> str:
    normalized = "".join(sequence.split()).upper().replace("U", "T")
    if not normalized:
        raise ValueError("测序序列为空")
    invalid = sorted(set(normalized) - set("ACGTN"))
    if invalid:
        raise ValueError(f"测序序列含不支持字符：{', '.join(invalid)}")
    return normalized


def _alignment_metrics(expected: str, observed: str):
    aligner = Align.PairwiseAligner(mode="local")
    aligner.match_score = 2
    aligner.mismatch_score = -1
    aligner.open_gap_score = -5
    aligner.extend_gap_score = -1
    alignment = aligner.align(expected, observed)[0]
    expected_blocks, observed_blocks = alignment.aligned
    substitutions = 0
    aligned_bases = 0
    matches = 0
    for (expected_start, expected_end), (observed_start, observed_end) in zip(
        expected_blocks,
        observed_blocks,
        strict=True,
    ):
        expected_block = expected[expected_start:expected_end]
        observed_block = observed[observed_start:observed_end]
        block_matches = sum(
            left == right
            for left, right in zip(expected_block, observed_block, strict=True)
        )
        matches += block_matches
        substitutions += len(expected_block) - block_matches
        aligned_bases += len(expected_block)

    insertion_events: list[int] = []
    deletion_events: list[int] = []
    for index in range(1, len(expected_blocks)):
        expected_gap = int(expected_blocks[index][0] - expected_blocks[index - 1][1])
        observed_gap = int(observed_blocks[index][0] - observed_blocks[index - 1][1])
        if expected_gap:
            deletion_events.append(expected_gap)
        if observed_gap:
            insertion_events.append(observed_gap)
    coverage = aligned_bases / len(expected)
    identity = matches / aligned_bases if aligned_bases else 0.0
    return (
        alignment.score,
        coverage,
        identity,
        substitutions,
        insertion_events,
        deletion_events,
        observed_blocks,
    )


def _premature_stop(
    expected_coding_sequence: str,
    observed: str,
    observed_blocks,
    coverage: float,
) -> bool:
    if coverage < 0.95 or not len(observed_blocks):
        return False
    start = int(observed_blocks[0][0])
    end = int(observed_blocks[-1][1])
    observed_region = observed[start:end]
    expected_has_stop = expected_coding_sequence[-3:] in {"TAA", "TAG", "TGA"}
    if len(observed_region) % 3:
        return False
    protein = str(Seq(observed_region).translate())
    allowed_terminal = expected_has_stop and protein.endswith("*")
    return "*" in (protein[:-1] if allowed_terminal else protein)


def judge_expression_read(
    *,
    clone_name: str,
    construct_id: str,
    read_sequence: str,
    expected_insert_sequence: str,
    expected_coding_sequence: str,
) -> ExpressionReadJudgment:
    read = _normalize_read(read_sequence)
    expected = "".join(expected_insert_sequence.split()).upper()
    coding = "".join(expected_coding_sequence.split()).upper()
    reverse_read = str(Seq(read).reverse_complement())
    for orientation, oriented_read in (("forward", read), ("reverse", reverse_read)):
        if expected in oriented_read:
            return ExpressionReadJudgment(
                clone_name=clone_name,
                construct_id=construct_id,
                status=ExpressionCloneJudgmentStatus.PASS,
                reason="预期 insert 完整且完全一致",
                orientation=orientation,
                coverage=1.0,
                identity=1.0,
                substitution_count=0,
                insertion_count=0,
                deletion_count=0,
                frameshift=False,
                premature_stop=False,
            )

    candidates = []
    for orientation, oriented_read in (("forward", read), ("reverse", reverse_read)):
        metrics = _alignment_metrics(expected, oriented_read)
        candidates.append((metrics[0], metrics[1], metrics[2], orientation, oriented_read, metrics))
    _, coverage, identity, orientation, oriented_read, metrics = max(
        candidates,
        key=lambda item: (item[0], item[1], item[2]),
    )
    (
        _,
        coverage,
        identity,
        substitutions,
        insertion_events,
        deletion_events,
        observed_blocks,
    ) = metrics
    insertion_count = sum(insertion_events)
    deletion_count = sum(deletion_events)
    frameshift = any(length % 3 for length in insertion_events + deletion_events)
    premature_stop = _premature_stop(coding, oriented_read, observed_blocks, coverage)
    details = [f"覆盖 {coverage:.1%}", f"一致性 {identity:.1%}"]
    if substitutions:
        details.append(f"{substitutions} 个碱基替换")
    if insertion_count:
        details.append(f"插入 {insertion_count} bp")
    if deletion_count:
        details.append(f"缺失 {deletion_count} bp")
    if frameshift:
        details.append("可能导致移码")
    if premature_stop:
        details.append("检测到提前终止风险")
    return ExpressionReadJudgment(
        clone_name=clone_name,
        construct_id=construct_id,
        status=ExpressionCloneJudgmentStatus.WARNING,
        reason="；".join(details) + "；需人工复核",
        orientation=orientation,
        coverage=coverage,
        identity=identity,
        substitution_count=substitutions,
        insertion_count=insertion_count,
        deletion_count=deletion_count,
        frameshift=frameshift,
        premature_stop=premature_stop,
    )
