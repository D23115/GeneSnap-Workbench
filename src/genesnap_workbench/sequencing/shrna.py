"""shRNA 的测序文件读取与完整正向 oligo 判读。"""

from __future__ import annotations

from dataclasses import dataclass
from enum import Enum
from pathlib import Path
import re

from Bio import SeqIO

from genesnap_workbench.sequence_core.dna import normalize_dna


IUPAC_READ_BASES = frozenset("ACGTRYSWKMBDHVN")


class CloneJudgmentStatus(str, Enum):
    PASS = "pass"
    FAIL = "fail"
    WARNING = "warning"


@dataclass(frozen=True, slots=True)
class SequencingRead:
    path: Path
    sample_name: str
    sequence: str
    format_name: str
    phred_quality: tuple[int, ...] = ()


@dataclass(frozen=True, slots=True)
class ShRNACloneJudgment:
    clone_name: str
    target_id: str
    status: CloneJudgmentStatus
    reason: str
    read_length: int
    match_start: int | None = None
    ambiguous_base_count: int = 0
    source_files: tuple[Path, ...] = ()


@dataclass(frozen=True, slots=True)
class ShRNAFileMatchPlan:
    matches: tuple[tuple[str, tuple[Path, ...]], ...]
    unmatched_files: tuple[Path, ...]
    ambiguous_files: tuple[Path, ...]

    def files_for(self, clone_name: str) -> tuple[Path, ...]:
        for name, paths in self.matches:
            if name == clone_name:
                return paths
        return ()


def _normalize_read_sequence(raw_sequence: str) -> str:
    lines = raw_sequence.lstrip("\ufeff").splitlines()
    sequence = "".join(
        line.strip()
        for line in lines
        if line.strip() and not line.lstrip().startswith(">")
    ).upper().replace("U", "T")
    if not sequence:
        raise ValueError("测序文件中没有可读取的碱基序列")
    invalid = sorted(set(sequence) - IUPAC_READ_BASES)
    if invalid:
        raise ValueError(f"测序序列包含不支持的字符：{', '.join(invalid)}")
    return sequence


def read_sequence_file(path: Path) -> SequencingRead:
    file_path = Path(path)
    suffix = file_path.suffix.casefold()
    if suffix in {".ab1", ".abi"}:
        record = SeqIO.read(file_path, "abi")
        sequence = _normalize_read_sequence(str(record.seq))
        qualities = tuple(record.letter_annotations.get("phred_quality", ()))
        format_name = "abi"
    elif suffix in {".seq", ".txt", ".fa", ".fasta"}:
        sequence = _normalize_read_sequence(
            file_path.read_text(encoding="utf-8", errors="replace"),
        )
        qualities = ()
        format_name = "fasta" if suffix in {".fa", ".fasta"} else "seq"
    else:
        raise ValueError(f"暂不支持的测序文件格式：{file_path.suffix or '(无扩展名)'}")
    return SequencingRead(
        path=file_path,
        sample_name=file_path.stem,
        sequence=sequence,
        format_name=format_name,
        phred_quality=qualities,
    )


def judge_shrna_read(
    *,
    clone_name: str,
    target_id: str,
    read_sequence: str,
    expected_forward_oligo: str,
    source_files: tuple[Path, ...] = (),
) -> ShRNACloneJudgment:
    read = _normalize_read_sequence(read_sequence)
    expected = normalize_dna(expected_forward_oligo)
    match_start = read.find(expected)
    ambiguous_count = sum(base not in "ACGT" for base in read)
    if match_start >= 0:
        return ShRNACloneJudgment(
            clone_name=clone_name,
            target_id=target_id,
            status=CloneJudgmentStatus.PASS,
            reason="完整正向 oligo 在测序 read 中连续且完全匹配",
            read_length=len(read),
            match_start=match_start,
            ambiguous_base_count=ambiguous_count,
            source_files=source_files,
        )
    if len(read) < len(expected):
        status = CloneJudgmentStatus.WARNING
        reason = "测序读长不足，无法覆盖完整正向 oligo"
    elif ambiguous_count / len(read) > 0.05:
        status = CloneJudgmentStatus.WARNING
        reason = "测序 read 含较多不确定碱基，需人工复核"
    else:
        status = CloneJudgmentStatus.FAIL
        reason = "未找到完整、连续且完全匹配的正向 oligo"
    return ShRNACloneJudgment(
        clone_name=clone_name,
        target_id=target_id,
        status=status,
        reason=reason,
        read_length=len(read),
        match_start=None,
        ambiguous_base_count=ambiguous_count,
        source_files=source_files,
    )


def match_shrna_sequence_files(
    sequencing_root: Path,
    clone_names: tuple[str, ...],
) -> ShRNAFileMatchPlan:
    root = Path(sequencing_root)
    supported = {".seq", ".ab1", ".abi", ".txt", ".fa", ".fasta"}
    files = tuple(
        sorted(
            (
                path
                for path in root.rglob("*")
                if path.is_file() and path.suffix.casefold() in supported
            ),
            key=lambda item: str(item).casefold(),
        ),
    )
    by_clone: dict[str, list[Path]] = {name: [] for name in clone_names}
    unmatched: list[Path] = []
    ambiguous: list[Path] = []
    patterns = {
        name: re.compile(
            rf"(?<![A-Za-z0-9]){re.escape(name)}(?![A-Za-z0-9])",
            re.IGNORECASE,
        )
        for name in clone_names
    }
    for path in files:
        matched = [name for name, pattern in patterns.items() if pattern.search(path.stem)]
        if len(matched) == 1:
            by_clone[matched[0]].append(path)
        elif len(matched) > 1:
            ambiguous.append(path)
        else:
            unmatched.append(path)
    return ShRNAFileMatchPlan(
        matches=tuple((name, tuple(by_clone[name])) for name in clone_names),
        unmatched_files=tuple(unmatched),
        ambiguous_files=tuple(ambiguous),
    )
