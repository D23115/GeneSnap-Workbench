"""Sequence-level comparison helpers for circular plasmid vectors."""

from __future__ import annotations

from collections import Counter
from dataclasses import dataclass
from heapq import nsmallest
import hashlib
from io import StringIO
from pathlib import Path
import re

from Bio import SeqIO
from Bio.Seq import Seq
from rapidfuzz.distance import Levenshtein


IUPAC_DNA = frozenset("ACGTRYSWKMBDHVN")


@dataclass(frozen=True)
class DifferenceBlock:
    kind: str
    internal_start: int
    internal_end: int
    public_start: int
    public_end: int
    internal_sequence: str
    public_sequence: str


@dataclass(frozen=True)
class CircularComparison:
    edit_distance: int
    identity_percent: float
    substitutions: int
    internal_extra_bp: int
    internal_missing_bp: int
    difference_block_count: int
    difference_blocks: tuple[DifferenceBlock, ...]
    difference_blocks_truncated: bool
    public_orientation: str
    public_rotation: int


@dataclass(frozen=True)
class ParsedSequence:
    path: Path
    record_index: int
    display_name: str
    sequence: str
    format_name: str
    topology: str
    feature_labels: tuple[str, ...]


def normalize_sequence(sequence: str) -> str:
    normalized = "".join(sequence.split()).upper().replace("U", "T")
    invalid = set(normalized) - IUPAC_DNA
    if invalid:
        chars = "".join(sorted(invalid))
        raise ValueError(f"Sequence contains unsupported characters: {chars}")
    if not normalized:
        raise ValueError("Sequence is empty")
    return normalized


def reverse_complement(sequence: str) -> str:
    return str(Seq(sequence).reverse_complement()).upper()


def _least_rotation_index(sequence: str) -> int:
    """Return the start index of the lexicographically least rotation."""
    if not sequence:
        return 0
    doubled = sequence + sequence
    length = len(sequence)
    first, second, offset = 0, 1, 0
    while first < length and second < length and offset < length:
        left = doubled[first + offset]
        right = doubled[second + offset]
        if left == right:
            offset += 1
            continue
        if left > right:
            first = first + offset + 1
            if first <= second:
                first = second + 1
        else:
            second = second + offset + 1
            if second <= first:
                second = first + 1
        offset = 0
    return min(first, second)


def _rotate(sequence: str, offset: int) -> str:
    if not sequence:
        return sequence
    offset %= len(sequence)
    return sequence[offset:] + sequence[:offset]


def canonical_circular_sequence(sequence: str) -> str:
    """Canonicalize a circular DNA sequence across origin and strand choices."""
    forward = normalize_sequence(sequence)
    reverse = reverse_complement(forward)
    forward = _rotate(forward, _least_rotation_index(forward))
    reverse = _rotate(reverse, _least_rotation_index(reverse))
    return min(forward, reverse)


def canonical_kmer_sketch(
    sequence: str,
    *,
    k: int = 17,
    sketch_size: int = 256,
) -> tuple[int, ...]:
    """Return a strand- and origin-independent bottom-k plasmid sketch."""
    normalized = normalize_sequence(sequence)
    if k < 1 or sketch_size < 1:
        raise ValueError("k and sketch_size must be positive")
    if len(normalized) < k:
        canonical = canonical_circular_sequence(normalized).encode("ascii")
        digest = hashlib.blake2b(canonical, digest_size=8).digest()
        return (int.from_bytes(digest, "big"),)

    encoded = {"A": 0, "C": 1, "G": 2, "T": 3}
    mask = (1 << (2 * k)) - 1
    reverse_shift = 2 * (k - 1)
    forward_value = 0
    reverse_value = 0
    valid_bases = 0
    kmers: set[int] = set()
    circular = normalized + normalized[: k - 1]
    emitted = 0
    for base in circular:
        value = encoded.get(base)
        if value is None:
            forward_value = 0
            reverse_value = 0
            valid_bases = 0
            continue
        forward_value = ((forward_value << 2) | value) & mask
        reverse_value = (reverse_value >> 2) | ((3 - value) << reverse_shift)
        valid_bases += 1
        if valid_bases < k:
            continue
        kmers.add(min(forward_value, reverse_value))
        emitted += 1
        if emitted >= len(normalized):
            break
    return tuple(sorted(nsmallest(sketch_size, kmers)))


def infer_internal_sequence_role(path: Path) -> str:
    suffix = path.suffix.lower()
    if suffix in {".dna", ".gb", ".gbk"}:
        return "载体图谱"
    text = str(path).lower()
    evidence_terms = (
        "测序",
        "结果",
        "consensus",
        "成功",
        "中抽",
        "大抽",
        "ban",
        "zt0",
        "_i",
    )
    if suffix in {".fasta", ".fa", ".fna"} and any(term in text for term in evidence_terms):
        return "测序共识序列"
    return "FASTA参考序列"


def _feature_labels(record) -> tuple[str, ...]:
    labels: list[str] = []
    for feature in record.features:
        for key in ("label", "gene", "product", "name"):
            values = feature.qualifiers.get(key, [])
            if isinstance(values, str):
                values = [values]
            for value in values:
                clean = str(value).strip()
                if clean and clean not in labels:
                    labels.append(clean)
    return tuple(labels)


def _decode_text_file(path: Path) -> str:
    data = path.read_bytes()
    for encoding in ("utf-8-sig", "gb18030"):
        try:
            return data.decode(encoding)
        except UnicodeDecodeError:
            continue
    return data.decode("latin-1")


def _read_malformed_genbank(path: Path, text: str) -> list[ParsedSequence]:
    origin_match = re.search(r"(?ms)^ORIGIN\s*\r?\n(.*?)^//", text)
    if origin_match is None:
        raise ValueError(f"Malformed GenBank file has no ORIGIN sequence: {path}")
    sequence = normalize_sequence("".join(re.findall(r"[A-Za-z]", origin_match.group(1))))
    labels: list[str] = []
    for value in re.findall(r'/(?:label|gene|product|name)="([^"]+)"', text):
        clean = value.strip()
        if clean and clean not in labels:
            labels.append(clean)
    locus_line = text.splitlines()[0] if text.splitlines() else ""
    topology = "circular" if "circular" in locus_line.casefold() else "unknown"
    return [
        ParsedSequence(
            path=path,
            record_index=1,
            display_name=path.stem,
            sequence=sequence,
            format_name="genbank-fallback",
            topology=topology,
            feature_labels=tuple(labels),
        )
    ]


def read_sequence_file(path: Path) -> list[ParsedSequence]:
    path = Path(path)
    suffix = path.suffix.lower()
    format_by_suffix = {
        ".dna": "snapgene",
        ".gb": "genbank",
        ".gbk": "genbank",
        ".fasta": "fasta",
        ".fa": "fasta",
        ".fna": "fasta",
    }
    try:
        format_name = format_by_suffix[suffix]
    except KeyError as exc:
        raise ValueError(f"Unsupported sequence file type: {suffix}") from exc

    if format_name == "genbank":
        genbank_text = _decode_text_file(path)
        try:
            bio_records = list(SeqIO.parse(StringIO(genbank_text), format_name))
        except Exception:
            return _read_malformed_genbank(path, genbank_text)
    elif format_name == "snapgene":
        with path.open("rb") as handle:
            bio_records = list(SeqIO.parse(handle, format_name))
    else:
        bio_records = list(SeqIO.parse(StringIO(_decode_text_file(path)), format_name))

    parsed: list[ParsedSequence] = []
    for index, record in enumerate(bio_records, start=1):
        sequence = normalize_sequence(str(record.seq))
        name = str(record.name or record.id or "").strip()
        if not name or name.startswith("<unknown"):
            name = path.stem if index == 1 else f"{path.stem}#{index}"
        parsed.append(
            ParsedSequence(
                path=path,
                record_index=index,
                display_name=name,
                sequence=sequence,
                format_name=format_name,
                topology=str(record.annotations.get("topology", "unknown")),
                feature_labels=_feature_labels(record),
            )
        )
    if not parsed:
        raise ValueError(f"No sequence records found: {path}")
    return parsed


def _anchor_offsets(query: str, target: str, max_offsets: int = 8) -> list[int]:
    """Find likely target origins from exact anchors before edit-distance scoring."""
    if len(target) <= 256:
        return list(range(len(target)))

    votes: Counter[int] = Counter()
    doubled = target + target
    sample_step = max(1, len(query) // 96)
    for anchor_size in (31, 21, 15):
        if min(len(query), len(target)) < anchor_size:
            continue
        for query_pos in range(0, len(query) - anchor_size + 1, sample_step):
            anchor = query[query_pos:query_pos + anchor_size]
            if set(anchor) - set("ACGT"):
                continue
            search_from = 0
            matches = 0
            while matches < 4:
                target_pos = doubled.find(anchor, search_from)
                if target_pos < 0 or target_pos >= len(target):
                    break
                votes[(target_pos - query_pos) % len(target)] += anchor_size
                search_from = target_pos + 1
                matches += 1
        if votes:
            break

    offsets = [offset for offset, _ in votes.most_common(max_offsets)]
    for fallback in (0, len(target) // 4, len(target) // 2, 3 * len(target) // 4):
        if fallback not in offsets:
            offsets.append(fallback)
    return offsets[: max_offsets + 4]


def _candidate_alignments(query: str, target: str):
    for orientation, oriented in (
        ("正向", target),
        ("反向互补", reverse_complement(target)),
    ):
        for rotation in _anchor_offsets(query, oriented):
            rotated = _rotate(oriented, rotation)
            yield orientation, rotation, rotated


def compare_circular_sequences(
    internal_sequence: str,
    public_sequence: str,
    *,
    max_difference_blocks: int = 200,
) -> CircularComparison:
    """Globally compare two plasmids while allowing circular origin/strand changes."""
    internal = normalize_sequence(internal_sequence)
    public = normalize_sequence(public_sequence)

    best = None
    for orientation, rotation, rotated_public in _candidate_alignments(internal, public):
        distance = Levenshtein.distance(internal, rotated_public)
        candidate = (distance, orientation != "正向", rotation, orientation, rotated_public)
        if best is None or candidate[:3] < best[:3]:
            best = candidate
            if distance == 0:
                break

    assert best is not None
    distance, _, rotation, orientation, aligned_public = best
    edit_operations = Levenshtein.editops(internal, aligned_public)
    substitutions = sum(operation.tag == "replace" for operation in edit_operations)
    internal_extra = sum(operation.tag == "delete" for operation in edit_operations)
    internal_missing = sum(operation.tag == "insert" for operation in edit_operations)

    raw_blocks = [
        opcode
        for opcode in Levenshtein.opcodes(internal, aligned_public)
        if opcode.tag != "equal"
    ]
    blocks = tuple(
        DifferenceBlock(
            kind=opcode.tag,
            internal_start=opcode.src_start,
            internal_end=opcode.src_end,
            public_start=opcode.dest_start,
            public_end=opcode.dest_end,
            internal_sequence=internal[opcode.src_start:opcode.src_end],
            public_sequence=aligned_public[opcode.dest_start:opcode.dest_end],
        )
        for opcode in raw_blocks[:max_difference_blocks]
    )
    denominator = max(len(internal), len(public))
    identity = max(0.0, (denominator - distance) / denominator * 100)
    return CircularComparison(
        edit_distance=distance,
        identity_percent=round(identity, 4),
        substitutions=substitutions,
        internal_extra_bp=internal_extra,
        internal_missing_bp=internal_missing,
        difference_block_count=len(raw_blocks),
        difference_blocks=blocks,
        difference_blocks_truncated=len(raw_blocks) > max_difference_blocks,
        public_orientation=orientation,
        public_rotation=rotation,
    )


def classify_similarity(identity_percent: float, edit_distance: int) -> str:
    if edit_distance == 0 and identity_percent == 100.0:
        return "完全一致"
    if identity_percent >= 99.5:
        return "高度一致（小改造）"
    if identity_percent >= 95.0:
        return "同骨架高度相似"
    if identity_percent >= 80.0:
        return "部分相似"
    return "未发现可靠近似匹配"
