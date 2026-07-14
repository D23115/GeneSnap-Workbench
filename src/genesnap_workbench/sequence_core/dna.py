"""Strict DNA normalization and deterministic sequence primitives."""

from __future__ import annotations

from collections import defaultdict
import hashlib


class InvalidDNASequenceError(ValueError):
    """Report unsupported DNA characters and their normalized coordinates."""

    def __init__(
        self,
        invalid_counts: dict[str, int],
        invalid_positions: dict[str, tuple[int, ...]],
    ) -> None:
        self.invalid_counts = invalid_counts
        self.invalid_positions = invalid_positions
        details = "; ".join(
            f"{base} x{invalid_counts[base]} at {invalid_positions[base]}"
            for base in sorted(invalid_counts)
        )
        super().__init__(f"Sequence contains unsupported DNA characters: {details}")


def normalize_dna(raw_sequence: str) -> str:
    """Return uppercase A/C/G/T DNA after removing one FASTA header and whitespace."""
    lines = raw_sequence.lstrip("\ufeff").splitlines()
    header_count = sum(line.lstrip().startswith(">") for line in lines)
    if header_count > 1:
        raise ValueError("A design must contain exactly one FASTA record")
    sequence_lines = [
        line
        for line in lines
        if not line.lstrip().startswith(">")
    ]
    normalized = "".join("".join(sequence_lines).split()).upper()
    if not normalized:
        raise ValueError("DNA sequence is empty")

    invalid_positions: dict[str, list[int]] = defaultdict(list)
    for position, base in enumerate(normalized):
        if base not in "ACGT":
            invalid_positions[base].append(position)
    if invalid_positions:
        raise InvalidDNASequenceError(
            {
                base: len(positions)
                for base, positions in invalid_positions.items()
            },
            {
                base: tuple(positions)
                for base, positions in invalid_positions.items()
            },
        )
    return normalized


def reverse_complement(sequence: str) -> str:
    """Return the reverse complement of strict DNA."""
    normalized = normalize_dna(sequence)
    return normalized.translate(str.maketrans("ACGT", "TGCA"))[::-1]


def sha256_sequence(sequence: str) -> str:
    """Return a stable SHA-256 checksum for normalized DNA."""
    normalized = normalize_dna(sequence)
    return hashlib.sha256(normalized.encode("ascii")).hexdigest()
