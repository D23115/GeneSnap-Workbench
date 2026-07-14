"""Encoding-stable GenBank output helpers."""

from __future__ import annotations

from pathlib import Path

from Bio import SeqIO
from Bio.SeqRecord import SeqRecord


def write_genbank_utf8(record: SeqRecord, path: Path) -> None:
    """Write one GenBank record without depending on the Windows code page."""

    output = Path(path)
    with output.open("w", encoding="utf-8", newline="\n") as handle:
        written = SeqIO.write(record, handle, "genbank")
    if written != 1:
        raise ValueError(f"GenBank 导出记录数异常：{written}")
