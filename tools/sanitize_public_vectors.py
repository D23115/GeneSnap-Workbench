"""将公开载体参考文件压缩为仅含序列、来源和最小 source 记录的 GenBank。"""

from __future__ import annotations

import hashlib
from dataclasses import dataclass
from pathlib import Path

from Bio import SeqIO
from Bio.SeqFeature import SeqFeature, SimpleLocation
from Bio.SeqRecord import SeqRecord


ROOT = Path(__file__).resolve().parents[1]
VECTOR_ROOT = ROOT / "src" / "genesnap_workbench" / "resources" / "vectors"


@dataclass(frozen=True)
class PublicVectorSource:
    filename: str
    record_id: str
    display_name: str
    source_url: str
    expected_sequence_sha256: str


SOURCES = (
    PublicVectorSource(
        filename="plko1_puro_snapgene_public.gb",
        record_id="pLKO1_puro_public",
        display_name="pLKO.1-puro sequence-only public reference",
        source_url=(
            "https://www.snapgene.com/plasmids/"
            "viral_expression_and_packaging_vectors/pLKO.1_puro"
        ),
        expected_sequence_sha256=(
            "13459c0789caffafc35b6f79f285da19c240a4c563c87d2da14951455d1e013e"
        ),
    ),
    PublicVectorSource(
        filename="puc57_snapgene_public.gb",
        record_id="pUC57_public",
        display_name="pUC57 sequence-only public reference",
        source_url="https://www.snapgene.com/plasmids/basic_cloning_vectors/pUC57",
        expected_sequence_sha256=(
            "3ee8c04a82fa801b4ce9b17b434208301c1b59091f7425299930f2ea96d7734f"
        ),
    ),
)


def _sequence_sha256(record: SeqRecord) -> str:
    return hashlib.sha256(str(record.seq).upper().encode("ascii")).hexdigest()


def sanitize(source: PublicVectorSource) -> None:
    path = VECTOR_ROOT / source.filename
    original = SeqIO.read(path, "genbank")
    original_hash = _sequence_sha256(original)
    if original_hash != source.expected_sequence_sha256:
        raise ValueError(f"{source.filename} 序列校验失败，拒绝重写")

    cleaned = SeqRecord(
        seq=original.seq,
        id=source.record_id,
        name=source.record_id,
        description=source.display_name,
    )
    cleaned.annotations = {
        "molecule_type": "DNA",
        "topology": "circular",
        "data_file_division": "SYN",
        "date": "01-JAN-1980",
        "comment": (
            "Sequence-only public reference. Third-party maps, notes, and "
            "annotations are not bundled. Verify the exact laboratory vector "
            f"before design. Source: {source.source_url}"
        ),
    }
    cleaned.features = [
        SeqFeature(
            SimpleLocation(0, len(cleaned.seq)),
            type="source",
            qualifiers={"label": ["sequence-only public reference"]},
        )
    ]

    temporary = path.with_name(path.name + ".tmp")
    SeqIO.write(cleaned, temporary, "genbank")
    rewritten = SeqIO.read(temporary, "genbank")
    if _sequence_sha256(rewritten) != original_hash:
        temporary.unlink(missing_ok=True)
        raise ValueError(f"{source.filename} 重写后序列发生变化")
    temporary.replace(path)
    print(f"sanitized: {source.filename} ({len(cleaned.seq)} bp)")


def main() -> None:
    for source in SOURCES:
        sanitize(source)


if __name__ == "__main__":
    main()
