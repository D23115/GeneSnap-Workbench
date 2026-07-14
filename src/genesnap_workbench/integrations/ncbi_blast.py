"""NCBI BLAST result classification for short shRNA targets."""

from __future__ import annotations

from dataclasses import dataclass
import re
from collections.abc import Callable, Iterable, Mapping
from io import StringIO

from Bio.Blast import NCBIWWW, NCBIXML

from genesnap_workbench.domain.shrna import BlastScreenStatus


@dataclass(frozen=True, slots=True)
class BlastAlignment:
    accession: str
    title: str
    identities: int
    aligned_query_bases: int


@dataclass(frozen=True, slots=True)
class BlastClassification:
    status: BlastScreenStatus
    note: str
    first_offtarget_gene: str | None = None
    first_offtarget_mismatches: int | None = None


class NCBIBlastError(RuntimeError):
    """Raised when the remote BLAST request or response cannot be completed."""


_GENE_PATTERNS = (
    re.compile(r"\[gene=([A-Za-z0-9_.-]+)\]", re.IGNORECASE),
    re.compile(r"\(([A-Z][A-Z0-9_.-]{1,30})\)"),
)


def _title_contains_gene(title: str, gene_symbol: str) -> bool:
    return re.search(
        rf"(?<![A-Za-z0-9]){re.escape(gene_symbol)}(?![A-Za-z0-9])",
        title,
        re.IGNORECASE,
    ) is not None


def _extract_gene(title: str) -> str | None:
    for pattern in _GENE_PATTERNS:
        matches = pattern.findall(title)
        if matches:
            return matches[-1]
    return None


def classify_blast_alignments(
    *,
    query_length: int,
    expected_gene_symbol: str,
    alignments: tuple[BlastAlignment, ...],
) -> BlastClassification:
    """Require more than three mismatches for every hit to another gene."""
    if query_length <= 0:
        raise ValueError("query_length must be positive")
    same_gene_hits = 0
    risky: list[tuple[int, str]] = []
    for alignment in alignments:
        extracted_gene = _extract_gene(alignment.title)
        if (
            extracted_gene is not None
            and extracted_gene.casefold() == expected_gene_symbol.casefold()
        ) or (
            extracted_gene is None
            and _title_contains_gene(alignment.title, expected_gene_symbol)
        ):
            same_gene_hits += 1
            continue
        mismatches = max(0, query_length - alignment.identities)
        if mismatches > 3:
            continue
        gene = extracted_gene
        if gene is None:
            gene = f"无法识别基因（{alignment.accession}）"
        risky.append((mismatches, gene))

    if risky:
        mismatches, gene = min(risky, key=lambda item: (item[0], item[1]))
        return BlastClassification(
            status=BlastScreenStatus.FAIL,
            note=f"命中其他基因 {gene}，仅 {mismatches} 个错配",
            first_offtarget_gene=gene,
            first_offtarget_mismatches=mismatches,
        )
    return BlastClassification(
        status=BlastScreenStatus.PASS,
        note=f"自动 BLAST 通过；同基因转录本命中 {same_gene_hits} 条",
    )


_SPECIES_ENTREZ = {
    "human": "Homo sapiens[Organism]",
    "mouse": "Mus musculus[Organism]",
    "rat": "Rattus norvegicus[Organism]",
}


QBlastRunner = Callable[..., object]
BlastParser = Callable[[object], Iterable[object]]


class NCBIBlastClient:
    """Batch short targets into one remote RefSeq RNA BLAST request."""

    def __init__(
        self,
        *,
        email: str | None = None,
        qblast_runner: QBlastRunner | None = None,
        blast_parser: BlastParser | None = None,
    ) -> None:
        self.email = (email or "").strip()
        self._qblast_runner = qblast_runner or NCBIWWW.qblast
        self._blast_parser = blast_parser or (lambda handle: tuple(NCBIXML.parse(handle)))

    def screen_sequences(
        self,
        sequences: Mapping[str, str],
        *,
        expected_gene_symbol: str,
        species: str,
    ) -> dict[str, BlastClassification]:
        if not sequences:
            return {}
        species_key = species.strip().lower()
        if species_key not in _SPECIES_ENTREZ:
            raise ValueError("物种必须是 human、mouse 或 rat")
        normalized: dict[str, str] = {}
        for query_id, sequence in sequences.items():
            target = "".join(sequence.split()).upper()
            if not 15 <= len(target) <= 30 or set(target) - set("ACGT"):
                raise ValueError(f"BLAST target 无效：{query_id}")
            normalized[query_id] = target
        fasta = "\n".join(f">{query_id}\n{sequence}" for query_id, sequence in normalized.items())
        if self.email:
            NCBIWWW.email = self.email
        NCBIWWW.tool = "GeneSnapWorkbench"
        try:
            handle = self._qblast_runner(
                program="blastn",
                database="refseq_rna",
                sequence=fasta,
                entrez_query=_SPECIES_ENTREZ[species_key],
                expect=1000,
                filter="F",
                hitlist_size=100,
                megablast=False,
                word_size=7,
                nucl_reward=1,
                nucl_penalty=-3,
                descriptions=100,
                alignments=100,
                format_type="XML",
            )
            records = tuple(self._blast_parser(handle))
        except Exception as error:
            raise NCBIBlastError(f"NCBI BLAST 请求或解析失败：{error}") from error
        finally:
            if "handle" in locals() and hasattr(handle, "close"):
                handle.close()

        query_ids = tuple(normalized)
        if len(records) != len(query_ids):
            raise NCBIBlastError(
                f"NCBI BLAST 返回 {len(records)} 个 query，预期 {len(query_ids)} 个",
            )
        results: dict[str, BlastClassification] = {}
        unmatched = list(query_ids)
        for record_index, record in enumerate(records):
            record_query = str(getattr(record, "query", "")).split()[0]
            query_id = (
                record_query
                if record_query in normalized and record_query in unmatched
                else unmatched[0]
            )
            if query_id in unmatched:
                unmatched.remove(query_id)
            alignments: list[BlastAlignment] = []
            for alignment in getattr(record, "alignments", ()):
                hsps = tuple(getattr(alignment, "hsps", ()))
                if not hsps:
                    continue
                hsp = max(hsps, key=lambda item: int(getattr(item, "identities", 0)))
                query_start = int(getattr(hsp, "query_start", 1))
                query_end = int(getattr(hsp, "query_end", query_start))
                alignments.append(
                    BlastAlignment(
                        accession=str(
                            getattr(alignment, "accession", None)
                            or getattr(alignment, "hit_id", "unknown")
                        ),
                        title=str(getattr(alignment, "hit_def", "")),
                        identities=int(getattr(hsp, "identities", 0)),
                        aligned_query_bases=abs(query_end - query_start) + 1,
                    ),
                )
            results[query_id] = classify_blast_alignments(
                query_length=len(normalized[query_id]),
                expected_gene_symbol=expected_gene_symbol,
                alignments=tuple(alignments),
            )
        return results

    def screen_targets(
        self,
        candidates,
        *,
        expected_gene_symbol: str,
        species: str,
    ) -> dict[str, BlastClassification]:
        by_id = {item.candidate_id: item.target_sequence for item in candidates}
        by_id_results = self.screen_sequences(
            by_id,
            expected_gene_symbol=expected_gene_symbol,
            species=species,
        )
        return {
            by_id[candidate_id]: result
            for candidate_id, result in by_id_results.items()
        }
