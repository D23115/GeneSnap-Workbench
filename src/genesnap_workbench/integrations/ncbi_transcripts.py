"""NCBI E-utilities transcript lookup with Biopython GenBank parsing."""

from __future__ import annotations

from collections.abc import Callable
from dataclasses import dataclass
from io import StringIO
from urllib.parse import urlencode
from urllib.request import Request, urlopen
from xml.etree import ElementTree

from Bio import SeqIO


_EUTILS_BASE = "https://eutils.ncbi.nlm.nih.gov/entrez/eutils"
_SPECIES_NAMES = {
    "human": "Homo sapiens",
    "mouse": "Mus musculus",
    "rat": "Rattus norvegicus",
}


@dataclass(frozen=True, slots=True)
class TranscriptCandidate:
    accession: str
    gene_symbol: str
    description: str
    cds_sequence: str
    protein_id: str | None
    is_mane_select: bool
    is_refseq_select: bool

    @property
    def display_label(self) -> str:
        flags: list[str] = []
        if self.is_mane_select:
            flags.append("MANE Select")
        elif self.is_refseq_select:
            flags.append("RefSeq Select")
        suffix = f" [{', '.join(flags)}]" if flags else ""
        return f"{self.accession} | CDS {len(self.cds_sequence)} bp{suffix}"


FetchText = Callable[[str, dict[str, str]], str]


def _default_fetch_text(endpoint: str, params: dict[str, str]) -> str:
    url = f"{_EUTILS_BASE}/{endpoint}?{urlencode(params)}"
    request = Request(url, headers={"User-Agent": "GeneSnapWorkbench/0.1"})
    try:
        with urlopen(request, timeout=30) as response:
            return response.read().decode("utf-8")
    except OSError as error:
        raise ConnectionError(f"NCBI 请求失败：{error}") from error


class NCBITranscriptClient:
    def __init__(
        self,
        *,
        fetch_text: FetchText | None = None,
        email: str | None = None,
        api_key: str | None = None,
    ) -> None:
        self._fetch_text = fetch_text or _default_fetch_text
        self._email = email
        self._api_key = api_key

    def _params(self, **values: str) -> dict[str, str]:
        params = {"tool": "GeneSnapWorkbench", **values}
        if self._email:
            params["email"] = self._email
        if self._api_key:
            params["api_key"] = self._api_key
        return params

    def find_candidates(
        self,
        gene_symbol: str,
        species: str,
        *,
        maximum: int = 100,
    ) -> tuple[TranscriptCandidate, ...]:
        symbol = gene_symbol.strip()
        species_key = species.strip().lower()
        if not symbol:
            raise ValueError("基因名不能为空")
        if species_key not in _SPECIES_NAMES:
            raise ValueError("物种必须是 human、mouse 或 rat")
        if maximum <= 0:
            raise ValueError("maximum 必须大于 0")

        organism = _SPECIES_NAMES[species_key]
        search_text = self._fetch_text(
            "esearch.fcgi",
            self._params(
                db="nuccore",
                term=(
                    f"{symbol}[Gene Name] AND \"{organism}\"[Organism] "
                    "AND srcdb_refseq[PROP] AND biomol_mrna[PROP]"
                ),
                retmode="xml",
                retmax=str(maximum),
            ),
        )
        try:
            root = ElementTree.fromstring(search_text)
        except ElementTree.ParseError as error:
            raise ConnectionError("NCBI 返回的转录本检索结果无法解析") from error
        ids = tuple(
            node.text.strip()
            for node in root.findall(".//IdList/Id")
            if node.text and node.text.strip()
        )
        if not ids:
            raise LookupError(f"NCBI 未找到 {symbol} 的 RefSeq mRNA 转录本")

        payload = self._fetch_text(
            "efetch.fcgi",
            self._params(
                db="nuccore",
                id=",".join(ids),
                rettype="gb",
                retmode="text",
            ),
        )
        candidates = self._parse_genbank_candidates(payload)
        if not candidates:
            raise LookupError(f"NCBI 返回的 {symbol} 转录本中没有可用 CDS")
        return tuple(
            sorted(
                candidates,
                key=lambda item: (
                    not item.is_mane_select,
                    not item.is_refseq_select,
                    item.accession,
                ),
            ),
        )

    def fetch_accession(self, accession: str) -> TranscriptCandidate:
        normalized = accession.strip()
        if not normalized:
            raise ValueError("转录本号不能为空")
        payload = self._fetch_text(
            "efetch.fcgi",
            self._params(
                db="nuccore",
                id=normalized,
                rettype="gb",
                retmode="text",
            ),
        )
        candidates = self._parse_genbank_candidates(payload)
        if not candidates:
            raise LookupError(f"转录本 {normalized} 没有可用 CDS")
        exact = next(
            (
                item
                for item in candidates
                if item.accession.upper() == normalized.upper()
                or item.accession.split(".", 1)[0].upper()
                == normalized.split(".", 1)[0].upper()
            ),
            candidates[0],
        )
        return exact

    @staticmethod
    def _parse_genbank_candidates(payload: str) -> tuple[TranscriptCandidate, ...]:
        candidates: list[TranscriptCandidate] = []
        for record in SeqIO.parse(StringIO(payload), "genbank"):
            searchable_parts = [record.description]
            searchable_parts.extend(str(value) for value in record.annotations.values())
            searchable = " ".join(searchable_parts).lower()
            for feature in record.features:
                if feature.type != "CDS":
                    continue
                qualifiers = feature.qualifiers
                qualifier_text = " ".join(
                    value
                    for values in qualifiers.values()
                    for value in values
                ).lower()
                combined = f"{searchable} {qualifier_text}"
                gene = (qualifiers.get("gene") or [""])[0]
                protein_id = (qualifiers.get("protein_id") or [None])[0]
                candidates.append(
                    TranscriptCandidate(
                        accession=record.id,
                        gene_symbol=gene,
                        description=record.description,
                        cds_sequence=str(feature.extract(record.seq)).upper(),
                        protein_id=protein_id,
                        is_mane_select="mane select" in combined,
                        is_refseq_select="refseq select" in combined,
                    ),
                )
                break
        return tuple(candidates)
