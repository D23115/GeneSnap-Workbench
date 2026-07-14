"""Broad GPP hairpin design integration using the public server-rendered forms."""

from __future__ import annotations

from dataclasses import dataclass
from decimal import Decimal, InvalidOperation
from html.parser import HTMLParser
import re
from urllib.parse import urlencode, urljoin
from urllib.request import Request, urlopen


BROAD_DESIGN_URL = "https://portals.broadinstitute.org/gpp/public/seq/search"


class BroadGPPError(RuntimeError):
    """Raised when Broad cannot be reached or its response cannot be parsed."""


@dataclass(frozen=True, slots=True)
class BroadHairpinCandidate:
    source_rank: int
    start_position: int
    intrinsic_score: Decimal
    target_sequence: str
    oligo_url: str


@dataclass(frozen=True, slots=True)
class BroadOligoPair:
    forward_sequence: str
    reverse_sequence: str


class _GridTableParser(HTMLParser):
    def __init__(self) -> None:
        super().__init__()
        self.rows: list[list[tuple[str, tuple[str, ...]]]] = []
        self._grid_depth = 0
        self._row: list[tuple[str, tuple[str, ...]]] | None = None
        self._cell_text: list[str] | None = None
        self._cell_links: list[str] = []

    def handle_starttag(self, tag: str, attrs: list[tuple[str, str | None]]) -> None:
        values = dict(attrs)
        if tag == "table":
            classes = (values.get("class") or "").split()
            if self._grid_depth or "grid" in classes:
                self._grid_depth += 1
            return
        if not self._grid_depth:
            return
        if tag == "tr":
            self._row = []
        elif tag in {"th", "td"} and self._row is not None:
            self._cell_text = []
            self._cell_links = []
        elif tag == "a" and self._cell_text is not None:
            href = values.get("href")
            if href:
                self._cell_links.append(href)

    def handle_data(self, data: str) -> None:
        if self._cell_text is not None:
            self._cell_text.append(data)

    def handle_endtag(self, tag: str) -> None:
        if not self._grid_depth:
            return
        if tag in {"th", "td"} and self._row is not None and self._cell_text is not None:
            text = " ".join("".join(self._cell_text).split())
            self._row.append((text, tuple(self._cell_links)))
            self._cell_text = None
            self._cell_links = []
        elif tag == "tr" and self._row is not None:
            if self._row:
                self.rows.append(self._row)
            self._row = None
        elif tag == "table":
            self._grid_depth -= 1


class _VisibleTextParser(HTMLParser):
    def __init__(self) -> None:
        super().__init__()
        self.parts: list[str] = []

    def handle_data(self, data: str) -> None:
        text = " ".join(data.split())
        if text:
            self.parts.append(text)


def _header_index(headers: list[str], prefix: str) -> int:
    for index, value in enumerate(headers):
        if value.casefold().startswith(prefix.casefold()):
            return index
    raise BroadGPPError(f"Broad 结果缺少字段：{prefix}")


def parse_hairpin_candidates(html: str) -> tuple[BroadHairpinCandidate, ...]:
    parser = _GridTableParser()
    parser.feed(html)
    if len(parser.rows) < 2:
        raise BroadGPPError("Broad 没有返回可解析的 hairpin 候选表")

    headers = [cell[0] for cell in parser.rows[0]]
    start_index = _header_index(headers, "Start Pos")
    score_index = _header_index(headers, "Intrinsic Score")
    target_index = _header_index(headers, "Target Sequence")
    oligo_index = _header_index(headers, "Oligo Design")
    candidates: list[BroadHairpinCandidate] = []
    for row in parser.rows[1:]:
        if max(start_index, score_index, target_index, oligo_index) >= len(row):
            continue
        target = "".join(row[target_index][0].split()).upper().replace("U", "T")
        if not 15 <= len(target) <= 30 or set(target) - set("ACGT"):
            continue
        try:
            start = int(row[start_index][0])
            score = Decimal(row[score_index][0])
        except (ValueError, InvalidOperation):
            continue
        links = row[oligo_index][1]
        oligo_url = next((urljoin(BROAD_DESIGN_URL, item) for item in links if "oligo" in item), "")
        candidates.append(
            BroadHairpinCandidate(
                source_rank=len(candidates) + 1,
                start_position=start,
                intrinsic_score=score,
                target_sequence=target,
                oligo_url=oligo_url,
            ),
        )
    if not candidates:
        raise BroadGPPError("Broad 候选表中没有有效 target")
    return tuple(candidates)


def parse_oligo_detail(html: str) -> BroadOligoPair:
    parser = _VisibleTextParser()
    parser.feed(html)
    visible_text = "\n".join(parser.parts)
    sequences = tuple(
        match.group(1).upper()
        for match in re.finditer(r"5'\s*-\s*([ACGTU\s]+?)\s*-\s*3'", visible_text, re.IGNORECASE)
    )
    normalized = tuple("".join(item.split()).replace("U", "T") for item in sequences)
    forward = next((item for item in normalized if item.startswith("CCGG") and item.endswith("TTTTTG")), "")
    reverse = next((item for item in normalized if item.startswith("AATTC")), "")
    if not forward or not reverse:
        raise BroadGPPError("Broad oligo 详情缺少正向或反向 Full sequence")
    return BroadOligoPair(forward_sequence=forward, reverse_sequence=reverse)


class BroadGPPClient:
    def __init__(self, *, timeout_seconds: int = 60) -> None:
        self.timeout_seconds = timeout_seconds

    @staticmethod
    def _request_text(request: Request, timeout_seconds: int) -> str:
        try:
            with urlopen(request, timeout=timeout_seconds) as response:
                return response.read().decode("utf-8", "replace")
        except OSError as error:
            raise BroadGPPError(f"Broad GPP 请求失败：{error}") from error

    def design_hairpins(self, transcript_sequence: str) -> tuple[BroadHairpinCandidate, ...]:
        sequence = "".join(transcript_sequence.split()).upper().replace("U", "T")
        if len(sequence) < 50 or set(sequence) - set("ACGT"):
            raise ValueError("提交 Broad 的转录本序列必须是至少 50 nt 的 A/C/G/T 序列")
        payload = urlencode({"data_entered": sequence, "Design": "Design »"}).encode()
        request = Request(
            BROAD_DESIGN_URL,
            data=payload,
            headers={"User-Agent": "GeneSnapWorkbench/0.3"},
        )
        return parse_hairpin_candidates(self._request_text(request, self.timeout_seconds))

    def fetch_oligos(self, candidate: BroadHairpinCandidate) -> BroadOligoPair:
        if not candidate.oligo_url:
            raise BroadGPPError("Broad 候选没有 oligo 详情链接")
        request = Request(candidate.oligo_url, headers={"User-Agent": "GeneSnapWorkbench/0.3"})
        return parse_oligo_detail(self._request_text(request, self.timeout_seconds))
