import io
import unittest

from Bio import SeqIO
from Bio.Seq import Seq
from Bio.SeqFeature import FeatureLocation, SeqFeature
from Bio.SeqRecord import SeqRecord

from genesnap_workbench.integrations.ncbi_transcripts import (
    NCBITranscriptClient,
    TranscriptCandidate,
)


def genbank_record(
    accession: str,
    gene: str,
    *,
    comment: str = "",
    cds: str = "ATG" + "GCT" * 30 + "TAA",
    gene_id: str | None = None,
) -> str:
    record = SeqRecord(Seq(cds), id=accession, name=accession, description=f"{gene} mRNA")
    record.annotations["molecule_type"] = "mRNA"
    record.annotations["comment"] = comment
    qualifiers = {"gene": [gene], "protein_id": ["NP_TEST.1"]}
    if gene_id is not None:
        qualifiers["db_xref"] = ["CCDS:CCDS111.1", f"GeneID:{gene_id}"]
    record.features = [
        SeqFeature(
            FeatureLocation(0, len(cds)),
            type="CDS",
            qualifiers=qualifiers,
        ),
    ]
    stream = io.StringIO()
    SeqIO.write(record, stream, "genbank")
    return stream.getvalue()


class NCBITranscriptClientTests(unittest.TestCase):
    def test_human_gene_candidates_put_mane_select_first(self):
        payload = genbank_record("NM_OTHER.1", "TP53") + genbank_record(
            "NM_MANE.2",
            "TP53",
            comment="This RefSeq is included in MANE Select.",
        )

        def fetch(endpoint, params):
            if endpoint == "esearch.fcgi":
                self.assertIn("TP53[Gene Name]", params["term"])
                self.assertIn("Homo sapiens", params["term"])
                self.assertEqual(params["retmax"], "100")
                return "<eSearchResult><IdList><Id>1</Id><Id>2</Id></IdList></eSearchResult>"
            return payload

        candidates = NCBITranscriptClient(fetch_text=fetch).find_candidates("TP53", "human")

        self.assertEqual(tuple(item.accession for item in candidates), ("NM_MANE.2", "NM_OTHER.1"))
        self.assertTrue(candidates[0].is_mane_select)
        self.assertEqual(candidates[0].cds_sequence, "ATG" + "GCT" * 30 + "TAA")

    def test_accession_lookup_returns_exact_record(self):
        payload = genbank_record(
            "NM_000546.6",
            "TP53",
            comment="MANE Select",
            gene_id="7157",
        )

        def fetch(endpoint, params):
            self.assertEqual(endpoint, "efetch.fcgi")
            self.assertEqual(params["id"], "NM_000546.6")
            return payload

        candidate = NCBITranscriptClient(fetch_text=fetch).fetch_accession("NM_000546.6")

        self.assertEqual(candidate.accession, "NM_000546.6")
        self.assertEqual(candidate.gene_symbol, "TP53")
        self.assertEqual(candidate.gene_id, "7157")
        self.assertTrue(candidate.is_mane_select)

    def test_transcript_without_geneid_xref_keeps_compatible_default(self):
        candidate = NCBITranscriptClient._parse_genbank_candidates(
            genbank_record("NM_NO_GENEID.1", "TEST"),
        )[0]

        self.assertIsNone(candidate.gene_id)

    def test_legacy_positional_candidate_construction_defaults_gene_id(self):
        candidate = TranscriptCandidate(
            "NM_LEGACY.1",
            "TEST",
            "legacy transcript",
            "ATGTAA",
            None,
            False,
            False,
        )

        self.assertIsNone(candidate.gene_id)

    def test_record_without_cds_is_ignored_and_empty_result_is_explicit(self):
        record = SeqRecord(Seq("A" * 80), id="NR_TEST.1", name="NR_TEST.1", description="RNA")
        record.annotations["molecule_type"] = "RNA"
        stream = io.StringIO()
        SeqIO.write(record, stream, "genbank")

        def fetch(endpoint, params):
            if endpoint == "esearch.fcgi":
                return "<eSearchResult><IdList><Id>1</Id></IdList></eSearchResult>"
            return stream.getvalue()

        with self.assertRaisesRegex(LookupError, "CDS"):
            NCBITranscriptClient(fetch_text=fetch).find_candidates("TEST", "mouse")


if __name__ == "__main__":
    unittest.main()
