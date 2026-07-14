import unittest

from genesnap_workbench.domain.shrna import BlastScreenStatus
from genesnap_workbench.integrations.ncbi_blast import (
    BlastAlignment,
    NCBIBlastClient,
    classify_blast_alignments,
)


class NCBIBlastClassificationTests(unittest.TestCase):
    def test_same_gene_isoforms_and_other_gene_with_four_mismatches_pass(self):
        result = classify_blast_alignments(
            query_length=21,
            expected_gene_symbol="TP53",
            alignments=(
                BlastAlignment(
                    accession="NM_000546.6",
                    title="Homo sapiens tumor protein p53 (TP53), transcript variant 1, mRNA",
                    identities=21,
                    aligned_query_bases=21,
                ),
                BlastAlignment(
                    accession="NM_OTHER.1",
                    title="Homo sapiens unrelated protein (OTHER1), mRNA",
                    identities=17,
                    aligned_query_bases=21,
                ),
            ),
        )

        self.assertEqual(result.status, BlastScreenStatus.PASS)
        self.assertIn("同基因", result.note)

    def test_other_gene_with_three_mismatches_fails(self):
        result = classify_blast_alignments(
            query_length=21,
            expected_gene_symbol="TP53",
            alignments=(
                BlastAlignment(
                    accession="NM_OTHER.1",
                    title="Homo sapiens unrelated protein (OTHER1), mRNA",
                    identities=18,
                    aligned_query_bases=21,
                ),
            ),
        )

        self.assertEqual(result.status, BlastScreenStatus.FAIL)
        self.assertEqual(result.first_offtarget_gene, "OTHER1")
        self.assertEqual(result.first_offtarget_mismatches, 3)

    def test_unidentified_full_match_is_not_silently_accepted(self):
        result = classify_blast_alignments(
            query_length=21,
            expected_gene_symbol="TP53",
            alignments=(
                BlastAlignment(
                    accession="NC_000001.11",
                    title="Homo sapiens chromosome 1 genomic sequence",
                    identities=21,
                    aligned_query_bases=21,
                ),
            ),
        )

        self.assertEqual(result.status, BlastScreenStatus.FAIL)
        self.assertIn("无法识别", result.first_offtarget_gene)

    def test_related_title_does_not_hide_a_different_gene_symbol(self):
        result = classify_blast_alignments(
            query_length=21,
            expected_gene_symbol="TP53",
            alignments=(
                BlastAlignment(
                    accession="NR_000001.1",
                    title="Homo sapiens TP53 target 1 (TP53TG1), transcript",
                    identities=21,
                    aligned_query_bases=21,
                ),
            ),
        )

        self.assertEqual(result.status, BlastScreenStatus.FAIL)
        self.assertEqual(result.first_offtarget_gene, "TP53TG1")

    def test_client_batches_targets_into_one_refseq_rna_request(self):
        calls = []

        class FakeHSP:
            identities = 21
            query_start = 1
            query_end = 21

        class FakeAlignment:
            accession = "NM_000546.6"
            hit_def = "Homo sapiens tumor protein p53 (TP53), mRNA"
            hsps = (FakeHSP(),)

        class FakeRecord:
            def __init__(self, query):
                self.query = query
                self.alignments = (FakeAlignment(),)

        def fake_qblast(**kwargs):
            calls.append(kwargs)
            return object()

        client = NCBIBlastClient(
            qblast_runner=fake_qblast,
            blast_parser=lambda handle: (
                FakeRecord("target-1"),
                FakeRecord("target-2"),
            ),
        )

        results = client.screen_sequences(
            {"target-1": "GACTCCAGTGGTAATCTACTG", "target-2": "GCCTCAACACGTCCCTAACTT"},
            expected_gene_symbol="TP53",
            species="human",
        )

        self.assertEqual(len(calls), 1)
        self.assertEqual(calls[0]["database"], "refseq_rna")
        self.assertIn(">target-1", calls[0]["sequence"])
        self.assertEqual(results["target-1"].status, BlastScreenStatus.PASS)

    def test_client_uses_request_order_when_ncbi_query_labels_are_rewritten(self):
        class FakeRecord:
            def __init__(self, query):
                self.query = query
                self.alignments = ()

        client = NCBIBlastClient(
            qblast_runner=lambda **kwargs: object(),
            blast_parser=lambda handle: (FakeRecord("Query_1"), FakeRecord("Query_2")),
        )

        results = client.screen_sequences(
            {"target-1": "GACTCCAGTGGTAATCTACTG", "target-2": "GCCTCAACACGTCCCTAACTT"},
            expected_gene_symbol="TP53",
            species="human",
        )

        self.assertEqual(tuple(results), ("target-1", "target-2"))


if __name__ == "__main__":
    unittest.main()
