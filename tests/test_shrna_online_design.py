import unittest
from decimal import Decimal

from genesnap_workbench.domain.shrna import BlastScreenStatus
from genesnap_workbench.integrations.broad_gpp import (
    BroadHairpinCandidate,
    BroadOligoPair,
)
from genesnap_workbench.integrations.ncbi_blast import BlastClassification
from genesnap_workbench.integrations.shrna_online import ShRNAOnlineDesigner
from genesnap_workbench.sequence_core.shrna import build_shrna_oligo_pair


TARGETS = (
    "GACTCCAGTGGTAATCTACTG",
    "GCCTCAACACGTCCCTAACTT",
    "GTGTGGTCCTGTGCACGTTTA",
    "GCTGACCTGATCGTACGATCA",
    "GATCGACCTGTTACGACTGAA",
)


class FakeBroadClient:
    def __init__(self):
        self.candidates = tuple(
            BroadHairpinCandidate(
                source_rank=index,
                start_position=position,
                intrinsic_score=Decimal(str(16 - index)),
                target_sequence=target,
                oligo_url=f"https://example.test/oligo/{index}",
            )
            for index, (position, target) in enumerate(
                zip((100, 300, 500, 520, 800), TARGETS, strict=True),
                start=1,
            )
        )

    def design_hairpins(self, sequence):
        return self.candidates

    def fetch_oligos(self, candidate):
        pair = build_shrna_oligo_pair(
            gene_symbol="TEST",
            target_no=1,
            target_id="test",
            target_sequence=candidate.target_sequence,
        )
        return BroadOligoPair(pair.forward_sequence, pair.reverse_sequence)


class FakeBlastClient:
    def screen_targets(self, candidates, *, expected_gene_symbol, species):
        results = {}
        for item in candidates:
            if item.target_sequence == TARGETS[1]:
                results[item.target_sequence] = BlastClassification(
                    status=BlastScreenStatus.FAIL,
                    note="命中其他基因 OTHER1，仅 2 个错配",
                    first_offtarget_gene="OTHER1",
                    first_offtarget_mismatches=2,
                )
            else:
                results[item.target_sequence] = BlastClassification(
                    status=BlastScreenStatus.PASS,
                    note="自动 BLAST 通过",
                )
        return results


class ShRNAOnlineDesignerTests(unittest.TestCase):
    def test_failed_initial_target_is_replaced_by_next_broad_score(self):
        designer = ShRNAOnlineDesigner(
            broad_client=FakeBroadClient(),
            blast_client=FakeBlastClient(),
            blast_batch_size=5,
        )

        result = designer.design(
            cds_sequence="ATG" * 300,
            gene_symbol="TP53",
            species="human",
            target_count=3,
        )

        self.assertEqual(
            tuple(item.target_sequence for item in result.selected_candidates),
            (TARGETS[0], TARGETS[2], TARGETS[3]),
        )
        self.assertTrue(
            all(item.blast_status is BlastScreenStatus.PASS for item in result.selected_candidates),
        )
        self.assertTrue(
            all(item.oligo_source == "broad_gpp" for item in result.selected_candidates),
        )
        self.assertFalse(result.requires_manual_confirmation)

    def test_broad_local_oligo_mismatch_requires_manual_confirmation(self):
        broad = FakeBroadClient()

        def mismatched_oligos(candidate):
            pair = build_shrna_oligo_pair(
                gene_symbol="TEST",
                target_no=1,
                target_id="test",
                target_sequence=candidate.target_sequence,
            )
            return BroadOligoPair("A" + pair.forward_sequence[1:], pair.reverse_sequence)

        broad.fetch_oligos = mismatched_oligos
        designer = ShRNAOnlineDesigner(
            broad_client=broad,
            blast_client=FakeBlastClient(),
            blast_batch_size=5,
        )

        result = designer.design(
            cds_sequence="ATG" * 300,
            gene_symbol="TP53",
            species="human",
            target_count=3,
        )

        self.assertTrue(result.requires_manual_confirmation)
        self.assertTrue(
            all("WARNING" in (item.oligo_comparison_note or "") for item in result.selected_candidates),
        )


if __name__ == "__main__":
    unittest.main()
