import unittest

from genesnap_workbench.domain.shrna import BlastScreenStatus
from genesnap_workbench.sequence_core.shrna_candidates import generate_shrna_candidates


class ShRNACandidateGeneratorTests(unittest.TestCase):
    def test_generates_ranked_candidates_with_positions_and_qc(self):
        cds = ("ATGCGTACGTTAGCTACGATCGTACCTGAGTCAAGT" * 20)[:700]

        candidates = generate_shrna_candidates(cds, max_candidates=30)

        self.assertGreaterEqual(len(candidates), 3)
        self.assertEqual(
            tuple(item.intrinsic_score for item in candidates),
            tuple(sorted((item.intrinsic_score for item in candidates), reverse=True)),
        )
        self.assertTrue(all(len(item.target_sequence) == 21 for item in candidates))
        self.assertTrue(all(item.start_position and item.start_position >= 26 for item in candidates))
        self.assertTrue(all(item.blast_status is BlastScreenStatus.UNAVAILABLE for item in candidates))
        self.assertTrue(all(30 <= item.gc_percent <= 60 for item in candidates))

    def test_extreme_gc_and_homopolymers_are_filtered(self):
        with self.assertRaisesRegex(ValueError, "候选"):
            generate_shrna_candidates("G" * 400)

    def test_duplicate_target_sequences_are_not_returned_twice(self):
        cds = "ATG" + ("ACGT" * 150)

        candidates = generate_shrna_candidates(cds)

        sequences = tuple(item.target_sequence for item in candidates)
        self.assertEqual(len(sequences), len(set(sequences)))

    def test_accepts_single_record_fasta(self):
        candidates = generate_shrna_candidates(">NM_TEST example\n" + "ACGT" * 150)

        self.assertGreaterEqual(len(candidates), 3)


if __name__ == "__main__":
    unittest.main()
