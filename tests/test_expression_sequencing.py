import unittest

from genesnap_workbench.sequencing.expression import (
    ExpressionCloneJudgmentStatus,
    judge_expression_read,
)


class ExpressionSequencingTests(unittest.TestCase):
    def test_complete_insert_inside_whole_plasmid_read_is_pass(self):
        expected = "GCCACC" + "ATG" + "GCT" * 30
        result = judge_expression_read(
            clone_name="TP53-FL-1",
            construct_id="construct-1",
            read_sequence="T" * 80 + expected + "C" * 80,
            expected_insert_sequence=expected,
            expected_coding_sequence=expected[6:],
        )

        self.assertEqual(result.status, ExpressionCloneJudgmentStatus.PASS)
        self.assertEqual(result.substitution_count, 0)
        self.assertEqual(result.deletion_count, 0)
        self.assertFalse(result.frameshift)

    def test_reverse_complement_whole_plasmid_read_is_pass(self):
        expected = "GCCACC" + "ATG" + "GCT" * 30
        from genesnap_workbench.sequence_core.dna import reverse_complement

        result = judge_expression_read(
            clone_name="TP53-FL-1",
            construct_id="construct-1",
            read_sequence=reverse_complement("T" * 80 + expected + "C" * 80),
            expected_insert_sequence=expected,
            expected_coding_sequence=expected[6:],
        )

        self.assertEqual(result.status, ExpressionCloneJudgmentStatus.PASS)
        self.assertEqual(result.orientation, "reverse")

    def test_single_substitution_requires_manual_review_and_reports_count(self):
        expected = "GCCACC" + "ATG" + "GCT" * 30
        observed = list(expected)
        observed[35] = "A" if observed[35] != "A" else "C"

        result = judge_expression_read(
            clone_name="TP53-FL-1",
            construct_id="construct-1",
            read_sequence="T" * 80 + "".join(observed) + "C" * 80,
            expected_insert_sequence=expected,
            expected_coding_sequence=expected[6:],
        )

        self.assertEqual(result.status, ExpressionCloneJudgmentStatus.WARNING)
        self.assertEqual(result.substitution_count, 1)
        self.assertEqual(result.deletion_count, 0)
        self.assertIn("1 个碱基替换", result.reason)

    def test_one_base_deletion_reports_frameshift(self):
        expected = "GCCACC" + "ATG" + "GCT" * 30
        observed = expected[:45] + expected[46:]

        result = judge_expression_read(
            clone_name="TP53-FL-1",
            construct_id="construct-1",
            read_sequence="T" * 80 + observed + "C" * 80,
            expected_insert_sequence=expected,
            expected_coding_sequence=expected[6:],
        )

        self.assertEqual(result.status, ExpressionCloneJudgmentStatus.WARNING)
        self.assertEqual(result.deletion_count, 1)
        self.assertTrue(result.frameshift)
        self.assertIn("移码", result.reason)


if __name__ == "__main__":
    unittest.main()
