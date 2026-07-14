import tempfile
import unittest
from pathlib import Path

from genesnap_workbench.sequencing.shrna import (
    CloneJudgmentStatus,
    judge_shrna_read,
    read_sequence_file,
)
from genesnap_workbench.sequence_core.shrna import build_shrna_oligo_pair


class ShRNASequencingTests(unittest.TestCase):
    def setUp(self):
        self.oligos = build_shrna_oligo_pair(
            gene_symbol="TP53",
            target_no=1,
            target_id="TP53-target-1",
            target_sequence="GACTCCAGTGGTAATCTACTG",
        )

    def test_plain_seq_file_is_read_without_changing_bases(self):
        with tempfile.TemporaryDirectory() as directory:
            path = Path(directory) / "TP53-1-1_U6_vendor.seq"
            path.write_text("AACCGG\nTTAA\n", encoding="ascii")

            read = read_sequence_file(path)

        self.assertEqual(read.sequence, "AACCGGTTAA")
        self.assertEqual(read.format_name, "seq")
        self.assertEqual(read.sample_name, "TP53-1-1_U6_vendor")

    def test_complete_forward_oligo_exact_match_is_pass(self):
        read_sequence = "A" * 80 + self.oligos.forward_sequence + "C" * 40

        judgment = judge_shrna_read(
            clone_name="TP53-1-1",
            target_id="TP53-target-1",
            read_sequence=read_sequence,
            expected_forward_oligo=self.oligos.forward_sequence,
        )

        self.assertEqual(judgment.status, CloneJudgmentStatus.PASS)
        self.assertEqual(judgment.match_start, 80)
        self.assertIn("完全匹配", judgment.reason)

    def test_single_base_mutation_is_fail_when_read_quality_is_decidable(self):
        mutated = list(self.oligos.forward_sequence)
        mutated[20] = "A" if mutated[20] != "A" else "C"
        read_sequence = "A" * 80 + "".join(mutated) + "C" * 40

        judgment = judge_shrna_read(
            clone_name="TP53-1-1",
            target_id="TP53-target-1",
            read_sequence=read_sequence,
            expected_forward_oligo=self.oligos.forward_sequence,
        )

        self.assertEqual(judgment.status, CloneJudgmentStatus.FAIL)
        self.assertIn("未找到完整", judgment.reason)

    def test_read_shorter_than_expected_oligo_is_warning(self):
        judgment = judge_shrna_read(
            clone_name="TP53-1-1",
            target_id="TP53-target-1",
            read_sequence=self.oligos.forward_sequence[:30],
            expected_forward_oligo=self.oligos.forward_sequence,
        )

        self.assertEqual(judgment.status, CloneJudgmentStatus.WARNING)
        self.assertIn("读长不足", judgment.reason)

    def test_ambiguous_bases_are_warning_when_exact_match_is_absent(self):
        read_sequence = "N" * 100 + "ACGT" * 20

        judgment = judge_shrna_read(
            clone_name="TP53-1-1",
            target_id="TP53-target-1",
            read_sequence=read_sequence,
            expected_forward_oligo=self.oligos.forward_sequence,
        )

        self.assertEqual(judgment.status, CloneJudgmentStatus.WARNING)
        self.assertIn("不确定碱基", judgment.reason)


if __name__ == "__main__":
    unittest.main()
