import unittest
from pathlib import Path
import struct
from tempfile import TemporaryDirectory

from genesnap_workbench.vector_library.comparison import (
    canonical_circular_sequence,
    canonical_kmer_sketch,
    classify_similarity,
    compare_circular_sequences,
    infer_internal_sequence_role,
    read_sequence_file,
)


class CircularSequenceComparisonTests(unittest.TestCase):
    def test_canonical_sequence_ignores_origin_and_strand(self):
        sequence = "AACCGTTAAGTC"
        rotated = sequence[5:] + sequence[:5]
        reverse_complement = "GACTTAACGGTT"

        expected = canonical_circular_sequence(sequence)

        self.assertEqual(canonical_circular_sequence(rotated), expected)
        self.assertEqual(canonical_circular_sequence(reverse_complement), expected)

    def test_exact_rotated_sequence_has_zero_differences(self):
        internal = "AACCGTTAAGTC"
        public = internal[4:] + internal[:4]

        result = compare_circular_sequences(internal, public)

        self.assertEqual(result.edit_distance, 0)
        self.assertEqual(result.identity_percent, 100.0)
        self.assertEqual(result.substitutions, 0)
        self.assertEqual(result.internal_extra_bp, 0)
        self.assertEqual(result.internal_missing_bp, 0)

    def test_reports_substitution_and_internal_insertion(self):
        public = "AACCGGTTTACG"
        internal = "AATCGGATTTACG"

        result = compare_circular_sequences(internal, public)

        self.assertEqual(result.edit_distance, 2)
        self.assertEqual(result.substitutions, 1)
        self.assertEqual(result.internal_extra_bp, 1)
        self.assertEqual(result.internal_missing_bp, 0)

    def test_similarity_classification_uses_stable_thresholds(self):
        self.assertEqual(classify_similarity(100.0, 0), "完全一致")
        self.assertEqual(classify_similarity(99.8, 2), "高度一致（小改造）")
        self.assertEqual(classify_similarity(97.0, 100), "同骨架高度相似")
        self.assertEqual(classify_similarity(85.0, 500), "部分相似")
        self.assertEqual(classify_similarity(70.0, 1000), "未发现可靠近似匹配")

    def test_kmer_sketch_ignores_circular_origin_and_strand(self):
        sequence = "AACCGTTAAGTCGGATCC"
        rotated = sequence[7:] + sequence[:7]
        reverse_complement = "GGATCCGACTTAACGGTT"

        expected = canonical_kmer_sketch(sequence, k=5, sketch_size=100)

        self.assertEqual(canonical_kmer_sketch(rotated, k=5, sketch_size=100), expected)
        self.assertEqual(
            canonical_kmer_sketch(reverse_complement, k=5, sketch_size=100),
            expected,
        )

    def test_short_sequence_sketch_is_deterministic(self):
        self.assertEqual(
            canonical_kmer_sketch("AAC", k=5, sketch_size=100),
            (8856309073986437891,),
        )

    def test_internal_fasta_under_result_folder_is_sequencing_consensus(self):
        path = Path("慢病毒") / "LV-037" / "二代测序结果" / "result.fasta"
        vendor_path = Path("GL014") / "GL014--ZT007788_FX001_i28.1.fasta"

        self.assertEqual(infer_internal_sequence_role(path), "测序共识序列")
        self.assertEqual(infer_internal_sequence_role(vendor_path), "测序共识序列")
        self.assertEqual(infer_internal_sequence_role(Path("LV-037.dna")), "载体图谱")

    def test_reads_fasta_sequence_with_metadata(self):
        with TemporaryDirectory() as tmpdir:
            path = Path(tmpdir) / "example.fasta"
            path.write_text(">verified_vector\naaccttgg\n", encoding="ascii")

            records = read_sequence_file(path)

        self.assertEqual(len(records), 1)
        self.assertEqual(records[0].display_name, "verified_vector")
        self.assertEqual(records[0].sequence, "AACCTTGG")
        self.assertEqual(records[0].format_name, "fasta")

    def test_reads_utf8_fasta_with_chinese_display_name(self):
        with TemporaryDirectory() as tmpdir:
            path = Path(tmpdir) / "中文载体.fasta"
            path.write_text(">已确认载体\naaccttgg\n", encoding="utf-8")

            records = read_sequence_file(path)

        self.assertEqual(records[0].display_name, "已确认载体")
        self.assertEqual(records[0].sequence, "AACCTTGG")

    def test_reads_snapgene_dna_in_binary_mode(self):
        cookie = struct.pack(">8sHHH", b"SnapGene", 1, 1, 1)
        sequence = b"AACCGGTTAACC"
        dna_packet = b"\x01" + sequence
        content = (
            struct.pack(">BI", 0x09, len(cookie))
            + cookie
            + struct.pack(">BI", 0x00, len(dna_packet))
            + dna_packet
        )
        with TemporaryDirectory() as tmpdir:
            path = Path(tmpdir) / "LV-037.dna"
            path.write_bytes(content)

            records = read_sequence_file(path)

        self.assertEqual(len(records), 1)
        self.assertEqual(records[0].sequence, sequence.decode("ascii"))
        self.assertEqual(records[0].topology, "circular")
        self.assertEqual(records[0].format_name, "snapgene")

    def test_reads_sequence_from_vendor_genbank_with_malformed_locus(self):
        genbank = """LOCUS       Exported File           12 bp ds-DNA    circular SYN 10-03-2021
FEATURES             Location/Qualifiers
     misc_feature    1..4
                     /label="MCS"
ORIGIN
        1 aacc ggtt aacc
//
"""
        with TemporaryDirectory() as tmpdir:
            path = Path(tmpdir) / "vendor_export.gb"
            path.write_text(genbank, encoding="ascii")

            records = read_sequence_file(path)

        self.assertEqual(len(records), 1)
        self.assertEqual(records[0].display_name, "vendor_export")
        self.assertEqual(records[0].sequence, "AACCGGTTAACC")
        self.assertEqual(records[0].topology, "circular")
        self.assertEqual(records[0].feature_labels, ("MCS",))

    def test_preserves_chinese_labels_in_gb18030_genbank(self):
        genbank = """LOCUS       Exported File           12 bp ds-DNA    circular SYN 10-03-2021
FEATURES             Location/Qualifiers
     promoter        1..4
                     /label="启动子"
ORIGIN
        1 aacc ggtt aacc
//
"""
        with TemporaryDirectory() as tmpdir:
            path = Path(tmpdir) / "vendor_chinese.gb"
            path.write_bytes(genbank.encode("gb18030"))

            records = read_sequence_file(path)

        self.assertEqual(records[0].feature_labels, ("启动子",))

    def test_missing_genbank_preserves_file_not_found_error(self):
        with TemporaryDirectory() as tmpdir:
            path = Path(tmpdir) / "missing.gb"

            with self.assertRaises(FileNotFoundError):
                read_sequence_file(path)


if __name__ == "__main__":
    unittest.main()
