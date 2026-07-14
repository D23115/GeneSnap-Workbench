import unittest

from genesnap_workbench.sequence_core.dna import (
    InvalidDNASequenceError,
    normalize_dna,
    reverse_complement,
    sha256_sequence,
)
from genesnap_workbench.sequence_core.syn_qc import (
    SYNQCRules,
    classify_local_gc_percent,
    classify_repeat_length,
    evaluate_syn_sequence,
)


class DNASequenceNormalizationTests(unittest.TestCase):
    def test_normalizes_single_record_fasta_and_whitespace(self):
        raw = ">example sequence\n ac gt\nAC\tGT \n"

        self.assertEqual(normalize_dna(raw), "ACGTACGT")

    def test_rejects_empty_sequence(self):
        with self.assertRaisesRegex(ValueError, "empty"):
            normalize_dna("  \n\t")

    def test_rejects_multiple_fasta_records_in_one_design(self):
        with self.assertRaisesRegex(ValueError, "one FASTA record"):
            normalize_dna(">first\nACGT\n>second\nTGCA")

    def test_rejects_n_and_reports_invalid_positions(self):
        with self.assertRaises(InvalidDNASequenceError) as context:
            normalize_dna("ACNTX")

        self.assertEqual(context.exception.invalid_counts, {"N": 1, "X": 1})
        self.assertEqual(context.exception.invalid_positions, {"N": (2,), "X": (4,)})

    def test_checksum_is_stable_after_normalization(self):
        self.assertEqual(sha256_sequence("ac gt"), sha256_sequence("ACGT"))
        self.assertEqual(len(sha256_sequence("ACGT")), 64)

    def test_reverse_complement_uses_strict_normalization(self):
        self.assertEqual(reverse_complement("AAGCTT"), "AAGCTT")
        with self.assertRaises(InvalidDNASequenceError):
            reverse_complement("ACGN")


class SYNQCClassificationTests(unittest.TestCase):
    def test_gc_boundary_classification(self):
        expected = {
            19.9: "high_risk",
            20.0: "warning",
            24.9: "warning",
            25.0: "pass",
            75.0: "pass",
            75.1: "warning",
            80.0: "warning",
            80.1: "high_risk",
        }

        for percent, severity in expected.items():
            with self.subTest(percent=percent):
                self.assertEqual(classify_local_gc_percent(percent), severity)

    def test_adjacent_gc_windows_are_merged_into_one_risk_region(self):
        result = evaluate_syn_sequence(
            "A" * 1000,
            SYNQCRules(),
            design_version_id="design-v1",
        )

        gc_risks = [risk for risk in result.risks if risk.rule_key == "local_gc"]
        self.assertEqual(len(gc_risks), 1)
        self.assertEqual((gc_risks[0].start, gc_risks[0].end), (0, 1000))
        self.assertEqual(gc_risks[0].severity, "high_risk")

    def test_repeat_boundary_classification(self):
        expected = {
            14: "pass",
            15: "warning",
            19: "warning",
            20: "high_risk",
        }

        for length, severity in expected.items():
            with self.subTest(length=length):
                self.assertEqual(classify_repeat_length(length), severity)

    def test_homopolymer_boundaries_preserve_coordinates(self):
        sequence = "CGTC" + "A" * 8 + "CGC" + "T" * 10 + "AC" + "G" * 6

        result = evaluate_syn_sequence(
            sequence,
            SYNQCRules(),
            design_version_id="design-v1",
        )

        homopolymers = [risk for risk in result.risks if risk.rule_key == "homopolymer"]
        observed = {
            (risk.start, risk.end): (risk.severity, risk.observed_value)
            for risk in homopolymers
        }
        self.assertEqual(observed[(4, 12)], ("warning", "A x 8"))
        self.assertEqual(observed[(15, 25)], ("high_risk", "T x 10"))
        self.assertEqual(observed[(27, 33)], ("high_risk", "G x 6"))

    def test_detects_direct_and_reverse_complement_repeats(self):
        direct = "ACGTCAGTACGATCA"
        reverse = "GATCCGTAACCTGAA"
        sequence = (
            direct
            + "TTGGA"
            + direct
            + "CCATA"
            + reverse
            + "AGTCC"
            + reverse_complement(reverse)
        )

        result = evaluate_syn_sequence(
            sequence,
            SYNQCRules(),
            design_version_id="design-v1",
        )

        repeat_risks = [risk for risk in result.risks if risk.rule_key == "repeat"]
        direct_risk = next(
            risk for risk in repeat_risks if "direct" in risk.observed_value
        )
        self.assertEqual(direct_risk.severity, "warning")
        self.assertIn("count=2", direct_risk.observed_value)
        self.assertIn("positions=1,21", direct_risk.observed_value)
        self.assertTrue(
            any(
                "inverted" in risk.observed_value
                and risk.severity == "warning"
                for risk in repeat_risks
            ),
        )

    def test_internal_default_restriction_sites_are_informational(self):
        sequence = "ACGTGATATCACGTCCCGGGACGTGAATTCACGTAAGCTTACGT"

        result = evaluate_syn_sequence(
            sequence,
            SYNQCRules(),
            design_version_id="design-v1",
        )

        site_risks = [risk for risk in result.risks if risk.rule_key == "restriction_site"]
        self.assertEqual(
            {risk.observed_value for risk in site_risks},
            {"EcoRV:GATATC", "SmaI:CCCGGG", "EcoRI:GAATTC", "HindIII:AAGCTT"},
        )
        self.assertTrue(all(risk.severity == "info" for risk in site_risks))
        self.assertFalse(any(risk.requires_confirmation for risk in site_risks))

    def test_result_preserves_version_checksum_and_confirmation_reasons(self):
        sequence = "C" * 10 + "A" * 31
        rules = SYNQCRules(rules_version="syn-qc-test-v1")

        result = evaluate_syn_sequence(
            sequence,
            rules,
            design_version_id="design-v7",
        )

        self.assertEqual(result.design_version_id, "design-v7")
        self.assertEqual(result.rules_version, "syn-qc-test-v1")
        self.assertEqual(result.sequence_checksum, sha256_sequence(sequence))
        self.assertEqual(result.sequence_length, 41)
        self.assertEqual(str(result.overall_gc_percent), "24.39024390243902439024390244")
        self.assertEqual(result.blocked_reasons, ())
        self.assertTrue(result.confirmable_warnings)
        self.assertTrue(any(risk.requires_confirmation for risk in result.risks))


if __name__ == "__main__":
    unittest.main()
