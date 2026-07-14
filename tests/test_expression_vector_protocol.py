import unittest
from dataclasses import replace
from datetime import datetime, timezone

from genesnap_workbench.domain.expression import ExpressionDesignInput
from genesnap_workbench.sequence_core.dna import reverse_complement
from genesnap_workbench.sequence_core.expression import create_expression_design
from genesnap_workbench.vector_library.expression import (
    apply_expression_protocol,
    expression_rules_from_protocol,
    validate_expression_protocol,
)
from genesnap_workbench.vector_library.expression_import import (
    resolve_restriction_insertion,
    scan_expression_restriction_sites,
)
from genesnap_workbench.vector_library.models import (
    ExpressionVectorProtocol,
    VectorRecord,
)


NOW = datetime(2026, 7, 12, 22, 0, tzinfo=timezone.utc)
LEFT_ARM = "ACGTCAGTACGATCGTACGATCGTAGCA"
RIGHT_TOP = "TGCATCGATGCTAGCTACGTA"


def vector_and_protocol():
    sequence = "A" * 100 + LEFT_ARM + "AGT" + RIGHT_TOP + "C" * 100
    vector = VectorRecord.from_sequence(
        vector_record_id="artificial-lv037-like",
        structural_display_name="artificial LV-037-like test vector",
        sequence=sequence,
    )
    left_boundary = 100 + len(LEFT_ARM)
    right_boundary = left_boundary + 3
    protocol = ExpressionVectorProtocol(
        protocol_id="lv037-oe-test",
        protocol_version_id="lv037-oe-test-v1",
        display_name="LV-037 OE artificial test",
        status="enabled",
        experimental_validation_status="unverified",
        vector_record_id=vector.vector_record_id,
        vector_checksum=vector.normalized_circular_sha256,
        workflow_type="expression",
        insertion_mode="confirmed_interval_with_homology_prefixes",
        left_boundary=left_boundary,
        right_boundary=right_boundary,
        left_primer_homology=LEFT_ARM,
        right_primer_homology=reverse_complement(RIGHT_TOP),
        kozak_sequence="GCCACC",
        stop_codon_rule="remove_for_c_terminal_fusion",
        c_terminal_fusion_name="3xFLAG",
    )
    return vector, protocol


class ExpressionVectorProtocolTests(unittest.TestCase):
    def test_protocol_validates_exact_vector_flanks_and_checksum(self):
        vector, protocol = vector_and_protocol()

        validation = validate_expression_protocol(vector, protocol)

        self.assertTrue(validation.is_valid)
        self.assertEqual(validation.errors, ())

    def test_protocol_generates_lv037_style_primers_and_expected_plasmid(self):
        vector, protocol = vector_and_protocol()
        design = create_expression_design(
            ExpressionDesignInput(
                project_id="OE-001",
                gene_symbol="TP53",
                species="human",
                source_cds="ATG" + "GCT" * 20 + "TAA",
                construct_lines=("FL",),
            ),
            expression_rules_from_protocol(protocol),
            design_version_id="OE-001-v1",
            created_at=NOW,
        )

        result = apply_expression_protocol(design, vector, protocol)
        plan = result.construct_plans[0]

        self.assertTrue(plan.forward_primer.startswith(LEFT_ARM + "GCCACC"))
        self.assertTrue(plan.reverse_primer.startswith(protocol.right_primer_homology))
        self.assertIn(design.constructs[0].insert_sequence, plan.expected_plasmid_sequence)
        self.assertNotIn("AGT" + RIGHT_TOP, plan.expected_plasmid_sequence)
        self.assertIn(RIGHT_TOP, plan.expected_plasmid_sequence)

    def test_single_restriction_site_resolution_validates_and_builds_expected_plasmid(self):
        vector = VectorRecord.from_sequence(
            vector_record_id="single-site-vector",
            structural_display_name="single-site expression test vector",
            sequence="A" * 30 + "ATGCATGCATGCATGCATG" + "GCTAGCGCCACC" + "A" * 50,
        )
        nhei = next(
            item
            for item in scan_expression_restriction_sites(vector.sequence)
            if item.enzyme_name == "NheI"
        )
        resolution = resolve_restriction_insertion(
            vector.sequence,
            left_occurrence=nhei,
            right_occurrence=nhei,
        )
        protocol = ExpressionVectorProtocol(
            protocol_id="single-site-oe",
            protocol_version_id="single-site-oe-v1",
            display_name="single-site OE test",
            status="enabled",
            experimental_validation_status="verified",
            vector_record_id=vector.vector_record_id,
            vector_checksum=vector.normalized_circular_sha256,
            workflow_type="expression",
            insertion_mode="single_restriction_site",
            left_boundary=resolution.left_boundary,
            right_boundary=resolution.right_boundary,
            left_primer_homology=resolution.left_primer_homology,
            right_primer_homology=resolution.right_primer_homology,
            kozak_sequence="GCCACC",
            stop_codon_rule="preserve",
            c_terminal_fusion_name=None,
        )
        design = create_expression_design(
            ExpressionDesignInput(
                project_id="OE-SINGLE-001",
                gene_symbol="TEST",
                species="human",
                source_cds="ATG" + "GCT" * 10 + "TAA",
                construct_lines=("FL",),
            ),
            expression_rules_from_protocol(protocol),
            design_version_id="OE-SINGLE-001-v1",
            created_at=NOW,
        )

        self.assertTrue(validate_expression_protocol(vector, protocol).is_valid)
        expected = apply_expression_protocol(design, vector, protocol)

        self.assertEqual(resolution.left_boundary, resolution.right_boundary)
        self.assertIn(
            design.constructs[0].insert_sequence,
            expected.construct_plans[0].expected_plasmid_sequence,
        )

    def test_protocol_checksum_mismatch_blocks_design(self):
        vector, protocol = vector_and_protocol()
        tampered = VectorRecord.from_sequence(
            vector_record_id=vector.vector_record_id,
            structural_display_name=vector.structural_display_name,
            sequence=vector.sequence[:-1] + "G",
        )

        validation = validate_expression_protocol(tampered, protocol)

        self.assertFalse(validation.is_valid)
        self.assertIn("vector_checksum_mismatch", validation.error_codes)

    def test_protocol_reports_invalid_boundaries_in_chinese(self):
        _, protocol = vector_and_protocol()

        with self.assertRaisesRegex(ValueError, "表达载体插入边界无效"):
            replace(protocol, right_boundary=protocol.left_boundary - 1)

    def test_protocol_reports_invalid_annealing_range_in_chinese(self):
        _, protocol = vector_and_protocol()

        with self.assertRaisesRegex(ValueError, "退火长度范围无效"):
            replace(protocol, anneal_min_bp=protocol.anneal_max_bp + 1)

    def test_insert_over_7000_bp_generates_four_primers_for_two_fragments(self):
        vector, protocol = vector_and_protocol()
        design = create_expression_design(
            ExpressionDesignInput(
                project_id="OE-LONG-001",
                gene_symbol="LONG1",
                species="human",
                source_cds="ATG" + "GCT" * 2400 + "TAA",
                construct_lines=("FL",),
            ),
            expression_rules_from_protocol(protocol),
            design_version_id="OE-LONG-001-v1",
            created_at=NOW,
        )

        plan = apply_expression_protocol(design, vector, protocol).construct_plans[0]

        self.assertEqual(len(plan.primers), 4)
        self.assertEqual(
            tuple(tuple(item.name.rsplit("-", 2)[-2:]) for item in plan.primers),
            (("P1", "F"), ("P1", "R"), ("P2", "F"), ("P2", "R")),
        )
        self.assertEqual(plan.primers[1].overlap_length, plan.primers[2].anneal_length)
        self.assertTrue(plan.primers[0].sequence.startswith(LEFT_ARM + "GCCACC"))
        self.assertTrue(plan.primers[-1].sequence.startswith(protocol.right_primer_homology))


if __name__ == "__main__":
    unittest.main()
