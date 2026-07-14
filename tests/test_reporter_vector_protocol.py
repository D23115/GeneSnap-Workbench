import unittest
from datetime import datetime, timezone

from genesnap_workbench.domain.reporter import ReporterDesignInput
from genesnap_workbench.sequence_core.dna import reverse_complement
from genesnap_workbench.sequence_core.reporter import create_reporter_design
from genesnap_workbench.vector_library.models import ReporterVectorProtocol, VectorRecord
from genesnap_workbench.vector_library.reporter import (
    apply_reporter_protocol,
    validate_reporter_protocol,
)
from tests.test_reporter_engine import promoter_sequence


NOW = datetime(2026, 7, 13, 2, 30, tzinfo=timezone.utc)
LEFT_ARM = "TCGATCGTACGTAGCTAGCTACGTA"
RIGHT_TOP = "GATCCGATCGTAGCTACGATC"


def vector_and_protocol():
    sequence = "A" * 100 + LEFT_ARM + "AGCT" + RIGHT_TOP + "C" * 100
    vector = VectorRecord.from_sequence(
        vector_record_id="artificial-gl002-like",
        structural_display_name="artificial GL002-like test vector",
        sequence=sequence,
    )
    left_boundary = 100 + len(LEFT_ARM)
    right_boundary = left_boundary + 4
    protocol = ReporterVectorProtocol(
        protocol_id="gl002-reporter-test",
        protocol_version_id="gl002-reporter-test-v1",
        display_name="GL002 reporter artificial test",
        status="enabled",
        experimental_validation_status="unverified",
        vector_record_id=vector.vector_record_id,
        vector_checksum=vector.normalized_circular_sha256,
        workflow_type="promoter_luciferase_reporter",
        insertion_mode="confirmed_interval_with_homology_prefixes",
        left_boundary=left_boundary,
        right_boundary=right_boundary,
        left_primer_homology=LEFT_ARM,
        right_primer_homology=reverse_complement(RIGHT_TOP),
        default_sequencing_method="Nanopore",
    )
    return vector, protocol


class ReporterVectorProtocolTests(unittest.TestCase):
    def test_zero_length_insertion_boundary_is_supported(self):
        sequence = "A" * 40 + LEFT_ARM + RIGHT_TOP + "C" * 40
        vector = VectorRecord.from_sequence(
            vector_record_id="zero-interval-vector",
            structural_display_name="zero interval vector",
            sequence=sequence,
        )
        boundary = 40 + len(LEFT_ARM)
        protocol = ReporterVectorProtocol(
            protocol_id="zero-interval-reporter",
            protocol_version_id="zero-interval-reporter-v1",
            display_name="zero interval reporter",
            status="enabled",
            experimental_validation_status="unverified",
            vector_record_id=vector.vector_record_id,
            vector_checksum=vector.normalized_circular_sha256,
            workflow_type="promoter_luciferase_reporter",
            insertion_mode="confirmed_interval_with_homology_prefixes",
            left_boundary=boundary,
            right_boundary=boundary,
            left_primer_homology=LEFT_ARM,
            right_primer_homology=reverse_complement(RIGHT_TOP),
        )

        self.assertTrue(validate_reporter_protocol(vector, protocol).is_valid)

    def test_protocol_validates_exact_vector_context(self):
        vector, protocol = vector_and_protocol()

        validation = validate_reporter_protocol(vector, protocol)

        self.assertTrue(validation.is_valid)
        self.assertEqual(validation.errors, ())

    def test_progressive_deletions_share_reverse_primer(self):
        vector, protocol = vector_and_protocol()
        design = create_reporter_design(
            ReporterDesignInput(
                project_id="RPT-VECTOR-001",
                gene_symbol="SGK1",
                species="human",
                promoter_sequence=promoter_sequence(2000),
                construct_lines=("WT", "P1500", "P1000", "P500"),
            ),
            protocol_version_id=protocol.protocol_version_id,
            design_version_id="RPT-VECTOR-001-v1",
            created_at=NOW,
        )

        result = apply_reporter_protocol(design, vector, protocol)

        self.assertEqual(len(result.construct_plans), 4)
        self.assertEqual(len({item.reverse_primer.sequence for item in result.construct_plans}), 1)
        self.assertEqual(
            tuple(item.forward_primer.name for item in result.construct_plans),
            ("SGK1-P2000-F", "SGK1-P1500-F", "SGK1-P1000-F", "SGK1-P500-F"),
        )
        self.assertEqual(result.construct_plans[0].reverse_primer.name, "SGK1-P-R")
        self.assertIn(
            design.constructs[0].insert_sequence,
            result.construct_plans[0].expected_plasmid_sequence,
        )

    def test_checksum_mismatch_blocks_reporter_design(self):
        vector, protocol = vector_and_protocol()
        tampered = VectorRecord.from_sequence(
            vector_record_id=vector.vector_record_id,
            structural_display_name=vector.structural_display_name,
            sequence=vector.sequence[:-1] + "G",
        )

        validation = validate_reporter_protocol(tampered, protocol)

        self.assertFalse(validation.is_valid)
        self.assertIn("vector_checksum_mismatch", validation.error_codes)


if __name__ == "__main__":
    unittest.main()
