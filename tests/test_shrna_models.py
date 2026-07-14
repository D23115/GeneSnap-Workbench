import unittest
from dataclasses import FrozenInstanceError
from datetime import datetime, timezone
from decimal import Decimal

from genesnap_workbench.domain.shrna import (
    BlastScreenStatus,
    ShRNACandidate,
    ShRNADesignInput,
    ShRNAOligoPair,
    ShRNATargetDesign,
)


NOW = datetime(2026, 7, 12, 16, 0, tzinfo=timezone.utc)


class ShRNAModelTests(unittest.TestCase):
    def test_design_input_defaults_to_three_targets_and_five_clones(self):
        design_input = ShRNADesignInput(
            project_id="KD-001",
            gene_symbol="TP53",
            species="human",
            cds_sequence="ATG" * 100,
            vector_protocol_version_id="plko-public-v1",
        )

        self.assertEqual(design_input.target_count, 3)
        self.assertEqual(design_input.clones_per_target, 5)

    def test_target_count_is_limited_to_current_business_range(self):
        with self.assertRaisesRegex(ValueError, "target_count"):
            ShRNADesignInput(
                project_id="KD-001",
                gene_symbol="TP53",
                species="human",
                cds_sequence="ATG" * 100,
                vector_protocol_version_id="plko-public-v1",
                target_count=4,
            )

    def test_candidate_normalizes_sequence_and_is_immutable(self):
        candidate = ShRNACandidate(
            candidate_id="cand-1",
            target_sequence="acgtacgtacgtacgtacgta",
            start_position=120,
            intrinsic_score=Decimal("8.7"),
            source_rank=1,
            blast_status=BlastScreenStatus.PENDING,
        )

        self.assertEqual(candidate.target_sequence, "ACGTACGTACGTACGTACGTA")
        with self.assertRaises(FrozenInstanceError):
            candidate.start_position = 300

    def test_target_requires_matching_oligo_and_clone_ownership(self):
        candidate = ShRNACandidate(
            candidate_id="cand-1",
            target_sequence="ACGTACGTACGTACGTACGTA",
            start_position=120,
            intrinsic_score=Decimal("8.7"),
            source_rank=1,
        )
        pair = ShRNAOligoPair(
            target_id="target-1",
            forward_name="TP53-1-F",
            forward_sequence="CCGGACGTACGTACGTACGTACGTACTCGAGTACGTACGTACGTACGTACGTTTTTTG",
            reverse_name="TP53-1-R",
            reverse_sequence="AATTCAAAAATACGTACGTACGTACGTACGTCTCGAGACGTACGTACGTACGTACGTA",
        )

        target = ShRNATargetDesign(
            target_id="target-1",
            target_no=1,
            candidate=candidate,
            oligos=pair,
            clone_names=("TP53-1-1", "TP53-1-2"),
        )

        self.assertEqual(target.clone_names[-1], "TP53-1-2")


if __name__ == "__main__":
    unittest.main()
