import unittest
from datetime import datetime, timezone
from decimal import Decimal

from genesnap_workbench.domain.shrna import (
    BlastScreenStatus,
    ShRNACandidate,
    ShRNADesignInput,
)
from genesnap_workbench.sequence_core.shrna import (
    advance_blast_selection,
    build_shrna_oligo_pair,
    create_shrna_design,
    select_initial_candidates,
)
from genesnap_workbench.vector_library.starters import load_public_plko1_puro_starter


NOW = datetime(2026, 7, 12, 16, 30, tzinfo=timezone.utc)
VECTOR, PROTOCOL = load_public_plko1_puro_starter(user_confirmed=True)


def candidate(
    ordinal: int,
    position: int,
    score: str,
    *,
    blast_status: BlastScreenStatus = BlastScreenStatus.PENDING,
) -> ShRNACandidate:
    bases = ("ACGT" * 8)[ordinal : ordinal + 21]
    return ShRNACandidate(
        candidate_id=f"candidate-{ordinal}",
        target_sequence=bases,
        start_position=position,
        intrinsic_score=Decimal(score),
        source_rank=ordinal,
        blast_status=blast_status,
    )


class ShRNADesignTests(unittest.TestCase):
    def test_initial_selection_sorts_by_score_and_requires_strictly_over_100_bp(self):
        candidates = (
            candidate(1, 100, "8.0"),
            candidate(2, 200, "9.5"),
            candidate(3, 301, "9.0"),
            candidate(4, 450, "8.5"),
        )

        result = select_initial_candidates(candidates, target_count=3)

        self.assertEqual(
            tuple(item.candidate_id for item in result.selected),
            ("candidate-2", "candidate-3", "candidate-4"),
        )
        self.assertFalse(result.spacing_relaxed)

    def test_selection_relaxes_spacing_only_when_strict_set_is_not_enough(self):
        candidates = (
            candidate(1, 100, "9.9"),
            candidate(2, 150, "9.8"),
            candidate(3, 220, "9.7"),
            candidate(4, 280, "9.6"),
        )

        result = select_initial_candidates(candidates, target_count=3)

        self.assertEqual(len(result.selected), 3)
        self.assertEqual(result.selected[0].candidate_id, "candidate-1")
        self.assertTrue(result.spacing_relaxed)

    def test_oligo_pair_matches_confirmed_plko_hairpin_rule(self):
        pair = build_shrna_oligo_pair(
            gene_symbol="TP53",
            target_no=1,
            target_id="TP53-target-1",
            target_sequence="GACTCCAGTGGTAATCTACTG",
        )

        self.assertEqual(pair.forward_name, "TP53-1-F")
        self.assertEqual(
            pair.forward_sequence,
            "CCGGGACTCCAGTGGTAATCTACTGCTCGAGCAGTAGATTACCACTGGAGTCTTTTTG",
        )
        self.assertEqual(
            pair.reverse_sequence,
            "AATTCAAAAAGACTCCAGTGGTAATCTACTGCTCGAGCAGTAGATTACCACTGGAGTC",
        )

    def test_design_keeps_unavailable_blast_targets_but_requires_confirmation(self):
        design_input = ShRNADesignInput(
            project_id="KD-001",
            gene_symbol="TP53",
            species="human",
            cds_sequence="ATG" * 300,
            vector_protocol_version_id=PROTOCOL.protocol_version_id,
        )
        selected = tuple(
            candidate(
                ordinal,
                position,
                str(10 - ordinal),
                blast_status=BlastScreenStatus.UNAVAILABLE,
            )
            for ordinal, position in ((1, 100), (2, 250), (3, 500))
        )

        design = create_shrna_design(
            design_input,
            selected,
            VECTOR,
            PROTOCOL,
            design_version_id="KD-001-v1",
            created_at=NOW,
        )

        self.assertTrue(design.requires_confirmation)
        self.assertIn("BLAST", " ".join(design.design_warnings))
        self.assertEqual(len(design.targets), 3)
        self.assertTrue(
            all(len(target.clone_names) == 5 for target in design.targets),
        )

    def test_failed_blast_candidate_cannot_enter_formal_design(self):
        design_input = ShRNADesignInput(
            project_id="KD-001",
            gene_symbol="TP53",
            species="human",
            cds_sequence="ATG" * 300,
            vector_protocol_version_id=PROTOCOL.protocol_version_id,
            target_count=1,
        )

        with self.assertRaisesRegex(ValueError, "BLAST"):
            create_shrna_design(
                design_input,
                (candidate(1, 100, "9.0", blast_status=BlastScreenStatus.FAIL),),
                VECTOR,
                PROTOCOL,
                design_version_id="KD-001-v1",
                created_at=NOW,
            )

    def test_failed_blast_target_is_replaced_by_next_score_without_spacing_rule(self):
        current = (
            candidate(1, 100, "9.9", blast_status=BlastScreenStatus.PASS),
            candidate(2, 350, "9.8", blast_status=BlastScreenStatus.FAIL),
            candidate(3, 700, "9.7", blast_status=BlastScreenStatus.PASS),
        )
        pool = current + (
            candidate(4, 720, "9.6", blast_status=BlastScreenStatus.PENDING),
            candidate(5, 950, "9.5", blast_status=BlastScreenStatus.PENDING),
        )

        resolution = advance_blast_selection(
            current,
            pool,
            target_count=3,
        )

        self.assertEqual(
            tuple(item.candidate_id for item in resolution.selected),
            ("candidate-1", "candidate-3", "candidate-4"),
        )
        self.assertEqual(
            tuple(item.candidate_id for item in resolution.needs_screening),
            ("candidate-4",),
        )
        self.assertFalse(resolution.completed)

    def test_blast_service_unavailable_keeps_current_targets_for_manual_confirmation(self):
        current = (
            candidate(1, 100, "9.9", blast_status=BlastScreenStatus.UNAVAILABLE),
            candidate(2, 350, "9.8", blast_status=BlastScreenStatus.UNAVAILABLE),
            candidate(3, 700, "9.7", blast_status=BlastScreenStatus.UNAVAILABLE),
        )

        resolution = advance_blast_selection(
            current,
            current,
            target_count=3,
        )

        self.assertEqual(resolution.selected, current)
        self.assertTrue(resolution.requires_confirmation)
        self.assertFalse(resolution.completed)
        self.assertEqual(resolution.needs_screening, ())

    def test_all_passed_targets_complete_blast_selection(self):
        current = tuple(
            candidate(
                ordinal,
                position,
                str(10 - ordinal),
                blast_status=BlastScreenStatus.PASS,
            )
            for ordinal, position in ((1, 100), (2, 350), (3, 700))
        )

        resolution = advance_blast_selection(
            current,
            current,
            target_count=3,
        )

        self.assertTrue(resolution.completed)
        self.assertFalse(resolution.requires_confirmation)
        self.assertEqual(resolution.needs_screening, ())


if __name__ == "__main__":
    unittest.main()
