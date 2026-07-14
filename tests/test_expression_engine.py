import unittest
from datetime import datetime, timezone

from genesnap_workbench.domain.expression import (
    ExpressionDesignInput,
    ExpressionVectorRules,
)
from genesnap_workbench.sequence_core.expression import create_expression_design


NOW = datetime(2026, 7, 12, 21, 30, tzinfo=timezone.utc)


def coding_sequence(amino_acids: int, *, terminal_stop: bool = True) -> str:
    sequence = "ATG" + "GCT" * (amino_acids - 1)
    return sequence + ("TAA" if terminal_stop else "")


class ExpressionEngineTests(unittest.TestCase):
    def rules(self, *, remove_stop: bool = True):
        return ExpressionVectorRules(
            protocol_version_id="lv037-test-v1",
            kozak_sequence="GCCACC",
            stop_codon_rule=(
                "remove_for_c_terminal_fusion" if remove_stop else "preserve"
            ),
            c_terminal_fusion_name="3xFLAG" if remove_stop else None,
            single_fragment_max_bp=7000,
        )

    def test_one_project_parses_full_length_truncation_and_deletion_lines(self):
        design = create_expression_design(
            ExpressionDesignInput(
                project_id="OE-001",
                gene_symbol="TP53",
                species="human",
                source_cds=coding_sequence(700),
                construct_lines=("FL", "1-300aa", "Δ301-600"),
            ),
            self.rules(),
            design_version_id="OE-001-v1",
            created_at=NOW,
        )

        self.assertEqual(
            tuple(item.construct_name for item in design.constructs),
            ("TP53-FL", "TP53-1-300aa", "TP53-Δ301-600"),
        )
        self.assertEqual(len(design.constructs[1].coding_sequence), 900)
        self.assertEqual(len(design.constructs[2].coding_sequence), 1200)

    def test_c_terminal_fusion_removes_stop_and_adds_editable_kozak(self):
        design = create_expression_design(
            ExpressionDesignInput(
                project_id="OE-001",
                gene_symbol="TP53",
                species="human",
                source_cds=coding_sequence(20),
                construct_lines=("FL",),
            ),
            self.rules(),
            design_version_id="OE-001-v1",
            created_at=NOW,
        )

        construct = design.constructs[0]
        self.assertFalse(construct.coding_sequence.endswith("TAA"))
        self.assertTrue(construct.insert_sequence.startswith("GCCACCATG"))
        self.assertEqual(construct.c_terminal_fusion_name, "3xFLAG")

    def test_deleting_amino_acid_one_reintroduces_start_codon(self):
        design = create_expression_design(
            ExpressionDesignInput(
                project_id="OE-001",
                gene_symbol="TP53",
                species="human",
                source_cds=coding_sequence(400),
                construct_lines=("Δ1-200",),
            ),
            self.rules(),
            design_version_id="OE-001-v1",
            created_at=NOW,
        )

        construct = design.constructs[0]
        self.assertTrue(construct.coding_sequence.startswith("ATG"))
        self.assertTrue(construct.start_codon_reintroduced)
        self.assertEqual(len(construct.coding_sequence), 3 + 200 * 3)

    def test_insert_over_7000_bp_uses_two_contiguous_pcr_fragments(self):
        design = create_expression_design(
            ExpressionDesignInput(
                project_id="OE-001",
                gene_symbol="LONG1",
                species="human",
                source_cds=coding_sequence(2500),
                construct_lines=("FL",),
            ),
            self.rules(),
            design_version_id="OE-001-v1",
            created_at=NOW,
        )

        fragments = design.constructs[0].fragments
        self.assertEqual(len(fragments), 2)
        self.assertEqual(fragments[0].start, 0)
        self.assertEqual(fragments[0].end, fragments[1].start)
        self.assertEqual(fragments[-1].end, len(design.constructs[0].insert_sequence))

    def test_unparseable_construct_line_is_preserved_for_manual_review(self):
        design = create_expression_design(
            ExpressionDesignInput(
                project_id="OE-001",
                gene_symbol="TP53",
                species="human",
                source_cds=coding_sequence(400),
                construct_lines=("特殊标签版本A",),
            ),
            self.rules(),
            design_version_id="OE-001-v1",
            created_at=NOW,
        )

        self.assertTrue(design.requires_confirmation)
        self.assertEqual(design.unparsed_lines, ("特殊标签版本A",))

    def test_point_mutation_checks_reference_amino_acid_and_generates_sequence(self):
        source = "ATG" + "AAG" + "GCT" * 18 + "TAA"
        design = create_expression_design(
            ExpressionDesignInput(
                project_id="OE-001",
                gene_symbol="TP53",
                species="human",
                source_cds=source,
                construct_lines=("K2R",),
            ),
            self.rules(),
            design_version_id="OE-001-v1",
            created_at=NOW,
        )

        construct = design.constructs[0]
        self.assertEqual(construct.construct_name, "TP53-K2R")
        self.assertEqual(construct.mutations[0].original_codon, "AAG")
        self.assertEqual(construct.mutations[0].new_amino_acid, "R")
        self.assertTrue(design.requires_confirmation)

    def test_point_mutation_rejects_wrong_reference_amino_acid(self):
        with self.assertRaisesRegex(ValueError, "observed"):
            create_expression_design(
                ExpressionDesignInput(
                    project_id="OE-001",
                    gene_symbol="TP53",
                    species="human",
                    source_cds=coding_sequence(20),
                    construct_lines=("K2R",),
                ),
                self.rules(),
                design_version_id="OE-001-v1",
                created_at=NOW,
            )


if __name__ == "__main__":
    unittest.main()
