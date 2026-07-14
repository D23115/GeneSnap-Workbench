import unittest
from datetime import datetime, timezone

from genesnap_workbench.domain.reporter import ReporterDesignInput
from genesnap_workbench.sequence_core.reporter import create_reporter_design


NOW = datetime(2026, 7, 13, 2, 0, tzinfo=timezone.utc)


def promoter_sequence(length=2000):
    return ("ACGT" * ((length + 3) // 4))[:length]


class ReporterEngineTests(unittest.TestCase):
    def test_wt_and_progressive_deletions_retain_sequence_nearest_tss(self):
        source = promoter_sequence(2000)
        design = create_reporter_design(
            ReporterDesignInput(
                project_id="RPT-001",
                gene_symbol="SGK1",
                species="human",
                promoter_sequence=source,
                construct_lines=("WT", "P1500", "P1000", "P500"),
            ),
            protocol_version_id="gl002-reporter-v1",
            design_version_id="RPT-001-v1",
            created_at=NOW,
        )

        self.assertEqual(
            tuple(item.construct_name for item in design.constructs),
            (
                "SGK1-promoter-WT",
                "SGK1-promoter-1500",
                "SGK1-promoter-1000",
                "SGK1-promoter-500",
            ),
        )
        self.assertEqual(design.constructs[1].insert_sequence, source[-1500:])
        self.assertEqual(design.constructs[-1].insert_sequence, source[-500:])

    def test_region_mutation_and_combined_mutation_use_final_replacement_sequence(self):
        source = promoter_sequence(2000)
        design = create_reporter_design(
            ReporterDesignInput(
                project_id="RPT-001",
                gene_symbol="SGK1",
                species="human",
                promoter_sequence=source,
                mutation_definitions=(
                    "mut1:101-104=TTTT",
                    "mut2:501-503=AAA",
                ),
                construct_lines=("mut1", "mut1+mut2"),
            ),
            protocol_version_id="gl002-reporter-v1",
            design_version_id="RPT-001-v1",
            created_at=NOW,
        )

        self.assertEqual(design.constructs[0].insert_sequence[100:104], "TTTT")
        self.assertEqual(design.constructs[1].insert_sequence[500:503], "AAA")
        self.assertEqual(design.constructs[1].mutation_names, ("mut1", "mut2"))
        self.assertTrue(design.requires_confirmation)

    def test_deletion_and_mutation_can_be_combined_in_one_construct(self):
        source = promoter_sequence(2000)
        design = create_reporter_design(
            ReporterDesignInput(
                project_id="RPT-001",
                gene_symbol="SGK1",
                species="human",
                promoter_sequence=source,
                mutation_definitions=("mut1:1601-1604=TTTT",),
                construct_lines=("P500+mut1",),
            ),
            protocol_version_id="gl002-reporter-v1",
            design_version_id="RPT-001-v1",
            created_at=NOW,
        )

        construct = design.constructs[0]
        self.assertEqual(len(construct.insert_sequence), 500)
        self.assertEqual(construct.insert_sequence[100:104], "TTTT")
        self.assertEqual(construct.construct_name, "SGK1-promoter-500-mut1")

    def test_mutation_outside_retained_promoter_is_rejected(self):
        with self.assertRaisesRegex(ValueError, "不在 P500"):
            create_reporter_design(
                ReporterDesignInput(
                    project_id="RPT-001",
                    gene_symbol="SGK1",
                    species="human",
                    promoter_sequence=promoter_sequence(2000),
                    mutation_definitions=("mut1:100-103=TTTT",),
                    construct_lines=("P500+mut1",),
                ),
                protocol_version_id="gl002-reporter-v1",
                design_version_id="RPT-001-v1",
                created_at=NOW,
            )


if __name__ == "__main__":
    unittest.main()
