import tempfile
import unittest
from dataclasses import replace
from datetime import date, datetime, timezone
from pathlib import Path

from genesnap_workbench.app.application import (
    DesignConfirmationRequired,
    GeneSnapApplicationService,
    NewReporterProjectCommand,
)
from tests.test_reporter_exports import export_promoter_sequence
from tests.test_reporter_vector_protocol import vector_and_protocol


NOW = datetime(2026, 7, 13, 3, 30, tzinfo=timezone.utc)


class ReporterApplicationTests(unittest.TestCase):
    def setUp(self):
        self.temp_dir = tempfile.TemporaryDirectory()
        self.addCleanup(self.temp_dir.cleanup)
        self.service = GeneSnapApplicationService(Path(self.temp_dir.name))

    def command(self, *, mutation=False):
        vector, protocol = vector_and_protocol()
        return NewReporterProjectCommand(
            project_id="RPT-APP-001",
            gene_symbol="SGK1",
            species="human",
            promoter_sequence=export_promoter_sequence(),
            mutation_definitions=("mut1:101-104=TTTT",) if mutation else (),
            construct_lines=("mut1",) if mutation else ("WT", "P1500", "P1000", "P500"),
            received_date=date(2026, 7, 13),
            due_date=date(2026, 7, 24),
            actor="tester",
            vector=vector,
            protocol=protocol,
            clones_per_construct=5,
            gene_id="6446",
        )

    def test_create_reporter_project_saves_outputs_and_unified_summary(self):
        stored = self.service.create_reporter_project(self.command(), created_at=NOW)

        self.assertEqual(stored.snapshot.status, "design_completed")
        self.assertTrue(stored.project_folder.name.endswith("-RPT"))
        self.assertEqual(len(stored.design.constructs), 4)
        self.assertEqual(stored.design.gene_id, "6446")
        artifacts = self.service.reporter_repository.list_artifacts(stored.project_id)
        self.assertEqual(len(artifacts), 8)
        self.assertTrue(all(item.path.exists() for item in artifacts))
        summary = next(
            item for item in self.service.list_all_projects() if item.project_id == stored.project_id
        )
        self.assertEqual(summary.project_category, "报告/检测类")
        self.assertEqual(summary.workflow_type, "promoter_luciferase_reporter")
        self.assertEqual(
            summary.design_summary,
            "报告类：WT + P1500 + P1000 + P500，共 4 个构建",
        )

        reopened = GeneSnapApplicationService(Path(self.temp_dir.name))
        restored = reopened.load_reporter_project(stored.project_id)
        self.assertEqual(restored.design, stored.design)
        self.assertEqual(restored.vector_design, stored.vector_design)

    def test_reporter_mutation_requires_confirmation_reason(self):
        command = self.command(mutation=True)

        with self.assertRaises(DesignConfirmationRequired):
            self.service.create_reporter_project(command, created_at=NOW)

        stored = self.service.create_reporter_project(
            replace(command, design_confirmation_reason="已核对突变区间和最终替换序列"),
            created_at=NOW,
        )
        self.assertFalse(stored.design.requires_confirmation)
        self.assertEqual(len(stored.design.confirmation_history), 1)

    def test_reporter_protocol_profile_reopens(self):
        vector, protocol = vector_and_protocol()

        saved = self.service.save_reporter_profile(vector, protocol)
        reopened = GeneSnapApplicationService(Path(self.temp_dir.name))

        self.assertEqual(reopened.list_reporter_profiles(), (saved,))
        self.assertEqual(reopened.load_reporter_profile(saved.profile_id), (vector, protocol))

    def test_reporter_sequencing_matches_clone_and_analyzes_insert(self):
        command = replace(
            self.command(),
            construct_lines=("WT",),
        )
        stored = self.service.create_reporter_project(command, created_at=NOW)
        construct = stored.design.constructs[0]
        sequencing_dir = stored.project_folder / "03_sequencing"
        (sequencing_dir / "vendor_SGK1-promoter-WT-1_result.seq").write_text(
            "T" * 80 + construct.insert_sequence + "C" * 80,
            encoding="ascii",
        )

        outcome = self.service.analyze_reporter_sequencing(
            stored.project_id,
            actor="tester",
            analyzed_at=NOW,
        )

        latest = {item.clone_name: item for item in outcome.project.snapshot.clone_results}
        self.assertEqual(latest["SGK1-promoter-WT-1"].status, "pass")
        self.assertEqual(latest["SGK1-promoter-WT-2"].status, "warning")
        self.assertEqual(outcome.project.snapshot.status, "analysis_completed")
        self.assertTrue(outcome.analysis_report.exists())


if __name__ == "__main__":
    unittest.main()
