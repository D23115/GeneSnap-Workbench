import tempfile
import unittest
from dataclasses import replace
from datetime import date, datetime, timezone
from pathlib import Path

from genesnap_workbench.app.application import (
    DesignConfirmationRequired,
    GeneSnapApplicationService,
    NewExpressionProjectCommand,
)
from tests.test_expression_vector_protocol import vector_and_protocol


NOW = datetime(2026, 7, 13, 0, 0, tzinfo=timezone.utc)


class ExpressionApplicationTests(unittest.TestCase):
    def setUp(self):
        self.temp_dir = tempfile.TemporaryDirectory()
        self.addCleanup(self.temp_dir.cleanup)
        self.service = GeneSnapApplicationService(Path(self.temp_dir.name))

    def command(self, construct_lines=("FL", "1-80aa")):
        vector, protocol = vector_and_protocol()
        return NewExpressionProjectCommand(
            project_id="OE-APP-001",
            gene_symbol="TP53",
            species="human",
            source_cds="ATG" + "GCT" * 120 + "TAA",
            construct_lines=construct_lines,
            received_date=date(2026, 7, 13),
            due_date=date(2026, 7, 24),
            actor="tester",
            vector=vector,
            protocol=protocol,
            clones_per_construct=5,
            primer_vendor_name="标准",
            sequencing_vendor_name="标准",
        )

    def test_create_expression_project_saves_multi_construct_outputs_and_summary(self):
        stored = self.service.create_expression_project(self.command(), created_at=NOW)

        self.assertEqual(stored.snapshot.status, "design_completed")
        self.assertTrue(stored.project_folder.name.endswith("-OE"))
        self.assertTrue((stored.project_folder / "03_sequencing").is_dir())
        artifacts = self.service.expression_repository.list_artifacts(stored.project_id)
        self.assertEqual(len(artifacts), 6)
        self.assertTrue(all(item.path.exists() for item in artifacts))
        self.assertEqual(len(stored.design.constructs), 2)
        self.assertEqual(stored.snapshot.clones_per_construct, 5)

        summary = next(
            item for item in self.service.list_all_projects() if item.project_id == stored.project_id
        )
        self.assertEqual(summary.project_category, "表达类")
        self.assertEqual(summary.workflow_type, "expression")
        self.assertEqual(summary.design_summary, "表达类：FL + 1-80aa，共 2 个构建")

        reopened = GeneSnapApplicationService(Path(self.temp_dir.name))
        restored = reopened.load_expression_project(stored.project_id)
        self.assertEqual(restored.design, stored.design)
        self.assertEqual(restored.vector_design, stored.vector_design)

    def test_expression_protocol_profile_is_available_after_restart(self):
        vector, protocol = vector_and_protocol()

        saved = self.service.save_expression_profile(vector, protocol)
        reopened = GeneSnapApplicationService(Path(self.temp_dir.name))

        self.assertEqual(reopened.list_expression_profiles(), (saved,))
        self.assertEqual(reopened.load_expression_profile(saved.profile_id), (vector, protocol))

    def test_mutation_requires_reason_before_formal_export(self):
        command = self.command(("A2G",))

        with self.assertRaises(DesignConfirmationRequired):
            self.service.create_expression_project(command, created_at=NOW)
        self.assertFalse(self.service.projects_root.exists())

        confirmed = replace(command, design_confirmation_reason="已核对突变位点和最终序列")
        stored = self.service.create_expression_project(confirmed, created_at=NOW)

        self.assertFalse(stored.design.requires_confirmation)
        self.assertEqual(len(stored.design.confirmation_history), 1)
        self.assertEqual(stored.design.confirmation_history[0].actor, "tester")

    def test_expression_sequencing_scans_multiple_files_and_preserves_warnings(self):
        stored = self.service.create_expression_project(
            self.command(("FL",)),
            created_at=NOW,
        )
        construct = stored.design.constructs[0]
        sequencing_dir = stored.project_folder / "03_sequencing"
        exact_read = "T" * 80 + construct.insert_sequence + "C" * 80
        (sequencing_dir / "vendor_TP53-FL-1_run1.seq").write_text(
            exact_read,
            encoding="ascii",
        )
        (sequencing_dir / "vendor_TP53-FL-1_run2.seq").write_text(
            exact_read,
            encoding="ascii",
        )
        changed = list(construct.insert_sequence)
        changed[40] = "A" if changed[40] != "A" else "C"
        (sequencing_dir / "vendor_TP53-FL-2_result.seq").write_text(
            "T" * 80 + "".join(changed) + "C" * 80,
            encoding="ascii",
        )

        outcome = self.service.analyze_expression_sequencing(
            stored.project_id,
            actor="tester",
            analyzed_at=NOW,
        )

        latest = {item.clone_name: item for item in outcome.project.snapshot.clone_results}
        self.assertEqual(latest["TP53-FL-1"].status, "pass")
        self.assertEqual(len(latest["TP53-FL-1"].source_files), 2)
        self.assertEqual(latest["TP53-FL-2"].status, "warning")
        self.assertEqual(latest["TP53-FL-2"].substitution_count, 1)
        self.assertEqual(latest["TP53-FL-3"].status, "warning")
        self.assertEqual(outcome.project.snapshot.status, "analysis_completed")
        self.assertTrue(outcome.analysis_report.exists())

        reviewed = self.service.confirm_expression_clone_review(
            stored.project_id,
            clone_name="TP53-FL-2",
            usable=True,
            note="BLAST 和蛋白翻译复核后确认是可接受的转录本差异",
            actor="tester",
            occurred_at=NOW,
        )
        latest = {item.clone_name: item for item in reviewed.snapshot.clone_results}
        self.assertTrue(latest["TP53-FL-2"].manually_confirmed_usable)
        self.assertEqual(latest["TP53-FL-2"].manual_review_status, "usable")
        self.assertEqual(len(reviewed.snapshot.clone_results), 6)
        summary = next(
            item for item in self.service.list_all_projects() if item.project_id == stored.project_id
        )
        self.assertEqual(summary.usable_clone_names, ("TP53-FL-1", "TP53-FL-2"))


if __name__ == "__main__":
    unittest.main()
