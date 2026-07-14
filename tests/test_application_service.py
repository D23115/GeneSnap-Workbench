import tempfile
import unittest
from datetime import date, datetime, timezone
from pathlib import Path

from genesnap_workbench.app.application import (
    DesignConfirmationRequired,
    GeneSnapApplicationService,
    NewSYNProjectCommand,
)
from tests.test_syn_design_engine import artificial_sequence


NOW = datetime(2026, 7, 12, 14, 0, tzinfo=timezone.utc)


class ApplicationServiceTests(unittest.TestCase):
    def setUp(self):
        self.temp_dir = tempfile.TemporaryDirectory()
        self.addCleanup(self.temp_dir.cleanup)
        self.data_root = Path(self.temp_dir.name)
        self.service = GeneSnapApplicationService(self.data_root)

    def command(self):
        return NewSYNProjectCommand(
            project_id="SYN-APP-001",
            target_name="artificial-600",
            raw_sequence=artificial_sequence(600, seed=7),
            input_format="plain",
            linearization_site="EcoRV",
            received_date=date(2026, 7, 12),
            due_date=date(2026, 8, 3),
            actor="tester",
            vector_sequence_confirmed=True,
        )

    def test_prepare_design_has_no_database_or_project_folder_side_effect(self):
        prepared = self.service.prepare_syn_project(self.command(), created_at=NOW)

        self.assertEqual(prepared.design.project_id, "SYN-APP-001")
        self.assertTrue(prepared.design.requires_confirmation)
        self.assertEqual(self.service.list_projects(), ())
        self.assertFalse((self.data_root / "projects").exists())

    def test_save_requires_reason_when_design_has_confirmation_items(self):
        prepared = self.service.prepare_syn_project(self.command(), created_at=NOW)

        with self.assertRaises(DesignConfirmationRequired):
            self.service.save_prepared_syn_project(
                self.command(),
                prepared,
                design_confirmation_reason=None,
                created_at=NOW,
            )
        self.assertEqual(self.service.list_projects(), ())

    def test_saved_project_reopens_with_exports_and_audit_records(self):
        command = self.command()
        prepared = self.service.prepare_syn_project(command, created_at=NOW)
        stored = self.service.save_prepared_syn_project(
            command,
            prepared,
            design_confirmation_reason="已核对模块和 overlap 风险",
            created_at=NOW,
        )

        self.assertEqual(stored.snapshot.status, "design_completed")
        self.assertTrue(stored.project_folder.is_dir())
        self.assertTrue((stored.project_folder / "03_sequencing").is_dir())
        artifacts = self.service.repository.list_artifacts(command.project_id)
        self.assertGreaterEqual(len(artifacts), 5)
        self.assertTrue(all(item.path.exists() for item in artifacts))
        self.assertTrue(
            any(
                item.field_path == "vector_reference_confirmation"
                for item in stored.design.manual_overrides
            ),
        )

        reopened = GeneSnapApplicationService(self.data_root)
        restored = reopened.load_project(command.project_id)
        self.assertEqual(restored.design, stored.design)
        self.assertEqual(restored.snapshot, stored.snapshot)


if __name__ == "__main__":
    unittest.main()
