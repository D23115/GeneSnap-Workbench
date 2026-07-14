import tempfile
import unittest
from dataclasses import replace
from datetime import date, datetime, timezone
from pathlib import Path

from genesnap_workbench.project_workflow.project_folders import (
    create_project_folder,
)
from genesnap_workbench.storage.syn_repository import (
    DuplicateProjectError,
    SQLiteSYNProjectRepository,
    StorageRevisionConflict,
)
from genesnap_workbench.template_engine.syn_exports import GeneratedArtifact
from tests.test_syn_models import make_design_version
from tests.test_syn_sequencing_confirmation import make_waiting_snapshot


NOW = datetime(2026, 7, 12, 12, 0, tzinfo=timezone.utc)


class SYNProjectFolderTests(unittest.TestCase):
    def test_project_folder_sanitizes_name_and_creates_standard_subfolders(self):
        with tempfile.TemporaryDirectory() as temp_dir:
            workspace = create_project_folder(
                Path(temp_dir),
                project_id="SYN:001",
                target_name="target/A*",
            )

            self.assertEqual(workspace.root.name, "SYN_001-target_A_-SYN")
            self.assertEqual(
                tuple(path.name for path in workspace.subfolders),
                (
                    "01_design",
                    "02_orders",
                    "03_sequencing",
                    "04_reports",
                    "05_experiment_records",
                    "99_archive",
                ),
            )
            self.assertTrue(all(path.is_dir() for path in workspace.subfolders))


class SYNPersistenceTests(unittest.TestCase):
    def setUp(self):
        self.temp_dir = tempfile.TemporaryDirectory()
        self.addCleanup(self.temp_dir.cleanup)
        self.root = Path(self.temp_dir.name)
        self.repository = SQLiteSYNProjectRepository(self.root / "genesnap.db")
        self.repository.initialize()

    def create_project(self):
        workspace = create_project_folder(
            self.root / "projects",
            project_id="SYN-001",
            target_name="SYN-target",
        )
        design = make_design_version()
        snapshot = make_waiting_snapshot()
        self.repository.create_project(
            project_id="SYN-001",
            target_name="SYN-target",
            received_date=date(2026, 7, 12),
            due_date=date(2026, 8, 3),
            project_folder=workspace.root,
            design=design,
            snapshot=snapshot,
            created_at=NOW,
        )
        return workspace, design, snapshot

    def test_project_id_is_unique_and_list_uses_saved_metadata(self):
        workspace, design, snapshot = self.create_project()

        with self.assertRaises(DuplicateProjectError):
            self.repository.create_project(
                project_id="SYN-001",
                target_name="another",
                received_date=date(2026, 7, 12),
                due_date=date(2026, 8, 3),
                project_folder=workspace.root,
                design=design,
                snapshot=snapshot,
                created_at=NOW,
            )

        summaries = self.repository.list_projects()
        self.assertEqual(len(summaries), 1)
        self.assertEqual(summaries[0].project_id, "SYN-001")
        self.assertEqual(summaries[0].status, snapshot.status)
        self.assertEqual(summaries[0].folder_suffix, "SYN")

    def test_design_and_snapshot_round_trip_after_repository_reopen(self):
        workspace, design, snapshot = self.create_project()
        reopened = SQLiteSYNProjectRepository(self.root / "genesnap.db")

        stored = reopened.load_project("SYN-001")

        self.assertEqual(stored.design, design)
        self.assertEqual(stored.snapshot, snapshot)
        self.assertEqual(stored.project_folder, workspace.root)
        self.assertEqual(stored.received_date, date(2026, 7, 12))
        self.assertEqual(stored.due_date, date(2026, 8, 3))

    def test_save_snapshot_uses_optimistic_revision(self):
        _, _, snapshot = self.create_project()
        updated = replace(snapshot, revision=snapshot.revision + 1, status="project_completed")
        self.repository.save_snapshot(
            "SYN-001",
            updated,
            expected_revision=snapshot.revision,
            updated_at=NOW,
        )

        with self.assertRaises(StorageRevisionConflict):
            self.repository.save_snapshot(
                "SYN-001",
                replace(updated, revision=updated.revision + 1),
                expected_revision=snapshot.revision,
                updated_at=NOW,
            )
        self.assertEqual(
            self.repository.load_project("SYN-001").snapshot,
            updated,
        )

    def test_artifact_history_is_append_only_and_reopens(self):
        workspace, _, _ = self.create_project()
        artifact_path = workspace.root / "01_design" / "design.json"
        artifact_path.write_text("{}", encoding="utf-8")
        artifact = GeneratedArtifact(
            artifact_type="design_json",
            design_version_id="design-v1",
            generated_at=NOW,
            path=artifact_path,
            content_sha256="a" * 64,
        )

        self.repository.append_artifacts("SYN-001", (artifact,))
        self.repository.append_artifacts(
            "SYN-001",
            (replace(artifact, path=artifact_path.with_name("design_2.json")),),
        )

        reopened = SQLiteSYNProjectRepository(self.root / "genesnap.db")
        artifacts = reopened.list_artifacts("SYN-001")
        self.assertEqual(len(artifacts), 2)
        self.assertEqual(artifacts[0].path, artifact_path)
        self.assertEqual(artifacts[1].path.name, "design_2.json")


if __name__ == "__main__":
    unittest.main()
