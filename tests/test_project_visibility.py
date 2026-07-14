import tempfile
import unittest
from datetime import date, datetime, timezone
from pathlib import Path

from genesnap_workbench.app.application import GeneSnapApplicationService, NewExpressionProjectCommand
from tests.test_expression_vector_protocol import vector_and_protocol


class ProjectVisibilityTests(unittest.TestCase):
    def test_manual_hidden_is_independent_from_completed_archive(self):
        with tempfile.TemporaryDirectory() as temp_dir:
            service = GeneSnapApplicationService(Path(temp_dir))
            vector, protocol = vector_and_protocol()
            stored = service.create_expression_project(
                NewExpressionProjectCommand(
                    project_id="OE-HIDE-001",
                    gene_symbol="TP53",
                    species="human",
                    source_cds="ATG" + "GCT" * 120 + "TAA",
                    construct_lines=("FL",),
                    received_date=date(2026, 7, 13),
                    due_date=date(2026, 7, 24),
                    actor="tester",
                    vector=vector,
                    protocol=protocol,
                ),
                created_at=datetime(2026, 7, 13, tzinfo=timezone.utc),
            )

            service.set_project_hidden(
                stored.project_id,
                hidden=True,
                reason="项目取消，暂不显示",
                actor="tester",
                occurred_at=datetime(2026, 7, 14, tzinfo=timezone.utc),
            )
            hidden = next(item for item in service.list_all_projects() if item.project_id == stored.project_id)
            self.assertTrue(hidden.is_manually_hidden)
            self.assertNotEqual(hidden.status, "project_completed")

            service.set_project_hidden(
                stored.project_id,
                hidden=False,
                reason="恢复项目显示",
                actor="tester",
                occurred_at=datetime(2026, 7, 15, tzinfo=timezone.utc),
            )
            visible = next(item for item in service.list_all_projects() if item.project_id == stored.project_id)
            self.assertFalse(visible.is_manually_hidden)


if __name__ == "__main__":
    unittest.main()
