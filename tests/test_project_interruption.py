import tempfile
import unittest
from datetime import date, datetime, timezone
from pathlib import Path

from genesnap_workbench.app.application import GeneSnapApplicationService, NewExpressionProjectCommand
from tests.test_expression_vector_protocol import vector_and_protocol


class ProjectInterruptionTests(unittest.TestCase):
    def setUp(self):
        self.temp_dir = tempfile.TemporaryDirectory()
        self.addCleanup(self.temp_dir.cleanup)
        self.service = GeneSnapApplicationService(Path(self.temp_dir.name))
        vector, protocol = vector_and_protocol()
        self.stored = self.service.create_expression_project(
            NewExpressionProjectCommand(
                project_id="OE-PAUSE-001",
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
            created_at=datetime(2026, 7, 13, 9, tzinfo=timezone.utc),
        )

    def test_pause_freezes_remaining_and_resume_extends_due_date_once(self):
        paused_at = datetime(2026, 7, 13, 10, tzinfo=timezone.utc)
        paused = self.service.mark_molecular_interrupted(
            self.stored.project_id,
            workflow_type="expression",
            interruption_type="pause",
            note="等待客户确认标签",
            actor="tester",
            occurred_at=paused_at,
        )

        self.assertEqual(paused.snapshot.status, "abnormal_or_paused")
        self.assertEqual(paused.snapshot.interruption_type, "pause")
        self.assertEqual(paused.snapshot.interrupted_previous_status, "design_completed")
        self.assertEqual(paused.snapshot.frozen_remaining_workdays, 9)

        resumed = self.service.resume_molecular_project(
            paused.project_id,
            workflow_type="expression",
            actor="tester",
            occurred_at=datetime(2026, 7, 15, 10, tzinfo=timezone.utc),
            note="客户已确认，继续项目",
        )

        self.assertEqual(resumed.snapshot.status, "design_completed")
        self.assertEqual(resumed.snapshot.accumulated_paused_workdays, 2)
        self.assertEqual(resumed.snapshot.effective_due_date, date(2026, 7, 28))
        self.assertIsNone(resumed.snapshot.interruption_type)
        self.assertEqual(resumed.snapshot.interruption_history[-1].paused_workdays, 2)

        summary = next(
            item for item in self.service.list_all_projects() if item.project_id == resumed.project_id
        )
        self.assertEqual(summary.due_date, date(2026, 7, 28))

    def test_interruption_requires_type_and_reason(self):
        with self.assertRaisesRegex(ValueError, "类型"):
            self.service.mark_molecular_interrupted(
                self.stored.project_id,
                workflow_type="expression",
                interruption_type="other",
                note="reason",
                actor="tester",
                occurred_at=datetime(2026, 7, 13, tzinfo=timezone.utc),
            )
        with self.assertRaisesRegex(ValueError, "原因"):
            self.service.mark_molecular_interrupted(
                self.stored.project_id,
                workflow_type="expression",
                interruption_type="abnormal",
                note="",
                actor="tester",
                occurred_at=datetime(2026, 7, 13, tzinfo=timezone.utc),
            )


if __name__ == "__main__":
    unittest.main()
