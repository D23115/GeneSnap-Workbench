import tempfile
import unittest
from datetime import date, datetime, timezone
from pathlib import Path

from genesnap_workbench.app.application import (
    GeneSnapApplicationService,
    NewExpressionProjectCommand,
)
from tests.test_expression_vector_protocol import vector_and_protocol


NOW = datetime(2026, 7, 13, 4, 0, tzinfo=timezone.utc)


class SequencingSubmissionTests(unittest.TestCase):
    def setUp(self):
        self.temp_dir = tempfile.TemporaryDirectory()
        self.addCleanup(self.temp_dir.cleanup)
        self.service = GeneSnapApplicationService(Path(self.temp_dir.name))
        vector, protocol = vector_and_protocol()
        self.stored = self.service.create_expression_project(
            NewExpressionProjectCommand(
                project_id="OE-SUB-001",
                gene_symbol="TP53",
                species="human",
                source_cds="ATG" + "GCT" * 120 + "TAA",
                construct_lines=("FL",),
                received_date=date(2026, 7, 13),
                due_date=date(2026, 7, 24),
                actor="tester",
                vector=vector,
                protocol=protocol,
                clones_per_construct=5,
            ),
            created_at=NOW,
        )
        for action in ("mark_primers_ordered", "mark_primers_arrived", "start_cloning"):
            self.stored = self.service.transition_molecular_project(
                self.stored.project_id,
                workflow_type="expression",
                action=action,
                actor="tester",
                occurred_at=NOW,
            )

    def test_mark_sent_creates_traceable_submission_with_all_expected_clones(self):
        sent = self.service.transition_molecular_project(
            self.stored.project_id,
            workflow_type="expression",
            action="mark_sent_for_sequencing",
            actor="tester",
            occurred_at=NOW,
            internal_submission_no="SEQ-20260713-01",
            vendor_order_no="ORDER-123456",
        )

        self.assertEqual(sent.snapshot.status, "sequencing_in_progress")
        self.assertEqual(len(sent.snapshot.sequencing_submissions), 1)
        submission = sent.snapshot.sequencing_submissions[0]
        self.assertEqual(submission.round_no, 1)
        self.assertEqual(submission.submission_kind, "initial")
        self.assertEqual(submission.internal_submission_no, "SEQ-20260713-01")
        self.assertEqual(submission.vendor_order_no, "ORDER-123456")
        self.assertEqual(
            submission.sample_names,
            tuple(f"TP53-FL-{index}" for index in range(1, 6)),
        )

        summary = next(
            item for item in self.service.list_all_projects() if item.project_id == sent.project_id
        )
        self.assertEqual(summary.latest_internal_submission_no, "SEQ-20260713-01")
        self.assertEqual(summary.latest_vendor_order_no, "ORDER-123456")

    def test_latest_tracking_numbers_can_be_corrected_without_losing_original_event(self):
        sent = self.service.transition_molecular_project(
            self.stored.project_id,
            workflow_type="expression",
            action="mark_sent_for_sequencing",
            actor="tester",
            occurred_at=NOW,
        )
        original_history = sent.snapshot.status_history

        updated = self.service.update_latest_sequencing_tracking(
            sent.project_id,
            workflow_type="expression",
            internal_submission_no="SEQ-CORRECTED",
            vendor_order_no="ORDER-CORRECTED",
            actor="tester",
            occurred_at=NOW,
            note="供应商返回订单号后补录",
        )

        latest = updated.snapshot.sequencing_submissions[-1]
        self.assertEqual(latest.internal_submission_no, "SEQ-CORRECTED")
        self.assertEqual(latest.vendor_order_no, "ORDER-CORRECTED")
        self.assertEqual(updated.snapshot.status_history[: len(original_history)], original_history)
        self.assertEqual(updated.snapshot.status_history[-1].event_type, "update_sequencing_tracking")

    def test_analysis_marks_latest_submission_analyzed(self):
        sent = self.service.transition_molecular_project(
            self.stored.project_id,
            workflow_type="expression",
            action="mark_sent_for_sequencing",
            actor="tester",
            occurred_at=NOW,
        )
        construct = sent.design.constructs[0]
        (sent.project_folder / "03_sequencing" / "TP53-FL-1_result.seq").write_text(
            "A" * 50 + construct.insert_sequence + "C" * 50,
            encoding="ascii",
        )

        analyzed = self.service.analyze_expression_sequencing(
            sent.project_id,
            actor="tester",
            analyzed_at=NOW,
        ).project

        self.assertEqual(analyzed.snapshot.sequencing_submissions[-1].status, "analyzed")


if __name__ == "__main__":
    unittest.main()
