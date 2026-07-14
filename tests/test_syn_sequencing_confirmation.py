import unittest
from datetime import datetime, timezone

from genesnap_workbench.domain.syn import (
    PlasmidPrepRecord,
    PlasmidPrepStatus,
    SYNAuditEvent,
    SYNAssemblyAttemptRecord,
    SYNAssemblyAttemptResult,
    SYNColonyPCRRecord,
    SYNColonyPCRResult,
    SYNProjectSnapshot,
    SYNSequencingResult,
)
from genesnap_workbench.project_workflow.syn_service import (
    SYNWorkflowRuleError,
    SYNWorkflowService,
)


NOW = datetime(2026, 7, 12, 10, 0, tzinfo=timezone.utc)


def make_waiting_snapshot(
    *,
    prep_statuses: tuple[PlasmidPrepStatus, ...] = (
        PlasmidPrepStatus.COMPLETED,
        PlasmidPrepStatus.COMPLETED,
        PlasmidPrepStatus.COMPLETED,
        PlasmidPrepStatus.COMPLETED,
    ),
) -> SYNProjectSnapshot:
    attempt = SYNAssemblyAttemptRecord(
        attempt_id="attempt-1",
        project_id="SYN-001",
        design_version_id="design-v1",
        syn_assembly_round_no=1,
        restart_from_substep="assembly_pcr",
        result=SYNAssemblyAttemptResult.SUCCESS,
        started_at=NOW,
        completed_at=NOW,
    )
    colonies = tuple(
        SYNColonyPCRRecord(
            clone_id=f"clone-{index}",
            attempt_id=attempt.attempt_id,
            clone_no=index,
            display_name=f"SYN-target-C{index:02d}",
            result=SYNColonyPCRResult.POSITIVE,
            recorded_at=NOW,
        )
        for index in range(1, len(prep_statuses) + 1)
    )
    prep_records = tuple(
        PlasmidPrepRecord(
            clone_id=f"clone-{index}",
            selected_for_prep=True,
            status=status,
            completed_at=(NOW if status is PlasmidPrepStatus.COMPLETED else None),
        )
        for index, status in enumerate(prep_statuses, start=1)
    )
    return SYNProjectSnapshot(
        project_id="SYN-001",
        revision=1,
        status="awaiting_sequencing_confirmation",
        resuspension_data_status="complete",
        syn_assembly_round_no=1,
        syn_assembly_substep=None,
        active_design_version_id="design-v1",
        attempts=(attempt,),
        colonies=colonies,
        prep_records=prep_records,
        sequencing_confirmations=(),
        status_history=(
            SYNAuditEvent(
                event_id="event-ready",
                event_type="finish_plasmid_prep",
                occurred_at=NOW,
                actor="tester",
            ),
        ),
        manual_override_history=(),
    )


class SYNSequencingConfirmationTests(unittest.TestCase):
    def setUp(self):
        self.service = SYNWorkflowService()

    def confirm(
        self,
        snapshot: SYNProjectSnapshot,
        clone_id: str,
        result: SYNSequencingResult,
        index: int,
    ) -> SYNProjectSnapshot:
        return self.service.confirm_sequencing(
            snapshot,
            clone_id=clone_id,
            result=result,
            note=None,
            expected_revision=snapshot.revision,
            confirmation_id=f"confirmation-{index}",
            event_id=f"event-confirmation-{index}",
            actor="tester",
            occurred_at=NOW,
        )

    def test_only_completed_prep_clone_can_be_confirmed(self):
        snapshot = make_waiting_snapshot(
            prep_statuses=(PlasmidPrepStatus.IN_PROGRESS,),
        )

        with self.assertRaisesRegex(SYNWorkflowRuleError, "完成小提"):
            self.confirm(
                snapshot,
                "clone-1",
                SYNSequencingResult.CORRECT,
                1,
            )

    def test_new_confirmation_supersedes_latest_without_deleting_history(self):
        snapshot = make_waiting_snapshot()
        snapshot = self.confirm(
            snapshot,
            "clone-1",
            SYNSequencingResult.INCORRECT,
            1,
        )
        snapshot = self.confirm(
            snapshot,
            "clone-1",
            SYNSequencingResult.CORRECT,
            2,
        )

        self.assertEqual(len(snapshot.sequencing_confirmations), 2)
        self.assertEqual(
            snapshot.sequencing_confirmations[-1].supersedes_confirmation_id,
            "confirmation-1",
        )
        summary = self.service.get_syn_sequencing_summary(snapshot)
        self.assertEqual(summary.correct_count, 1)
        self.assertEqual(summary.incorrect_count, 0)
        self.assertEqual(summary.unconfirmed_count, 3)
        self.assertEqual(summary.usable_clone_names, ("SYN-target-C01",))

    def test_summary_separates_correct_incorrect_uncertain_and_unconfirmed(self):
        snapshot = make_waiting_snapshot()
        snapshot = self.confirm(
            snapshot,
            "clone-1",
            SYNSequencingResult.CORRECT,
            1,
        )
        snapshot = self.confirm(
            snapshot,
            "clone-2",
            SYNSequencingResult.INCORRECT,
            2,
        )
        snapshot = self.confirm(
            snapshot,
            "clone-3",
            SYNSequencingResult.UNCERTAIN,
            3,
        )

        summary = self.service.get_syn_sequencing_summary(snapshot)
        self.assertEqual(summary.correct_count, 1)
        self.assertEqual(summary.incorrect_count, 1)
        self.assertEqual(summary.uncertain_count, 1)
        self.assertEqual(summary.unconfirmed_count, 1)
        self.assertEqual(
            summary.display_text,
            "人工确认：正确 1 / 错误 1 / 不确定 1 / 未确认 1",
        )

    def test_completion_requires_correct_clone_and_explicit_user_confirmation(self):
        snapshot = self.confirm(
            make_waiting_snapshot(),
            "clone-1",
            SYNSequencingResult.CORRECT,
            1,
        )

        with self.assertRaisesRegex(SYNWorkflowRuleError, "确认完成"):
            self.service.complete_project(
                snapshot,
                user_confirmed=False,
                expected_revision=snapshot.revision,
                event_id="event-complete",
                actor="tester",
                occurred_at=NOW,
            )

        completed = self.service.complete_project(
            snapshot,
            user_confirmed=True,
            expected_revision=snapshot.revision,
            event_id="event-complete",
            actor="tester",
            occurred_at=NOW,
        )
        self.assertEqual(completed.status, "project_completed")
        self.assertEqual(completed.actual_completed_at, NOW)
        self.assertEqual(
            completed.status_history[-1].event_type,
            "complete_project",
        )

    def test_no_correct_clone_cannot_complete_but_can_continue_screening(self):
        snapshot = self.confirm(
            make_waiting_snapshot(),
            "clone-1",
            SYNSequencingResult.INCORRECT,
            1,
        )

        with self.assertRaisesRegex(SYNWorkflowRuleError, "正确克隆"):
            self.service.complete_project(
                snapshot,
                user_confirmed=True,
                expected_revision=snapshot.revision,
                event_id="event-complete",
                actor="tester",
                occurred_at=NOW,
            )

        preview = self.service.preview_additional_screening(
            snapshot,
            target_name="SYN-target",
            colony_count=2,
            expected_revision=snapshot.revision,
            preview_id="preview-1",
        )
        self.assertEqual(preview.display_names[-1], "SYN-target-C06")

    def test_correct_clone_blocks_additional_screening_and_restart(self):
        snapshot = self.confirm(
            make_waiting_snapshot(),
            "clone-1",
            SYNSequencingResult.CORRECT,
            1,
        )

        with self.assertRaisesRegex(SYNWorkflowRuleError, "已有正确克隆"):
            self.service.preview_additional_screening(
                snapshot,
                target_name="SYN-target",
                expected_revision=snapshot.revision,
                preview_id="preview-1",
            )
        with self.assertRaisesRegex(SYNWorkflowRuleError, "已有正确克隆"):
            self.service.restart_assembly(
                snapshot,
                restart_from_substep="assembly_pcr",
                expected_revision=snapshot.revision,
                attempt_id="attempt-2",
                event_id="event-restart",
                actor="tester",
                occurred_at=NOW,
            )


if __name__ == "__main__":
    unittest.main()
