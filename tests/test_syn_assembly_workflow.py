import unittest
from dataclasses import replace
from datetime import datetime, timezone

from genesnap_workbench.domain.syn import (
    PlasmidPrepStatus,
    SYNAssemblyAttemptRecord,
    SYNAssemblyAttemptResult,
    SYNAuditEvent,
    SYNColonyPCRResult,
    SYNProjectSnapshot,
    SYNRoute,
)
from genesnap_workbench.project_workflow.syn_service import (
    SYNRevisionConflict,
    SYNWorkflowRuleError,
    SYNWorkflowService,
)
from genesnap_workbench.project_workflow.syn_state import (
    SYNStateTransitionService,
)


NOW = datetime(2026, 7, 12, 9, 0, tzinfo=timezone.utc)


def make_snapshot(
    *,
    status: str = "syn_assembly_in_progress",
    substep: str | None = "colony_pcr",
) -> SYNProjectSnapshot:
    return SYNProjectSnapshot(
        project_id="SYN-001",
        revision=1,
        status=status,
        resuspension_data_status="complete",
        syn_assembly_round_no=0,
        syn_assembly_substep=substep,
        active_design_version_id="design-v1",
        attempts=(),
        colonies=(),
        prep_records=(),
        sequencing_confirmations=(),
        status_history=(
            SYNAuditEvent(
                event_id="event-1",
                event_type="start_assembly",
                occurred_at=NOW,
                actor="tester",
            ),
        ),
        manual_override_history=(),
    )


class SYNAssemblyWorkflowTests(unittest.TestCase):
    def setUp(self):
        self.service = SYNWorkflowService()

    def start_screening(
        self,
        snapshot: SYNProjectSnapshot | None = None,
        *,
        route: SYNRoute = SYNRoute.SINGLE_POOL,
        high_risk: bool = False,
        count: int | None = None,
    ) -> SYNProjectSnapshot:
        current = snapshot or make_snapshot()
        return self.service.start_initial_colony_screening(
            current,
            target_name="SYN-target",
            route=route,
            high_risk=high_risk,
            colony_count=count,
            expected_revision=current.revision,
            attempt_id="attempt-1",
            event_id="event-screening",
            actor="tester",
            occurred_at=NOW,
        )

    def test_default_colony_count_is_8_or_12_by_risk_route(self):
        ordinary = self.start_screening()
        high_risk = self.start_screening(high_risk=True)
        modular = self.start_screening(route=SYNRoute.MODULAR)

        self.assertEqual(len(ordinary.colonies), 8)
        self.assertEqual(len(high_risk.colonies), 12)
        self.assertEqual(len(modular.colonies), 12)

    def test_first_round_names_and_custom_count_are_deterministic(self):
        updated = self.start_screening(count=3)

        self.assertEqual(updated.syn_assembly_round_no, 1)
        self.assertEqual(
            tuple(colony.display_name for colony in updated.colonies),
            ("SYN-target-C01", "SYN-target-C02", "SYN-target-C03"),
        )
        self.assertEqual(updated.attempts[0].restart_from_substep, "assembly_pcr")

    def test_optimistic_revision_blocks_duplicate_action(self):
        with self.assertRaises(SYNRevisionConflict):
            self.service.start_initial_colony_screening(
                make_snapshot(),
                target_name="SYN-target",
                route=SYNRoute.SINGLE_POOL,
                high_risk=False,
                colony_count=8,
                expected_revision=0,
                attempt_id="attempt-1",
                event_id="event-screening",
                actor="tester",
                occurred_at=NOW,
            )

    def test_substep_experiment_records_append_without_auto_advancing(self):
        snapshot = self.service.begin_initial_attempt(
            make_snapshot(substep="assembly_pcr"),
            expected_revision=1,
            attempt_id="attempt-1",
            event_id="event-attempt",
            actor="tester",
            occurred_at=NOW,
        )

        for index, result in enumerate(
            (
                SYNAssemblyAttemptResult.FAILED,
                SYNAssemblyAttemptResult.SUCCESS,
            ),
            start=1,
        ):
            snapshot = self.service.record_assembly_step(
                snapshot,
                result=result,
                note=f"第 {index} 次记录",
                expected_revision=snapshot.revision,
                record_id=f"step-record-{index}",
                event_id=f"event-step-{index}",
                actor="tester",
                occurred_at=NOW,
            )

        self.assertEqual(snapshot.syn_assembly_substep, "assembly_pcr")
        self.assertEqual(len(snapshot.step_records), 2)
        self.assertEqual(
            tuple(record.step_attempt_no for record in snapshot.step_records),
            (1, 2),
        )
        self.assertEqual(
            tuple(record.result for record in snapshot.step_records),
            (
                SYNAssemblyAttemptResult.FAILED,
                SYNAssemblyAttemptResult.SUCCESS,
            ),
        )

    def test_colony_results_are_appended_and_only_positive_can_be_selected(self):
        snapshot = self.start_screening(count=4)
        for index, result in enumerate(
            (
                SYNColonyPCRResult.POSITIVE,
                SYNColonyPCRResult.POSITIVE,
                SYNColonyPCRResult.NEGATIVE,
                SYNColonyPCRResult.UNCERTAIN,
            ),
            start=1,
        ):
            snapshot = self.service.record_colony_pcr(
                snapshot,
                clone_id=f"attempt-1-clone-{index}",
                result=result,
                observed_note=None,
                expected_revision=snapshot.revision,
                event_id=f"event-result-{index}",
                actor="tester",
                occurred_at=NOW,
            )

        self.assertEqual(len(snapshot.colonies), 8)
        with self.assertRaisesRegex(SYNWorkflowRuleError, "阳性"):
            self.service.select_clones_for_prep(
                snapshot,
                clone_ids=("attempt-1-clone-3",),
                expected_revision=snapshot.revision,
                event_id="event-select-bad",
                actor="tester",
                occurred_at=NOW,
            )

        selected = self.service.select_clones_for_prep(
            snapshot,
            clone_ids=None,
            expected_revision=snapshot.revision,
            event_id="event-select",
            actor="tester",
            occurred_at=NOW,
        )
        self.assertEqual(selected.status, "plasmid_prep_in_progress")
        self.assertEqual(
            tuple(record.clone_id for record in selected.prep_records),
            ("attempt-1-clone-1", "attempt-1-clone-2"),
        )

    def test_all_selected_preps_must_complete_before_waiting_for_sequencing(self):
        snapshot = self.start_screening(count=3)
        for index in range(1, 4):
            snapshot = self.service.record_colony_pcr(
                snapshot,
                clone_id=f"attempt-1-clone-{index}",
                result=SYNColonyPCRResult.POSITIVE,
                observed_note=None,
                expected_revision=snapshot.revision,
                event_id=f"event-result-{index}",
                actor="tester",
                occurred_at=NOW,
            )
        snapshot = self.service.select_clones_for_prep(
            snapshot,
            clone_ids=None,
            expected_revision=snapshot.revision,
            event_id="event-select",
            actor="tester",
            occurred_at=NOW,
        )

        for index in range(1, 3):
            snapshot = self.service.record_plasmid_prep(
                snapshot,
                clone_id=f"attempt-1-clone-{index}",
                status=PlasmidPrepStatus.COMPLETED,
                expected_revision=snapshot.revision,
                event_id=f"event-prep-{index}",
                actor="tester",
                occurred_at=NOW,
            )
            self.assertEqual(snapshot.status, "plasmid_prep_in_progress")

        snapshot = self.service.record_plasmid_prep(
            snapshot,
            clone_id="attempt-1-clone-3",
            status=PlasmidPrepStatus.COMPLETED,
            expected_revision=snapshot.revision,
            event_id="event-prep-3",
            actor="tester",
            occurred_at=NOW,
        )
        self.assertEqual(snapshot.status, "awaiting_sequencing_confirmation")

    def test_additional_screening_preview_is_side_effect_free_and_continues_numbers(self):
        snapshot = self.start_screening(count=3)
        snapshot = replace(
            snapshot,
            status="awaiting_sequencing_confirmation",
        )

        preview = self.service.preview_additional_screening(
            snapshot,
            target_name="SYN-target",
            expected_revision=snapshot.revision,
            preview_id="preview-1",
        )

        self.assertEqual(len(preview.display_names), 8)
        self.assertEqual(preview.display_names[0], "SYN-target-C04")
        self.assertEqual(preview.display_names[-1], "SYN-target-C11")
        self.assertEqual(len(snapshot.colonies), 3)

        confirmed = self.service.confirm_additional_screening(
            snapshot,
            preview,
            expected_revision=snapshot.revision,
            event_id="event-additional",
            actor="tester",
            occurred_at=NOW,
        )
        self.assertEqual(confirmed.status, "syn_assembly_in_progress")
        self.assertEqual(confirmed.syn_assembly_substep, "colony_pcr")
        self.assertEqual(len(confirmed.colonies), 11)

    def test_additional_screening_does_not_reselect_previously_prepped_clone(self):
        snapshot = self.start_screening(count=2)
        snapshot = self.service.record_colony_pcr(
            snapshot,
            clone_id="attempt-1-clone-1",
            result=SYNColonyPCRResult.POSITIVE,
            observed_note=None,
            expected_revision=snapshot.revision,
            event_id="event-result-old",
            actor="tester",
            occurred_at=NOW,
        )
        snapshot = self.service.select_clones_for_prep(
            snapshot,
            clone_ids=None,
            expected_revision=snapshot.revision,
            event_id="event-select-old",
            actor="tester",
            occurred_at=NOW,
        )
        snapshot = self.service.record_plasmid_prep(
            snapshot,
            clone_id="attempt-1-clone-1",
            status=PlasmidPrepStatus.COMPLETED,
            expected_revision=snapshot.revision,
            event_id="event-prep-old",
            actor="tester",
            occurred_at=NOW,
        )
        preview = self.service.preview_additional_screening(
            snapshot,
            target_name="SYN-target",
            colony_count=1,
            expected_revision=snapshot.revision,
            preview_id="preview-new",
        )
        snapshot = self.service.confirm_additional_screening(
            snapshot,
            preview,
            expected_revision=snapshot.revision,
            event_id="event-additional",
            actor="tester",
            occurred_at=NOW,
        )
        snapshot = self.service.record_colony_pcr(
            snapshot,
            clone_id="attempt-1-clone-3",
            result=SYNColonyPCRResult.POSITIVE,
            observed_note=None,
            expected_revision=snapshot.revision,
            event_id="event-result-new",
            actor="tester",
            occurred_at=NOW,
        )

        selected = self.service.select_clones_for_prep(
            snapshot,
            clone_ids=None,
            expected_revision=snapshot.revision,
            event_id="event-select-new",
            actor="tester",
            occurred_at=NOW,
        )
        self.assertEqual(selected.prep_records[-1].clone_id, "attempt-1-clone-3")

    def test_restart_creates_new_round_and_resets_clone_numbers(self):
        snapshot = self.start_screening(count=2)
        snapshot = replace(
            snapshot,
            status="awaiting_sequencing_confirmation",
        )
        restarted = self.service.restart_assembly(
            snapshot,
            restart_from_substep="vector_assembly_transformation",
            expected_revision=snapshot.revision,
            attempt_id="attempt-2",
            event_id="event-restart",
            actor="tester",
            occurred_at=NOW,
        )
        self.assertEqual(restarted.syn_assembly_round_no, 2)
        self.assertEqual(
            restarted.syn_assembly_substep,
            "vector_assembly_transformation",
        )

        with self.assertRaisesRegex(SYNWorkflowRuleError, "菌落 PCR"):
            self.service.create_colonies_for_active_round(
                restarted,
                target_name="SYN-target",
                colony_count=2,
                expected_revision=restarted.revision,
                event_id="event-too-early",
                actor="tester",
                occurred_at=NOW,
            )
        restarted = SYNStateTransitionService().advance_assembly_substep(
            restarted,
            to_substep="colony_pcr",
            event_id="event-to-colony-pcr",
            actor="tester",
            occurred_at=NOW,
        )

        screened = self.service.create_colonies_for_active_round(
            restarted,
            target_name="SYN-target",
            colony_count=2,
            expected_revision=restarted.revision,
            event_id="event-round-2-colonies",
            actor="tester",
            occurred_at=NOW,
        )
        self.assertEqual(
            tuple(colony.display_name for colony in screened.colonies[-2:]),
            ("SYN-target-R2-C01", "SYN-target-R2-C02"),
        )


if __name__ == "__main__":
    unittest.main()
