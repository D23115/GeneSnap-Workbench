import unittest
from datetime import datetime, timezone

from genesnap_workbench.domain.syn import SYNProjectSnapshot, SYNAuditEvent
from genesnap_workbench.project_workflow.syn_materials import MaterialReadiness
from genesnap_workbench.project_workflow.syn_state import (
    SYNMaterialOverrideRequired,
    SYNStateTransitionError,
    SYNStateTransitionService,
    display_status_label,
)


NOW = datetime(2026, 7, 12, 8, 0, tzinfo=timezone.utc)


def make_snapshot(status: str = "materials_ordered") -> SYNProjectSnapshot:
    return SYNProjectSnapshot(
        project_id="SYN-001",
        revision=1,
        status=status,
        resuspension_data_status="missing",
        syn_assembly_round_no=0,
        syn_assembly_substep=None,
        active_design_version_id="design-v1",
        attempts=(),
        colonies=(),
        prep_records=(),
        sequencing_confirmations=(),
        status_history=(
            SYNAuditEvent(
                event_id="event-1",
                event_type="mark_materials_ordered",
                occurred_at=NOW,
                actor="tester",
            ),
        ),
        manual_override_history=(),
    )


def complete_readiness() -> MaterialReadiness:
    return MaterialReadiness(
        is_ready=True,
        can_start_with_override=True,
        missing_oligo_ids=(),
        missing_fields=(),
        errors=(),
    )


class SYNStateTransitionTests(unittest.TestCase):
    def setUp(self):
        self.service = SYNStateTransitionService()

    def test_material_status_labels_are_workflow_aware(self):
        self.assertEqual(
            display_status_label("materials_ordered", "de_novo_gene_synthesis"),
            "oligo 已订购",
        )
        self.assertEqual(
            display_status_label("materials_arrived", "de_novo_gene_synthesis"),
            "oligo 已到货",
        )
        self.assertEqual(
            display_status_label("materials_ordered", "expression"),
            "引物已订购",
        )

    def test_mark_materials_arrived_updates_status_and_appends_audit(self):
        snapshot = make_snapshot()

        updated = self.service.mark_materials_arrived(
            snapshot,
            resuspension_complete=False,
            event_id="event-2",
            actor="tester",
            occurred_at=NOW,
        )

        self.assertEqual(updated.status, "materials_arrived")
        self.assertEqual(updated.resuspension_data_status, "missing")
        self.assertEqual(
            tuple(event.event_id for event in updated.status_history),
            ("event-1", "event-2"),
        )
        event = updated.status_history[-1]
        self.assertEqual(event.from_status, "materials_ordered")
        self.assertEqual(event.to_status, "materials_arrived")

    def test_mark_materials_ordered_uses_controlled_transition(self):
        snapshot = make_snapshot("design_completed")

        updated = self.service.mark_materials_ordered(
            snapshot,
            event_id="event-ordered",
            actor="tester",
            occurred_at=NOW,
        )

        self.assertEqual(updated.status, "materials_ordered")
        self.assertEqual(updated.status_history[-1].event_type, "mark_materials_ordered")

    def test_missing_materials_require_confirmation_and_reason(self):
        snapshot = make_snapshot("materials_arrived")
        readiness = MaterialReadiness(
            is_ready=False,
            can_start_with_override=True,
            missing_oligo_ids=("oligo-2",),
            missing_fields=("standard_volume_per_oligo_ul",),
            errors=(),
        )

        with self.assertRaises(SYNMaterialOverrideRequired):
            self.service.start_assembly(
                snapshot,
                readiness,
                confirm_missing=False,
                override_reason=None,
                event_id="event-2",
                override_id="override-1",
                actor="tester",
                occurred_at=NOW,
            )
        with self.assertRaisesRegex(ValueError, "reason"):
            self.service.start_assembly(
                snapshot,
                readiness,
                confirm_missing=True,
                override_reason=" ",
                event_id="event-2",
                override_id="override-1",
                actor="tester",
                occurred_at=NOW,
            )

        updated = self.service.start_assembly(
            snapshot,
            readiness,
            confirm_missing=True,
            override_reason="实际已按纸质记录完成复溶",
            event_id="event-2",
            override_id="override-1",
            actor="tester",
            occurred_at=NOW,
        )

        self.assertEqual(updated.status, "syn_assembly_in_progress")
        self.assertEqual(updated.syn_assembly_substep, "assembly_pcr")
        self.assertEqual(len(updated.manual_override_history), 1)
        self.assertEqual(
            updated.manual_override_history[0].reason,
            "实际已按纸质记录完成复溶",
        )

    def test_blocking_material_error_cannot_be_overridden(self):
        readiness = MaterialReadiness(
            is_ready=False,
            can_start_with_override=False,
            missing_oligo_ids=(),
            missing_fields=(),
            errors=("设计版本不一致",),
        )

        with self.assertRaisesRegex(SYNStateTransitionError, "设计版本不一致"):
            self.service.start_assembly(
                make_snapshot("materials_arrived"),
                readiness,
                confirm_missing=True,
                override_reason="尝试覆盖",
                event_id="event-2",
                override_id="override-1",
                actor="tester",
                occurred_at=NOW,
            )

    def test_assembly_substeps_only_advance_in_controlled_order(self):
        snapshot = self.service.start_assembly(
            make_snapshot("materials_arrived"),
            complete_readiness(),
            confirm_missing=False,
            override_reason=None,
            event_id="event-2",
            override_id="override-unused",
            actor="tester",
            occurred_at=NOW,
        )

        with self.assertRaises(SYNStateTransitionError):
            self.service.advance_assembly_substep(
                snapshot,
                to_substep="colony_pcr",
                event_id="event-jump",
                actor="tester",
                occurred_at=NOW,
            )

        for index, substep in enumerate(
            (
                "amplification_pcr",
                "vector_assembly_transformation",
                "colony_pcr",
            ),
            start=3,
        ):
            snapshot = self.service.advance_assembly_substep(
                snapshot,
                to_substep=substep,
                event_id=f"event-{index}",
                actor="tester",
                occurred_at=NOW,
            )
        self.assertEqual(snapshot.syn_assembly_substep, "colony_pcr")

    def test_status_correction_preserves_original_history(self):
        snapshot = make_snapshot("materials_arrived")

        corrected = self.service.correct_status(
            snapshot,
            to_status="design_completed",
            reason="误点了到货",
            event_id="event-correction",
            actor="tester",
            occurred_at=NOW,
        )

        self.assertEqual(corrected.status, "design_completed")
        self.assertEqual(tuple(event.event_id for event in snapshot.status_history), ("event-1",))
        self.assertEqual(
            tuple(event.event_id for event in corrected.status_history),
            ("event-1", "event-correction"),
        )
        correction = corrected.status_history[-1]
        self.assertEqual(correction.event_type, "correct_status")
        self.assertEqual(correction.from_status, "materials_arrived")
        self.assertEqual(correction.to_status, "design_completed")
        self.assertEqual(correction.note, "误点了到货")


if __name__ == "__main__":
    unittest.main()
