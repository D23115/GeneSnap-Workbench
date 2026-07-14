import tempfile
import unittest
from datetime import date, datetime, timezone
from pathlib import Path

from genesnap_workbench.app.application import (
    GeneSnapApplicationService,
    NewSYNProjectCommand,
)
from genesnap_workbench.domain.syn import (
    PlasmidPrepStatus,
    SYNAssemblyAttemptResult,
    SYNColonyPCRResult,
    SYNRoute,
    SYNSequencingResult,
)
from genesnap_workbench.project_workflow.syn_materials import MaterialReadiness
from genesnap_workbench.project_workflow.syn_service import SYNWorkflowService
from genesnap_workbench.project_workflow.syn_state import SYNStateTransitionService


NOW = datetime(2026, 7, 12, 16, 0, tzinfo=timezone.utc)
FIXTURES = Path(__file__).parent / "fixtures" / "public"


class SYNEndToEndTests(unittest.TestCase):
    def setUp(self):
        self.temp_dir = tempfile.TemporaryDirectory()
        self.addCleanup(self.temp_dir.cleanup)
        self.data_root = Path(self.temp_dir.name)
        self.app = GeneSnapApplicationService(self.data_root)
        self.state = SYNStateTransitionService()
        self.workflow = SYNWorkflowService()

    def run_case(self, filename: str, project_id: str):
        raw_sequence = (FIXTURES / filename).read_text(encoding="ascii")
        command = NewSYNProjectCommand(
            project_id=project_id,
            target_name=filename.removesuffix(".fasta"),
            raw_sequence=raw_sequence,
            input_format="fasta",
            linearization_site="EcoRV",
            received_date=date(2026, 7, 12),
            due_date=date(2026, 8, 3),
            actor="e2e",
            vector_sequence_confirmed=True,
        )
        prepared = self.app.prepare_syn_project(command, created_at=NOW)
        stored = self.app.save_prepared_syn_project(
            command,
            prepared,
            design_confirmation_reason="公开人工序列端到端复核",
            created_at=NOW,
        )
        original_revision = stored.snapshot.revision
        snapshot = self.state.mark_materials_ordered(
            stored.snapshot,
            event_id=f"{project_id}-ordered",
            actor="e2e",
            occurred_at=NOW,
        )
        snapshot = self.state.mark_materials_arrived(
            snapshot,
            resuspension_complete=True,
            event_id=f"{project_id}-arrived",
            actor="e2e",
            occurred_at=NOW,
        )
        snapshot = self.state.start_assembly(
            snapshot,
            MaterialReadiness(True, True, (), (), ()),
            confirm_missing=False,
            override_reason=None,
            event_id=f"{project_id}-start",
            override_id=f"{project_id}-unused",
            actor="e2e",
            occurred_at=NOW,
        )
        snapshot = self.workflow.begin_initial_attempt(
            snapshot,
            expected_revision=snapshot.revision,
            attempt_id=f"{project_id}-attempt-1",
            event_id=f"{project_id}-attempt-event",
            actor="e2e",
            occurred_at=NOW,
        )
        for step, next_step in (
            ("assembly_pcr", "amplification_pcr"),
            ("amplification_pcr", "vector_assembly_transformation"),
            ("vector_assembly_transformation", "colony_pcr"),
        ):
            self.assertEqual(snapshot.syn_assembly_substep, step)
            snapshot = self.workflow.record_assembly_step(
                snapshot,
                result=SYNAssemblyAttemptResult.SUCCESS,
                note="公开端到端案例",
                expected_revision=snapshot.revision,
                record_id=f"{project_id}-{step}-record",
                event_id=f"{project_id}-{step}-record-event",
                actor="e2e",
                occurred_at=NOW,
            )
            snapshot = self.state.advance_assembly_substep(
                snapshot,
                to_substep=next_step,
                event_id=f"{project_id}-{next_step}-advance-event",
                actor="e2e",
                occurred_at=NOW,
            )
        high_risk = any(
            risk.severity == "high_risk" for risk in stored.design.qc_result.risks
        )
        snapshot = self.workflow.start_initial_colony_screening(
            snapshot,
            target_name=stored.target_name,
            route=stored.design.module_plan.route,
            high_risk=high_risk,
            colony_count=None,
            expected_revision=snapshot.revision,
            attempt_id=f"{project_id}-attempt-unused",
            event_id=f"{project_id}-colonies",
            actor="e2e",
            occurred_at=NOW,
        )
        clone_id = snapshot.colonies[0].clone_id
        snapshot = self.workflow.record_colony_pcr(
            snapshot,
            clone_id=clone_id,
            result=SYNColonyPCRResult.POSITIVE,
            observed_note="预期条带",
            expected_revision=snapshot.revision,
            event_id=f"{project_id}-positive",
            actor="e2e",
            occurred_at=NOW,
        )
        snapshot = self.workflow.select_clones_for_prep(
            snapshot,
            clone_ids=(clone_id,),
            expected_revision=snapshot.revision,
            event_id=f"{project_id}-select-prep",
            actor="e2e",
            occurred_at=NOW,
        )
        snapshot = self.workflow.record_plasmid_prep(
            snapshot,
            clone_id=clone_id,
            status=PlasmidPrepStatus.COMPLETED,
            expected_revision=snapshot.revision,
            event_id=f"{project_id}-prep",
            actor="e2e",
            occurred_at=NOW,
        )
        snapshot = self.workflow.confirm_sequencing(
            snapshot,
            clone_id=clone_id,
            result=SYNSequencingResult.CORRECT,
            note="人工确认正确",
            expected_revision=snapshot.revision,
            confirmation_id=f"{project_id}-confirmation",
            event_id=f"{project_id}-confirmation-event",
            actor="e2e",
            occurred_at=NOW,
        )
        snapshot = self.workflow.complete_project(
            snapshot,
            user_confirmed=True,
            expected_revision=snapshot.revision,
            event_id=f"{project_id}-complete",
            actor="e2e",
            occurred_at=NOW,
        )
        self.app.repository.save_snapshot(
            project_id,
            snapshot,
            expected_revision=original_revision,
            updated_at=NOW,
        )
        return self.app.load_project(project_id)

    def test_public_cases_complete_and_survive_restart(self):
        low = self.run_case("syn_low_risk_600bp.fasta", "SYN-E2E-600")
        boundary = self.run_case("syn_boundary_1100bp.fasta", "SYN-E2E-1100")
        modular = self.run_case("syn_modular_1300bp.fasta", "SYN-E2E-1300")

        self.assertEqual(low.snapshot.status, "project_completed")
        self.assertEqual(low.design.module_plan.route, SYNRoute.SINGLE_POOL)
        self.assertEqual(len({item.clone_id for item in low.snapshot.colonies}), 8)
        self.assertTrue(boundary.design.module_plan.requires_confirmation)
        self.assertEqual(modular.design.module_plan.route, SYNRoute.MODULAR)
        self.assertEqual(
            len({item.clone_id for item in modular.snapshot.colonies}),
            12,
        )

        reopened = GeneSnapApplicationService(self.data_root)
        for project_id in ("SYN-E2E-600", "SYN-E2E-1100", "SYN-E2E-1300"):
            restored = reopened.load_project(project_id)
            self.assertEqual(restored.snapshot.status, "project_completed")
            artifacts = reopened.repository.list_artifacts(project_id)
            self.assertTrue(artifacts)
            self.assertTrue(all(item.path.exists() for item in artifacts))

    def test_fixtures_are_public_artificial_data_only(self):
        for path in FIXTURES.glob("syn_*bp.fasta"):
            content = path.read_text(encoding="ascii")
            self.assertIn("artificial_public_test_sequence", content)
            self.assertNotIn("private_reference", content)
            self.assertNotIn("E:\\", content)


if __name__ == "__main__":
    unittest.main()
