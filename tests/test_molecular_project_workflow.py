import tempfile
import unittest
from datetime import date, datetime, timezone
from pathlib import Path

from genesnap_workbench.app.application import (
    GeneSnapApplicationService,
    NewExpressionProjectCommand,
)
from tests.test_expression_vector_protocol import vector_and_protocol


NOW = datetime(2026, 7, 13, 1, 0, tzinfo=timezone.utc)


class MolecularProjectWorkflowTests(unittest.TestCase):
    def setUp(self):
        self.temp_dir = tempfile.TemporaryDirectory()
        self.addCleanup(self.temp_dir.cleanup)
        self.service = GeneSnapApplicationService(Path(self.temp_dir.name))
        vector, protocol = vector_and_protocol()
        self.stored = self.service.create_expression_project(
            NewExpressionProjectCommand(
                project_id="OE-STATE-001",
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

    def transition(self, action, **kwargs):
        self.stored = self.service.transition_molecular_project(
            self.stored.project_id,
            workflow_type="expression",
            action=action,
            actor="tester",
            occurred_at=NOW,
            **kwargs,
        )
        return self.stored

    def test_controlled_workflow_requires_prep_after_usable_clone(self):
        self.assertEqual(self.transition("mark_primers_ordered").snapshot.status, "primers_ordered")
        self.assertEqual(self.transition("mark_primers_arrived").snapshot.status, "primers_arrived")
        self.assertEqual(self.transition("start_cloning").snapshot.status, "cloning_in_progress")
        self.assertEqual(
            self.transition("mark_sent_for_sequencing").snapshot.status,
            "sequencing_in_progress",
        )

        construct = self.stored.design.constructs[0]
        exact = "T" * 80 + construct.insert_sequence + "C" * 80
        sequencing_dir = self.stored.project_folder / "03_sequencing"
        (sequencing_dir / "vendor_TP53-FL-1_result.seq").write_text(
            exact,
            encoding="ascii",
        )
        outcome = self.service.analyze_expression_sequencing(
            self.stored.project_id,
            actor="tester",
            analyzed_at=NOW,
        )
        self.stored = outcome.project

        with self.assertRaisesRegex(ValueError, "选择"):
            self.transition("start_plasmid_prep")

        prep = self.transition(
            "start_plasmid_prep",
            selected_clone_names=("TP53-FL-1",),
        )
        self.assertEqual(prep.snapshot.status, "plasmid_prep_in_progress")
        self.assertEqual(prep.snapshot.selected_prep_clone_names, ("TP53-FL-1",))
        self.assertEqual(
            self.transition("complete_plasmid_prep").snapshot.status,
            "plasmid_prep_completed",
        )
        completed = self.transition("complete_project")
        self.assertEqual(completed.snapshot.status, "project_completed")
        self.assertEqual(completed.snapshot.actual_completed_at, NOW)

    def test_invalid_status_jump_is_rejected_without_overwriting_history(self):
        original_history = self.stored.snapshot.status_history

        with self.assertRaises(ValueError):
            self.transition("start_cloning")

        restored = self.service.load_expression_project(self.stored.project_id)
        self.assertEqual(restored.snapshot.status, "design_completed")
        self.assertEqual(restored.snapshot.status_history, original_history)


if __name__ == "__main__":
    unittest.main()
