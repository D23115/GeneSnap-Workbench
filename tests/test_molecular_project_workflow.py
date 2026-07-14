import json
import tempfile
import unittest
from dataclasses import fields
from datetime import date, datetime, timezone
from decimal import Decimal
from pathlib import Path

from genesnap_workbench.app.application import (
    GeneSnapApplicationService,
    NewExpressionProjectCommand,
    NewReporterProjectCommand,
    NewShRNAProjectCommand,
)
from genesnap_workbench.domain.expression import ExpressionProjectSnapshot
from genesnap_workbench.domain.reporter import ReporterProjectSnapshot
from genesnap_workbench.domain.shrna import (
    BlastScreenStatus,
    ShRNACandidate,
    ShRNAProjectSnapshot,
)
from genesnap_workbench.storage.codec import loads_record
from tests.test_expression_vector_protocol import vector_and_protocol
from tests.test_reporter_exports import export_promoter_sequence
from tests.test_reporter_vector_protocol import (
    vector_and_protocol as reporter_vector_and_protocol,
)


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

    def create_shrna_project(self):
        return self.service.create_shrna_project(
            NewShRNAProjectCommand(
                project_id="KD-TRACK-001",
                gene_symbol="TP53",
                species="human",
                cds_sequence="ATG" * 300,
                candidates=(
                    ShRNACandidate(
                        candidate_id="manual-1",
                        target_sequence="GACTCCAGTGGTAATCTACTG",
                        start_position=None,
                        intrinsic_score=Decimal("0"),
                        source_rank=1,
                        blast_status=BlastScreenStatus.MANUALLY_ACCEPTED,
                    ),
                ),
                target_count=1,
                clones_per_target=5,
                received_date=date(2026, 7, 13),
                due_date=date(2026, 7, 24),
                actor="tester",
                vector_sequence_confirmed=True,
            ),
            created_at=NOW,
        )

    def create_reporter_project(self):
        vector, protocol = reporter_vector_and_protocol()
        return self.service.create_reporter_project(
            NewReporterProjectCommand(
                project_id="RPT-TRACK-001",
                gene_symbol="SGK1",
                species="human",
                promoter_sequence=export_promoter_sequence(),
                construct_lines=("WT",),
                mutation_definitions=(),
                received_date=date(2026, 7, 13),
                due_date=date(2026, 7, 24),
                actor="tester",
                vector=vector,
                protocol=protocol,
            ),
            created_at=NOW,
        )

    def test_old_molecular_snapshots_load_empty_tracking_number_defaults(self):
        tracking_fields = (
            "internal_project_no",
            "primer_submission_no",
            "primer_vendor_order_no",
        )
        snapshots = (
            (
                ShRNAProjectSnapshot,
                {
                    "project_id": "KD-OLD-001",
                    "revision": 1,
                    "status": "design_completed",
                    "active_design_version_id": "KD-OLD-001-v1",
                    "clone_results": {"__tuple__": []},
                    "status_history": {"__tuple__": []},
                },
            ),
            (
                ExpressionProjectSnapshot,
                {
                    "project_id": "OE-OLD-001",
                    "revision": 1,
                    "status": "design_completed",
                    "active_design_version_id": "OE-OLD-001-v1",
                    "status_history": {"__tuple__": []},
                },
            ),
            (
                ReporterProjectSnapshot,
                {
                    "project_id": "RPT-OLD-001",
                    "revision": 1,
                    "status": "design_completed",
                    "active_design_version_id": "RPT-OLD-001-v1",
                    "status_history": {"__tuple__": []},
                },
            ),
        )

        for snapshot_type, old_fields in snapshots:
            with self.subTest(snapshot_type=snapshot_type.__name__):
                payload = json.dumps(
                    {"__type__": snapshot_type.__name__, "fields": old_fields},
                )
                snapshot = loads_record(payload)

                self.assertEqual(
                    tuple(getattr(snapshot, name) for name in tracking_fields),
                    ("", "", ""),
                )
                self.assertEqual(
                    tuple(item.name for item in fields(snapshot_type)[-3:]),
                    tracking_fields,
                )

    def test_mark_primers_ordered_saves_three_optional_numbers_atomically(self):
        original = self.stored.snapshot

        ordered = self.transition(
            "mark_primers_ordered",
            internal_project_no=" 内 部/001 ",
            primer_submission_no=" PR-SUB-001 ",
            primer_vendor_order_no=" PR-ORDER-001 ",
        )

        self.assertEqual(ordered.snapshot.status, "primers_ordered")
        self.assertEqual(ordered.snapshot.internal_project_no, "内 部/001")
        self.assertEqual(ordered.snapshot.primer_submission_no, "PR-SUB-001")
        self.assertEqual(ordered.snapshot.primer_vendor_order_no, "PR-ORDER-001")
        self.assertEqual(ordered.snapshot.revision, original.revision + 1)
        self.assertEqual(ordered.snapshot.status_history[:-1], original.status_history)
        self.assertEqual(
            ordered.snapshot.status_history[-1].event_type,
            "mark_primers_ordered",
        )

    def test_tracking_number_update_requires_note_and_appends_history(self):
        original = self.stored.snapshot

        with self.assertRaisesRegex(ValueError, "修改说明"):
            self.service.update_molecular_tracking_numbers(
                self.stored.project_id,
                workflow_type="expression",
                internal_project_no="INT-002",
                primer_submission_no="PR-SUB-002",
                primer_vendor_order_no="PR-ORDER-002",
                actor="tester",
                occurred_at=NOW,
                note="  ",
            )

        unchanged = self.service.load_expression_project(self.stored.project_id)
        self.assertEqual(unchanged.snapshot, original)

        updated = self.service.update_molecular_tracking_numbers(
            self.stored.project_id,
            workflow_type="expression",
            internal_project_no=" INT-002 ",
            primer_submission_no=" PR-SUB-002 ",
            primer_vendor_order_no=" PR-ORDER-002 ",
            actor="tester",
            occurred_at=NOW,
            note="收到供应商回执后补录",
        )

        self.assertEqual(updated.snapshot.internal_project_no, "INT-002")
        self.assertEqual(updated.snapshot.primer_submission_no, "PR-SUB-002")
        self.assertEqual(updated.snapshot.primer_vendor_order_no, "PR-ORDER-002")
        self.assertEqual(updated.snapshot.status_history[:-1], original.status_history)
        self.assertEqual(
            updated.snapshot.status_history[-1].event_type,
            "update_molecular_tracking_numbers",
        )
        audit_note = updated.snapshot.status_history[-1].note
        self.assertIn("收到供应商回执后补录", audit_note)
        self.assertIn("内部编号：- -> INT-002", audit_note)
        self.assertIn("引物送单号：- -> PR-SUB-002", audit_note)
        self.assertIn("引物订单号：- -> PR-ORDER-002", audit_note)

    def test_all_three_molecular_summaries_expose_tracking_numbers(self):
        shrna = self.create_shrna_project()
        reporter = self.create_reporter_project()
        projects = (
            (self.stored.project_id, "expression", "INT-OE", "SUB-OE", "ORDER-OE"),
            (shrna.project_id, "shrna_knockdown", "INT-KD", "SUB-KD", "ORDER-KD"),
            (
                reporter.project_id,
                "promoter_luciferase_reporter",
                "INT-RPT",
                "SUB-RPT",
                "ORDER-RPT",
            ),
        )
        for project_id, workflow_type, internal_no, submission_no, order_no in projects:
            self.service.transition_molecular_project(
                project_id,
                workflow_type=workflow_type,
                action="mark_primers_ordered",
                actor="tester",
                occurred_at=NOW,
                internal_project_no=internal_no,
                primer_submission_no=submission_no,
                primer_vendor_order_no=order_no,
            )

        summaries = {item.project_id: item for item in self.service.list_all_projects()}
        for project_id, _, internal_no, submission_no, order_no in projects:
            with self.subTest(project_id=project_id):
                summary = summaries[project_id]
                self.assertEqual(summary.internal_project_no, internal_no)
                self.assertEqual(summary.primer_submission_no, submission_no)
                self.assertEqual(summary.primer_vendor_order_no, order_no)

    def test_controlled_workflow_requires_prep_after_usable_clone(self):
        ordered = self.transition("mark_primers_ordered")
        self.assertEqual(ordered.snapshot.status, "primers_ordered")
        self.assertEqual(ordered.snapshot.internal_project_no, "")
        self.assertEqual(ordered.snapshot.primer_submission_no, "")
        self.assertEqual(ordered.snapshot.primer_vendor_order_no, "")
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
