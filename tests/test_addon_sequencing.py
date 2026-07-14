import tempfile
import unittest
from datetime import date, datetime, timezone
from decimal import Decimal
from pathlib import Path

from genesnap_workbench.app.application import GeneSnapApplicationService, NewShRNAProjectCommand
from genesnap_workbench.domain.shrna import BlastScreenStatus, ShRNACandidate


NOW = datetime(2026, 7, 13, 5, 0, tzinfo=timezone.utc)


class AddOnSequencingTests(unittest.TestCase):
    def setUp(self):
        self.temp_dir = tempfile.TemporaryDirectory()
        self.addCleanup(self.temp_dir.cleanup)
        self.service = GeneSnapApplicationService(Path(self.temp_dir.name))
        candidates = tuple(
            ShRNACandidate(
                candidate_id=f"candidate-{index}",
                target_sequence=sequence,
                start_position=index * 150,
                intrinsic_score=Decimal(str(10 - index)),
                source_rank=index,
                blast_status=BlastScreenStatus.PASS,
            )
            for index, sequence in enumerate(
                (
                    "GACTCCAGTGGTAATCTACTG",
                    "GCAGTCACAGCACATGACGGA",
                    "GCTGCTCAGATAGCGATGGTC",
                ),
                start=1,
            )
        )
        self.stored = self.service.create_shrna_project(
            NewShRNAProjectCommand(
                project_id="KD-ADD-001",
                gene_symbol="TP53",
                species="human",
                cds_sequence="ATG" * 400,
                candidates=candidates,
                target_count=3,
                clones_per_target=5,
                received_date=date(2026, 7, 13),
                due_date=date(2026, 7, 24),
                actor="tester",
                vector_sequence_confirmed=True,
            ),
            created_at=NOW,
        )
        for action in (
            "mark_primers_ordered",
            "mark_primers_arrived",
            "start_cloning",
            "mark_sent_for_sequencing",
        ):
            self.stored = self.service.transition_molecular_project(
                self.stored.project_id,
                workflow_type="shrna_knockdown",
                action=action,
                actor="tester",
                occurred_at=NOW,
            )

    def _write_complete_round_with_only_target_two_failed(self):
        sequencing_dir = self.stored.project_folder / "03_sequencing"
        for target in self.stored.design.targets:
            for clone_name in target.clone_names:
                sequence = (
                    "A" * 80 + target.oligos.forward_sequence + "C" * 80
                    if target.target_no in {1, 3}
                    else "A" * 300
                )
                (sequencing_dir / f"{clone_name}_result.seq").write_text(
                    sequence,
                    encoding="ascii",
                )

    def test_preview_only_adds_failed_target_and_continues_clone_numbers(self):
        self._write_complete_round_with_only_target_two_failed()
        analyzed = self.service.analyze_shrna_sequencing(
            self.stored.project_id,
            actor="tester",
            analyzed_at=NOW,
        ).project

        preview = self.service.preview_addon_sequencing(
            analyzed.project_id,
            workflow_type="shrna_knockdown",
        )

        self.assertEqual(preview.affected_owner_labels, ("shRNA-2",))
        self.assertEqual(preview.clones_per_owner, 10)
        self.assertEqual(
            preview.sample_names,
            tuple(f"TP53-2-{index}" for index in range(6, 16)),
        )
        self.assertEqual(preview.round_no, 2)

        before = self.service.load_shrna_project(analyzed.project_id)
        confirmed = self.service.confirm_addon_sequencing(
            preview,
            actor="tester",
            occurred_at=NOW,
            sequencing_vendor_name="擎科",
        )

        self.assertEqual(len(before.snapshot.sequencing_submissions), 1)
        self.assertEqual(confirmed.snapshot.status, "add_on_in_progress")
        self.assertEqual(len(confirmed.snapshot.sequencing_submissions), 2)
        submission = confirmed.snapshot.sequencing_submissions[-1]
        self.assertEqual(submission.submission_kind, "add_on")
        self.assertEqual(submission.sample_names, preview.sample_names)
        self.assertTrue(Path(submission.form_path).exists())
        self.assertIn("TP53_加测2", Path(submission.form_path).name)

        target_two = confirmed.design.targets[1]
        sequencing_dir = confirmed.project_folder / "03_sequencing"
        for sample_name in preview.sample_names:
            (sequencing_dir / f"{sample_name}_result.seq").write_text(
                "A" * 80 + target_two.oligos.forward_sequence + "C" * 80,
                encoding="ascii",
            )
        reanalyzed = self.service.analyze_shrna_sequencing(
            confirmed.project_id,
            actor="tester",
            analyzed_at=NOW,
        ).project
        latest = {item.clone_name: item for item in reanalyzed.snapshot.clone_results}
        self.assertEqual(latest["TP53-2-6"].status, "pass")
        self.assertEqual(reanalyzed.snapshot.sequencing_submissions[-1].status, "analyzed")

    def test_rework_uses_new_attempt_and_n_suffix_instead_of_addon_numbering(self):
        self._write_complete_round_with_only_target_two_failed()
        analyzed = self.service.analyze_shrna_sequencing(
            self.stored.project_id,
            actor="tester",
            analyzed_at=NOW,
        ).project

        reworking = self.service.start_molecular_rework(
            analyzed.project_id,
            workflow_type="shrna_knockdown",
            actor="tester",
            occurred_at=NOW,
            note="第二个 target 重新连接转化",
        )
        preview = self.service.preview_rework_submission(
            reworking.project_id,
            workflow_type="shrna_knockdown",
        )

        self.assertEqual(reworking.snapshot.status, "rework_in_progress")
        self.assertEqual(reworking.snapshot.experiment_attempt_no, 2)
        self.assertEqual(preview.affected_owner_labels, ("shRNA-2",))
        self.assertEqual(
            preview.sample_names,
            tuple(f"TP53_2n_{index}" for index in range(1, 6)),
        )

        submitted = self.service.confirm_rework_submission(
            preview,
            actor="tester",
            occurred_at=NOW,
            sequencing_vendor_name="擎科",
        )
        submission = submitted.snapshot.sequencing_submissions[-1]
        self.assertEqual(submitted.snapshot.status, "sequencing_in_progress")
        self.assertEqual(submission.submission_kind, "post_rework")
        self.assertEqual(submission.experiment_attempt_no, 2)
        self.assertIn("TP53_加测2", Path(submission.form_path).name)

        target_two = submitted.design.targets[1]
        (submitted.project_folder / "03_sequencing" / "TP53_2n_1_result.seq").write_text(
            "A" * 80 + target_two.oligos.forward_sequence + "C" * 80,
            encoding="ascii",
        )
        analyzed_rework = self.service.analyze_shrna_sequencing(
            submitted.project_id,
            actor="tester",
            analyzed_at=NOW,
        ).project
        latest = {item.clone_name: item for item in analyzed_rework.snapshot.clone_results}
        self.assertEqual(latest["TP53_2n_1"].status, "pass")


if __name__ == "__main__":
    unittest.main()
