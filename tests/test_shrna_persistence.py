import tempfile
import unittest
from datetime import date, datetime, timezone
from decimal import Decimal
from pathlib import Path

from genesnap_workbench.domain.shrna import (
    BlastScreenStatus,
    ShRNAAuditEvent,
    ShRNACandidate,
    ShRNADesignInput,
    ShRNAProjectSnapshot,
)
from genesnap_workbench.sequence_core.shrna import create_shrna_design
from genesnap_workbench.storage.shrna_repository import SQLiteShRNAProjectRepository
from genesnap_workbench.vector_library.starters import load_public_plko1_puro_starter


NOW = datetime(2026, 7, 12, 18, 0, tzinfo=timezone.utc)


def make_design():
    vector, protocol = load_public_plko1_puro_starter(user_confirmed=True)
    design_input = ShRNADesignInput(
        project_id="KD-PERSIST-001",
        gene_symbol="TP53",
        species="human",
        cds_sequence="ATG" * 300,
        vector_protocol_version_id=protocol.protocol_version_id,
        target_count=1,
    )
    candidate = ShRNACandidate(
        candidate_id="candidate-1",
        target_sequence="GACTCCAGTGGTAATCTACTG",
        start_position=120,
        intrinsic_score=Decimal("9.1"),
        source_rank=1,
        blast_status=BlastScreenStatus.PASS,
    )
    return create_shrna_design(
        design_input,
        (candidate,),
        vector,
        protocol,
        design_version_id="KD-PERSIST-001-v1",
        created_at=NOW,
    )


class ShRNAPersistenceTests(unittest.TestCase):
    def setUp(self):
        self.temp_dir = tempfile.TemporaryDirectory()
        self.addCleanup(self.temp_dir.cleanup)
        self.database = Path(self.temp_dir.name) / "workbench.db"
        self.repository = SQLiteShRNAProjectRepository(self.database)
        self.repository.initialize()

    def test_design_and_snapshot_round_trip_after_reopen(self):
        design = make_design()
        snapshot = ShRNAProjectSnapshot(
            project_id=design.project_id,
            revision=1,
            status="design_completed",
            active_design_version_id=design.design_version_id,
            clone_results=(),
            status_history=(
                ShRNAAuditEvent(
                    event_id="event-1",
                    event_type="complete_design",
                    occurred_at=NOW,
                    actor="tester",
                    from_status="recorded",
                    to_status="design_completed",
                ),
            ),
        )
        project_folder = Path(self.temp_dir.name) / "projects" / "KD-PERSIST-001-TP53-KD"
        project_folder.mkdir(parents=True)

        self.repository.create_project(
            project_id=design.project_id,
            gene_symbol=design.gene_symbol,
            species=design.species,
            received_date=date(2026, 7, 12),
            due_date=date(2026, 7, 23),
            project_folder=project_folder,
            design=design,
            snapshot=snapshot,
            created_at=NOW,
        )

        reopened = SQLiteShRNAProjectRepository(self.database)
        reopened.initialize()
        stored = reopened.load_project(design.project_id)

        self.assertEqual(stored.design, design)
        self.assertEqual(stored.snapshot, snapshot)
        self.assertEqual(stored.folder_suffix, "KD")
        self.assertEqual(reopened.list_projects()[0].workflow_type, "shrna_knockdown")


if __name__ == "__main__":
    unittest.main()
