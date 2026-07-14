import tempfile
import unittest
from datetime import date, datetime, timezone
from pathlib import Path

from genesnap_workbench.domain.expression import (
    ExpressionAuditEvent,
    ExpressionDesignInput,
    ExpressionProjectSnapshot,
)
from genesnap_workbench.sequence_core.expression import create_expression_design
from genesnap_workbench.storage.expression_repository import (
    SQLiteExpressionProjectRepository,
)
from genesnap_workbench.vector_library.expression import (
    apply_expression_protocol,
    expression_rules_from_protocol,
)
from tests.test_expression_vector_protocol import vector_and_protocol


NOW = datetime(2026, 7, 12, 23, 0, tzinfo=timezone.utc)


class ExpressionPersistenceTests(unittest.TestCase):
    def setUp(self):
        self.temp_dir = tempfile.TemporaryDirectory()
        self.addCleanup(self.temp_dir.cleanup)
        self.database = Path(self.temp_dir.name) / "workbench.db"
        self.repository = SQLiteExpressionProjectRepository(self.database)
        self.repository.initialize()

    def test_multi_construct_design_and_vector_result_round_trip(self):
        vector, protocol = vector_and_protocol()
        design = create_expression_design(
            ExpressionDesignInput(
                project_id="OE-PERSIST-001",
                gene_symbol="TP53",
                species="human",
                source_cds="ATG" + "GCT" * 120 + "TAA",
                construct_lines=("FL", "1-80aa", "Δ81-100"),
            ),
            expression_rules_from_protocol(protocol),
            design_version_id="OE-PERSIST-001-v1",
            created_at=NOW,
        )
        vector_design = apply_expression_protocol(design, vector, protocol)
        snapshot = ExpressionProjectSnapshot(
            project_id=design.project_id,
            revision=1,
            status="design_completed",
            active_design_version_id=design.design_version_id,
            status_history=(
                ExpressionAuditEvent(
                    event_id="event-1",
                    event_type="complete_design",
                    occurred_at=NOW,
                    actor="tester",
                    from_status="recorded",
                    to_status="design_completed",
                ),
            ),
        )
        project_folder = (
            Path(self.temp_dir.name) / "projects" / "OE-PERSIST-001-TP53-OE"
        )
        project_folder.mkdir(parents=True)

        self.repository.create_project(
            project_id=design.project_id,
            gene_symbol=design.gene_symbol,
            species=design.species,
            vector_name="人工表达载体",
            received_date=date(2026, 7, 12),
            due_date=date(2026, 7, 23),
            project_folder=project_folder,
            design=design,
            vector_design=vector_design,
            snapshot=snapshot,
            created_at=NOW,
        )

        reopened = SQLiteExpressionProjectRepository(self.database)
        reopened.initialize()
        stored = reopened.load_project(design.project_id)

        self.assertEqual(stored.design, design)
        self.assertEqual(stored.vector_design, vector_design)
        self.assertEqual(stored.snapshot, snapshot)
        self.assertEqual(stored.folder_suffix, "OE")
        self.assertEqual(reopened.list_projects()[0].workflow_type, "expression")
        self.assertEqual(reopened.list_projects()[0].construct_count, 3)


if __name__ == "__main__":
    unittest.main()
