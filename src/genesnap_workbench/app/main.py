"""Windows desktop entry point."""

from __future__ import annotations

import os
import argparse
import json
import random
import sys
import time
from datetime import date, datetime
from decimal import Decimal
from pathlib import Path

from PySide6.QtCore import QStandardPaths
from PySide6.QtWidgets import QApplication, QFileDialog

from genesnap_workbench.app.application import (
    GeneSnapApplicationService,
    NewExpressionProjectCommand,
    NewReporterProjectCommand,
    NewSYNProjectCommand,
    NewShRNAProjectCommand,
)
from genesnap_workbench.domain.shrna import BlastScreenStatus, ShRNACandidate
from genesnap_workbench.sequence_core.dna import reverse_complement
from genesnap_workbench.vector_library.models import (
    ExpressionVectorProtocol,
    ReporterVectorProtocol,
    VectorRecord,
)
from genesnap_workbench.app.desktop import MainWindow
from genesnap_workbench.app.identity import configure_application_identity


def default_data_root() -> Path:
    override = os.environ.get("GENESNAP_DATA_DIR")
    if override:
        return Path(override)
    app_data = QStandardPaths.writableLocation(
        QStandardPaths.StandardLocation.AppLocalDataLocation,
    )
    return Path(app_data)


def ensure_visible_projects_root(service: GeneSnapApplicationService) -> Path:
    """Ask once where project folders should live; never default them to hidden AppData."""
    if service.has_custom_projects_root:
        return service.projects_root
    documents = QStandardPaths.writableLocation(
        QStandardPaths.StandardLocation.DocumentsLocation,
    )
    suggested = Path(documents or Path.home() / "Documents") / "GeneSnap Workbench Projects"
    chosen = QFileDialog.getExistingDirectory(
        None,
        "选择 GeneSnap 项目保存位置",
        str(suggested.parent),
    )
    return service.set_projects_root(Path(chosen) if chosen else suggested)


def _smoke_expression_profile() -> tuple[VectorRecord, ExpressionVectorProtocol]:
    left_arm = "ACGTCAGTACGATCGTACGATCGTAGCA"
    right_top = "TGCATCGATGCTAGCTACGTA"
    sequence = "A" * 100 + left_arm + "AGT" + right_top + "C" * 100
    vector = VectorRecord.from_sequence(
        vector_record_id="packaged-smoke-expression-vector",
        structural_display_name="Packaged smoke expression vector",
        sequence=sequence,
    )
    left_boundary = 100 + len(left_arm)
    protocol = ExpressionVectorProtocol(
        protocol_id="packaged-smoke-expression",
        protocol_version_id="packaged-smoke-expression-v1",
        display_name="Packaged smoke expression",
        status="enabled",
        experimental_validation_status="unverified",
        vector_record_id=vector.vector_record_id,
        vector_checksum=vector.normalized_circular_sha256,
        workflow_type="expression",
        insertion_mode="confirmed_interval_with_homology_prefixes",
        left_boundary=left_boundary,
        right_boundary=left_boundary + 3,
        left_primer_homology=left_arm,
        right_primer_homology=reverse_complement(right_top),
        kozak_sequence="GCCACC",
        stop_codon_rule="remove_for_c_terminal_fusion",
        c_terminal_fusion_name="3xFLAG",
    )
    return vector, protocol


def _smoke_reporter_profile() -> tuple[VectorRecord, ReporterVectorProtocol]:
    left_arm = "TCGATCGTACGTAGCTAGCTACGTA"
    right_top = "GATCCGATCGTAGCTACGATC"
    sequence = "A" * 100 + left_arm + "AGCT" + right_top + "C" * 100
    vector = VectorRecord.from_sequence(
        vector_record_id="packaged-smoke-reporter-vector",
        structural_display_name="Packaged smoke reporter vector",
        sequence=sequence,
    )
    left_boundary = 100 + len(left_arm)
    protocol = ReporterVectorProtocol(
        protocol_id="packaged-smoke-reporter",
        protocol_version_id="packaged-smoke-reporter-v1",
        display_name="Packaged smoke reporter",
        status="enabled",
        experimental_validation_status="unverified",
        vector_record_id=vector.vector_record_id,
        vector_checksum=vector.normalized_circular_sha256,
        workflow_type="promoter_luciferase_reporter",
        insertion_mode="confirmed_interval_with_homology_prefixes",
        left_boundary=left_boundary,
        right_boundary=left_boundary + 4,
        left_primer_homology=left_arm,
        right_primer_homology=reverse_complement(right_top),
        default_sequencing_method="Nanopore",
    )
    return vector, protocol


def run_packaged_smoke(data_root: Path, report_path: Path) -> int:
    """Exercise all packaged workflows, exports, SQLite, and reopen behavior."""
    try:
        generator = random.Random(20260712)
        sequence = "".join(generator.choice("ACGT") for _ in range(600))
        service = GeneSnapApplicationService(data_root)
        syn_command = NewSYNProjectCommand(
            project_id="PACKAGED-SMOKE-SYN",
            target_name="artificial-packaged-smoke",
            raw_sequence=sequence,
            input_format="plain",
            linearization_site="EcoRV",
            received_date=date(2026, 7, 12),
            due_date=date(2026, 8, 3),
            actor="packaged-smoke",
            vector_sequence_confirmed=True,
        )
        created_at = datetime.now().astimezone()
        prepared = service.prepare_syn_project(syn_command, created_at=created_at)
        service.save_prepared_syn_project(
            syn_command,
            prepared,
            design_confirmation_reason="成品自动验收人工序列",
            created_at=created_at,
        )
        service.create_shrna_project(
            NewShRNAProjectCommand(
                project_id="PACKAGED-SMOKE-KD",
                gene_symbol="TP53",
                species="human",
                cds_sequence="ATG" * 300,
                candidates=(
                    ShRNACandidate(
                        candidate_id="packaged-smoke-target",
                        target_sequence="GACTCCAGTGGTAATCTACTG",
                        start_position=120,
                        intrinsic_score=Decimal("9.1"),
                        source_rank=1,
                        blast_status=BlastScreenStatus.MANUALLY_ACCEPTED,
                    ),
                ),
                target_count=1,
                clones_per_target=5,
                received_date=date(2026, 7, 12),
                due_date=date(2026, 7, 23),
                actor="packaged-smoke",
                vector_sequence_confirmed=True,
            ),
            created_at=created_at,
        )
        expression_vector, expression_protocol = _smoke_expression_profile()
        service.create_expression_project(
            NewExpressionProjectCommand(
                project_id="PACKAGED-SMOKE-OE",
                gene_symbol="TP53",
                species="human",
                source_cds="ATG" + "GCT" * 120 + "TAA",
                construct_lines=("FL", "1-80aa"),
                received_date=date(2026, 7, 12),
                due_date=date(2026, 7, 23),
                actor="packaged-smoke",
                vector=expression_vector,
                protocol=expression_protocol,
            ),
            created_at=created_at,
        )
        reporter_vector, reporter_protocol = _smoke_reporter_profile()
        service.create_reporter_project(
            NewReporterProjectCommand(
                project_id="PACKAGED-SMOKE-RPT",
                gene_symbol="SGK1",
                species="human",
                promoter_sequence="ACGT" * 500,
                construct_lines=("WT", "P1000", "P500"),
                mutation_definitions=(),
                received_date=date(2026, 7, 12),
                due_date=date(2026, 7, 23),
                actor="packaged-smoke",
                vector=reporter_vector,
                protocol=reporter_protocol,
            ),
            created_at=created_at,
        )
        reopened = GeneSnapApplicationService(data_root)
        workflow_records = {
            "syn": (
                reopened.load_project("PACKAGED-SMOKE-SYN"),
                reopened.repository.list_artifacts("PACKAGED-SMOKE-SYN"),
            ),
            "shrna": (
                reopened.load_shrna_project("PACKAGED-SMOKE-KD"),
                reopened.shrna_repository.list_artifacts("PACKAGED-SMOKE-KD"),
            ),
            "expression": (
                reopened.load_expression_project("PACKAGED-SMOKE-OE"),
                reopened.expression_repository.list_artifacts("PACKAGED-SMOKE-OE"),
            ),
            "reporter": (
                reopened.load_reporter_project("PACKAGED-SMOKE-RPT"),
                reopened.reporter_repository.list_artifacts("PACKAGED-SMOKE-RPT"),
            ),
        }
        payload = {
            "ok": True,
            "workflows": {
                name: {
                    "project_status": stored.snapshot.status,
                    "artifact_count": len(artifacts),
                    "artifacts": [str(item.path) for item in artifacts],
                    "design_version_id": stored.design.design_version_id,
                }
                for name, (stored, artifacts) in workflow_records.items()
            },
        }
        exit_code = 0
    except Exception as error:
        payload = {
            "ok": False,
            "error_type": type(error).__name__,
            "error": str(error),
        }
        exit_code = 1
    report_path.parent.mkdir(parents=True, exist_ok=True)
    report_path.write_text(
        json.dumps(payload, ensure_ascii=False, indent=2) + "\n",
        encoding="utf-8",
    )
    return exit_code


def _parse_args(argv: list[str]):
    parser = argparse.ArgumentParser(add_help=False)
    parser.add_argument("--smoke-test", action="store_true")
    parser.add_argument("--data-dir", type=Path)
    parser.add_argument("--smoke-report", type=Path)
    parser.add_argument("--ui-screenshot", type=Path)
    return parser.parse_known_args(argv)[0]


def run_ui_screenshot(data_root: Path, screenshot_path: Path) -> int:
    """Render the packaged dashboard offscreen for artifact-level UI checks."""
    os.environ["QT_QPA_PLATFORM"] = "offscreen"
    app = QApplication.instance() or QApplication([])
    configure_application_identity(app)
    service = GeneSnapApplicationService(data_root)
    window = MainWindow(service, actor="packaged-ui-smoke")
    window.resize(1440, 860)
    window.show()
    for _ in range(20):
        app.processEvents()
        time.sleep(0.05)
    window.repaint()
    app.processEvents()
    screenshot_path.parent.mkdir(parents=True, exist_ok=True)
    saved = window.grab().save(str(screenshot_path))
    window.close()
    app.processEvents()
    return 0 if saved else 1


def main() -> int:
    args = _parse_args(sys.argv[1:])
    if args.smoke_test:
        if args.data_dir is None or args.smoke_report is None:
            return 2
        return run_packaged_smoke(args.data_dir, args.smoke_report)
    if args.ui_screenshot is not None:
        if args.data_dir is None:
            return 2
        return run_ui_screenshot(args.data_dir, args.ui_screenshot)
    app = QApplication.instance() or QApplication(sys.argv)
    configure_application_identity(app)
    service = GeneSnapApplicationService(default_data_root())
    ensure_visible_projects_root(service)
    window = MainWindow(service, actor=os.environ.get("USERNAME", "本机用户"))
    window.show()
    return app.exec()


if __name__ == "__main__":
    raise SystemExit(main())
