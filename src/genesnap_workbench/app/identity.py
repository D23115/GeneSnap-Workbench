"""Application name and icon shared by source and packaged entry points."""

from __future__ import annotations

from pathlib import Path

from PySide6.QtGui import QIcon
from PySide6.QtWidgets import QApplication


def application_icon_path() -> Path:
    return (
        Path(__file__).resolve().parents[1]
        / "resources"
        / "icons"
        / "genesnap_workbench.ico"
    )


def configure_application_identity(app: QApplication) -> None:
    app.setApplicationName("GeneSnap Workbench")
    app.setOrganizationName("GeneSnap")
    icon = QIcon(str(application_icon_path()))
    if not icon.isNull():
        app.setWindowIcon(icon)
