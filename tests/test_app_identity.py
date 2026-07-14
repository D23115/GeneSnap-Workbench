import os
import unittest
from pathlib import Path

os.environ.setdefault("QT_QPA_PLATFORM", "offscreen")

from PySide6.QtGui import QIcon
from PySide6.QtWidgets import QApplication, QWidget

from genesnap_workbench.app.identity import (
    application_icon_path,
    configure_application_identity,
)


class ApplicationIdentityTests(unittest.TestCase):
    @classmethod
    def setUpClass(cls):
        cls.qt_app = QApplication.instance() or QApplication([])

    def test_application_icon_has_small_and_large_windows_sizes(self):
        icon_path = application_icon_path()

        self.assertTrue(icon_path.is_file())
        self.assertEqual(icon_path.suffix.lower(), ".ico")
        icon = QIcon(str(icon_path))
        self.assertFalse(icon.isNull())
        sizes = {(size.width(), size.height()) for size in icon.availableSizes()}
        self.assertIn((16, 16), sizes)
        self.assertIn((256, 256), sizes)

    def test_application_identity_sets_name_and_inherited_window_icon(self):
        configure_application_identity(self.qt_app)
        window = QWidget()
        self.addCleanup(window.close)

        self.assertEqual(self.qt_app.applicationName(), "GeneSnap Workbench")
        self.assertEqual(self.qt_app.organizationName(), "GeneSnap")
        self.assertFalse(self.qt_app.windowIcon().isNull())
        self.assertFalse(window.windowIcon().isNull())

    def test_windows_packaging_uses_the_same_icon(self):
        project_root = Path(__file__).resolve().parents[1]
        spec_text = (project_root / "packaging" / "genesnap_workbench.spec").read_text(
            encoding="utf-8",
        )
        installer_text = (
            project_root / "packaging" / "GeneSnapWorkbench.iss"
        ).read_text(encoding="utf-8")

        self.assertIn('icon=str(resources_root / "icons" / "genesnap_workbench.ico")', spec_text)
        self.assertIn("SetupIconFile=", installer_text)
        self.assertIn("UninstallDisplayIcon=", installer_text)


if __name__ == "__main__":
    unittest.main()
