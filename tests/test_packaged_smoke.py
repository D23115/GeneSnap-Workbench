import json
import tempfile
import unittest
from pathlib import Path

from PySide6.QtGui import QImage

from genesnap_workbench.app.main import (
    run_packaged_smoke,
    run_ui_screenshot,
)


class PackagedSmokeTests(unittest.TestCase):
    def test_smoke_creates_reopens_project_and_writes_report(self):
        with tempfile.TemporaryDirectory() as temp_dir:
            root = Path(temp_dir)
            report = root / "smoke-report.json"

            exit_code = run_packaged_smoke(root / "data", report)

            self.assertEqual(exit_code, 0)
            payload = json.loads(report.read_text(encoding="utf-8"))
            self.assertTrue(payload["ok"])
            self.assertEqual(
                set(payload["workflows"]),
                {"shrna", "expression", "reporter", "syn"},
            )
            for workflow in payload["workflows"].values():
                self.assertEqual(workflow["project_status"], "design_completed")
                self.assertGreaterEqual(workflow["artifact_count"], 5)
                self.assertTrue(all(Path(path).exists() for path in workflow["artifacts"]))

            screenshot = root / "ui.png"
            self.assertEqual(run_ui_screenshot(root / "data", screenshot), 0)
            self.assertTrue(screenshot.exists())
            self.assertGreater(screenshot.stat().st_size, 5000)
            image = QImage(str(screenshot))
            sampled = 0
            near_black = 0
            for y in range(0, image.height(), 8):
                for x in range(0, image.width(), 8):
                    color = image.pixelColor(x, y)
                    sampled += 1
                    if color.red() < 8 and color.green() < 8 and color.blue() < 8:
                        near_black += 1
            self.assertLess(near_black / sampled, 0.02)


if __name__ == "__main__":
    unittest.main()
