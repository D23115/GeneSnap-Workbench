import tempfile
import unittest
from pathlib import Path

from genesnap_workbench.app.application import GeneSnapApplicationService


class ApplicationSettingsTests(unittest.TestCase):
    def test_projects_root_can_be_changed_and_is_persisted(self):
        with tempfile.TemporaryDirectory() as data_dir, tempfile.TemporaryDirectory() as output_dir:
            service = GeneSnapApplicationService(Path(data_dir))
            chosen = Path(output_dir) / "GeneSnap 项目"

            selected = service.set_projects_root(chosen)

            reloaded = GeneSnapApplicationService(Path(data_dir))
            self.assertEqual(selected, chosen.resolve())
            self.assertEqual(reloaded.projects_root, selected)
            self.assertTrue(reloaded.has_custom_projects_root)


if __name__ == "__main__":
    unittest.main()
