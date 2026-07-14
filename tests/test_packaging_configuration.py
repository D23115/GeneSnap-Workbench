from __future__ import annotations

import tomllib
import unittest
from pathlib import Path


ROOT = Path(__file__).resolve().parents[1]


class PackagingConfigurationTests(unittest.TestCase):
    def test_build_script_uses_repository_venv_or_path_python(self):
        script = (ROOT / "packaging" / "build_windows.ps1").read_text(
            encoding="utf-8-sig"
        )

        self.assertNotIn("..\\..\\..\\.venv", script)
        self.assertIn('Join-Path $ProjectRoot ".venv\\Scripts\\python.exe"', script)
        self.assertIn("Get-Command python", script)

    def test_pyinstaller_is_declared_as_build_dependency(self):
        with (ROOT / "pyproject.toml").open("rb") as handle:
            pyproject = tomllib.load(handle)

        build_dependencies = pyproject["project"]["optional-dependencies"]["build"]
        self.assertTrue(
            any(item.lower().startswith("pyinstaller") for item in build_dependencies)
        )

    def test_build_script_waits_for_and_validates_packaged_smoke(self):
        script = (ROOT / "packaging" / "build_windows.ps1").read_text(
            encoding="utf-8-sig"
        )

        self.assertIn("Start-Process", script)
        self.assertIn("--smoke-test", script)
        self.assertIn("-Wait", script)
        self.assertIn("-PassThru", script)
        self.assertIn("ConvertFrom-Json", script)
        self.assertIn("$SmokePayload.ok", script)

    def test_portable_and_installer_outputs_include_license_notices(self):
        build_script = (ROOT / "packaging" / "build_windows.ps1").read_text(
            encoding="utf-8-sig"
        )
        installer_script = (ROOT / "packaging" / "GeneSnapWorkbench.iss").read_text(
            encoding="utf-8-sig"
        )

        for filename in ("LICENSE", "NOTICE", "THIRD_PARTY_NOTICES.md"):
            with self.subTest(filename=filename):
                self.assertIn(filename, build_script)
                self.assertIn(filename, installer_script)


if __name__ == "__main__":
    unittest.main()
