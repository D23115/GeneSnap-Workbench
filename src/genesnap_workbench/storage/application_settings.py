"""Small local settings store for paths and one-time integration acknowledgements."""

from __future__ import annotations

import json
from pathlib import Path


class LocalApplicationSettingsStore:
    def __init__(self, path: Path) -> None:
        self.path = Path(path)

    def load(self) -> dict[str, object]:
        if not self.path.exists():
            return {}
        try:
            value = json.loads(self.path.read_text(encoding="utf-8"))
        except (OSError, json.JSONDecodeError) as error:
            raise ValueError(f"应用设置无法读取：{self.path}") from error
        if not isinstance(value, dict):
            raise ValueError("应用设置格式无效")
        return value

    def update(self, **values: object) -> None:
        settings = self.load()
        settings.update(values)
        self.path.parent.mkdir(parents=True, exist_ok=True)
        temporary = self.path.with_suffix(self.path.suffix + ".tmp")
        temporary.write_text(
            json.dumps(settings, ensure_ascii=False, indent=2),
            encoding="utf-8",
        )
        temporary.replace(self.path)
