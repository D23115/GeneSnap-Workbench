"""Independent manual visibility history for all project workflow types."""

from __future__ import annotations

from dataclasses import asdict, dataclass
from datetime import datetime
import json
from pathlib import Path


@dataclass(frozen=True, slots=True)
class ProjectVisibilityEvent:
    project_id: str
    hidden: bool
    reason: str
    actor: str
    occurred_at: datetime


class LocalProjectVisibilityStore:
    def __init__(self, path: Path) -> None:
        self.path = Path(path)

    def _load(self) -> list[dict]:
        if not self.path.exists():
            return []
        return json.loads(self.path.read_text(encoding="utf-8"))

    def append(self, event: ProjectVisibilityEvent) -> None:
        if not event.reason.strip():
            raise ValueError("隐藏或恢复显示时必须填写原因")
        rows = self._load()
        payload = asdict(event)
        payload["occurred_at"] = event.occurred_at.isoformat()
        rows.append(payload)
        self.path.parent.mkdir(parents=True, exist_ok=True)
        temporary = self.path.with_suffix(self.path.suffix + ".tmp")
        temporary.write_text(
            json.dumps(rows, ensure_ascii=False, indent=2, sort_keys=True) + "\n",
            encoding="utf-8",
        )
        temporary.replace(self.path)

    def is_hidden(self, project_id: str) -> bool:
        hidden = False
        for row in self._load():
            if row.get("project_id") == project_id:
                hidden = bool(row.get("hidden"))
        return hidden

    def history(self, project_id: str) -> tuple[ProjectVisibilityEvent, ...]:
        return tuple(
            ProjectVisibilityEvent(
                project_id=row["project_id"],
                hidden=bool(row["hidden"]),
                reason=row["reason"],
                actor=row["actor"],
                occurred_at=datetime.fromisoformat(row["occurred_at"]),
            )
            for row in self._load()
            if row.get("project_id") == project_id
        )
