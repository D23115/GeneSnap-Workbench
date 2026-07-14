"""Create predictable Windows-safe project folder structures."""

from __future__ import annotations

import re
from dataclasses import dataclass
from pathlib import Path


STANDARD_SUBFOLDERS = (
    "01_design",
    "02_orders",
    "03_sequencing",
    "04_reports",
    "05_experiment_records",
    "99_archive",
)


def sanitize_windows_name(value: str) -> str:
    cleaned = re.sub(r'[<>:"/\\|?*\x00-\x1f]', "_", value.strip())
    cleaned = cleaned.rstrip(". ")
    if not cleaned:
        raise ValueError("文件夹名称字段不能为空")
    return cleaned


@dataclass(frozen=True, slots=True)
class ProjectWorkspace:
    root: Path
    subfolders: tuple[Path, ...]

    def folder(self, name: str) -> Path:
        for path in self.subfolders:
            if path.name == name:
                return path
        raise KeyError(name)


def create_project_folder(
    projects_root: Path,
    *,
    project_id: str,
    target_name: str,
    folder_suffix: str = "SYN",
) -> ProjectWorkspace:
    root = Path(projects_root)
    root.mkdir(parents=True, exist_ok=True)
    folder_name = "-".join(
        (
            sanitize_windows_name(project_id),
            sanitize_windows_name(target_name),
            sanitize_windows_name(folder_suffix),
        ),
    )
    project_root = root / folder_name
    project_root.mkdir()
    subfolders = tuple(project_root / name for name in STANDARD_SUBFOLDERS)
    for subfolder in subfolders:
        subfolder.mkdir()
    return ProjectWorkspace(root=project_root, subfolders=subfolders)
