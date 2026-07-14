"""shRNA 项目和追加式输出历史的 SQLite 持久化。"""

from __future__ import annotations

import sqlite3
from contextlib import contextmanager
from dataclasses import dataclass
from datetime import date, datetime
from pathlib import Path

from genesnap_workbench.domain.shrna import ShRNADesignVersion, ShRNAProjectSnapshot
from genesnap_workbench.template_engine.syn_exports import GeneratedArtifact

from .codec import dumps_record, loads_record
from .syn_repository import DuplicateProjectError, StorageRevisionConflict


@dataclass(frozen=True, slots=True)
class ShRNAProjectSummary:
    project_id: str
    gene_symbol: str
    species: str
    status: str
    received_date: date
    due_date: date
    project_folder: Path
    folder_suffix: str
    revision: int
    workflow_type: str = "shrna_knockdown"


@dataclass(frozen=True, slots=True)
class StoredShRNAProject:
    project_id: str
    gene_symbol: str
    species: str
    received_date: date
    due_date: date
    project_folder: Path
    folder_suffix: str
    design: ShRNADesignVersion
    snapshot: ShRNAProjectSnapshot
    created_at: datetime
    updated_at: datetime
    workflow_type: str = "shrna_knockdown"


class SQLiteShRNAProjectRepository:
    def __init__(self, database_path: Path) -> None:
        self.database_path = Path(database_path)

    def _connect(self) -> sqlite3.Connection:
        connection = sqlite3.connect(self.database_path)
        connection.row_factory = sqlite3.Row
        connection.execute("PRAGMA foreign_keys = ON")
        return connection

    @contextmanager
    def _connection(self):
        connection = self._connect()
        try:
            yield connection
            connection.commit()
        except Exception:
            connection.rollback()
            raise
        finally:
            connection.close()

    def initialize(self) -> None:
        self.database_path.parent.mkdir(parents=True, exist_ok=True)
        with self._connection() as connection:
            connection.executescript(
                """
                CREATE TABLE IF NOT EXISTS shrna_projects (
                    project_id TEXT PRIMARY KEY,
                    gene_symbol TEXT NOT NULL,
                    species TEXT NOT NULL,
                    received_date TEXT NOT NULL,
                    due_date TEXT NOT NULL,
                    project_folder TEXT NOT NULL,
                    folder_suffix TEXT NOT NULL,
                    design_json TEXT NOT NULL,
                    snapshot_json TEXT NOT NULL,
                    revision INTEGER NOT NULL,
                    status TEXT NOT NULL,
                    created_at TEXT NOT NULL,
                    updated_at TEXT NOT NULL
                );

                CREATE TABLE IF NOT EXISTS shrna_generated_artifacts (
                    artifact_id INTEGER PRIMARY KEY AUTOINCREMENT,
                    project_id TEXT NOT NULL,
                    artifact_json TEXT NOT NULL,
                    generated_at TEXT NOT NULL,
                    FOREIGN KEY(project_id) REFERENCES shrna_projects(project_id)
                );

                CREATE INDEX IF NOT EXISTS idx_shrna_artifacts_project
                ON shrna_generated_artifacts(project_id, artifact_id);
                """,
            )

    def create_project(
        self,
        *,
        project_id: str,
        gene_symbol: str,
        species: str,
        received_date: date,
        due_date: date,
        project_folder: Path,
        design: ShRNADesignVersion,
        snapshot: ShRNAProjectSnapshot,
        created_at: datetime,
        folder_suffix: str = "KD",
    ) -> None:
        if design.project_id != project_id or snapshot.project_id != project_id:
            raise ValueError("项目、设计版本和状态快照的 project_id 必须一致")
        if design.design_version_id != snapshot.active_design_version_id:
            raise ValueError("状态快照没有绑定当前 shRNA 设计版本")
        try:
            with self._connection() as connection:
                connection.execute(
                    """
                    INSERT INTO shrna_projects (
                        project_id, gene_symbol, species, received_date, due_date,
                        project_folder, folder_suffix, design_json, snapshot_json,
                        revision, status, created_at, updated_at
                    ) VALUES (?, ?, ?, ?, ?, ?, ?, ?, ?, ?, ?, ?, ?)
                    """,
                    (
                        project_id,
                        gene_symbol,
                        species,
                        received_date.isoformat(),
                        due_date.isoformat(),
                        str(Path(project_folder)),
                        folder_suffix,
                        dumps_record(design),
                        dumps_record(snapshot),
                        snapshot.revision,
                        snapshot.status,
                        created_at.isoformat(),
                        created_at.isoformat(),
                    ),
                )
        except sqlite3.IntegrityError as error:
            raise DuplicateProjectError(f"项目号已存在：{project_id}") from error

    def load_project(self, project_id: str) -> StoredShRNAProject:
        with self._connection() as connection:
            row = connection.execute(
                "SELECT * FROM shrna_projects WHERE project_id = ?",
                (project_id,),
            ).fetchone()
        if row is None:
            raise KeyError(project_id)
        design = loads_record(row["design_json"])
        snapshot = loads_record(row["snapshot_json"])
        if not isinstance(design, ShRNADesignVersion):
            raise TypeError("数据库中的 shRNA 设计记录类型无效")
        if not isinstance(snapshot, ShRNAProjectSnapshot):
            raise TypeError("数据库中的 shRNA 项目快照类型无效")
        return StoredShRNAProject(
            project_id=row["project_id"],
            gene_symbol=row["gene_symbol"],
            species=row["species"],
            received_date=date.fromisoformat(row["received_date"]),
            due_date=date.fromisoformat(row["due_date"]),
            project_folder=Path(row["project_folder"]),
            folder_suffix=row["folder_suffix"],
            design=design,
            snapshot=snapshot,
            created_at=datetime.fromisoformat(row["created_at"]),
            updated_at=datetime.fromisoformat(row["updated_at"]),
        )

    def list_projects(self) -> tuple[ShRNAProjectSummary, ...]:
        with self._connection() as connection:
            rows = connection.execute(
                """
                SELECT project_id, gene_symbol, species, status, received_date,
                       due_date, project_folder, folder_suffix, revision
                FROM shrna_projects
                ORDER BY created_at DESC, project_id
                """,
            ).fetchall()
        return tuple(
            ShRNAProjectSummary(
                project_id=row["project_id"],
                gene_symbol=row["gene_symbol"],
                species=row["species"],
                status=row["status"],
                received_date=date.fromisoformat(row["received_date"]),
                due_date=date.fromisoformat(row["due_date"]),
                project_folder=Path(row["project_folder"]),
                folder_suffix=row["folder_suffix"],
                revision=row["revision"],
            )
            for row in rows
        )

    def save_snapshot(
        self,
        project_id: str,
        snapshot: ShRNAProjectSnapshot,
        *,
        expected_revision: int,
        updated_at: datetime,
    ) -> None:
        if snapshot.project_id != project_id:
            raise ValueError("snapshot project_id 与待保存项目不一致")
        with self._connection() as connection:
            cursor = connection.execute(
                """
                UPDATE shrna_projects
                SET snapshot_json = ?, revision = ?, status = ?, updated_at = ?
                WHERE project_id = ? AND revision = ?
                """,
                (
                    dumps_record(snapshot),
                    snapshot.revision,
                    snapshot.status,
                    updated_at.isoformat(),
                    project_id,
                    expected_revision,
                ),
            )
            if cursor.rowcount != 1:
                raise StorageRevisionConflict(
                    f"项目 {project_id} 已被其他窗口更新，请重新载入",
                )

    def append_artifacts(
        self,
        project_id: str,
        artifacts: tuple[GeneratedArtifact, ...],
    ) -> None:
        with self._connection() as connection:
            exists = connection.execute(
                "SELECT 1 FROM shrna_projects WHERE project_id = ?",
                (project_id,),
            ).fetchone()
            if exists is None:
                raise KeyError(project_id)
            connection.executemany(
                """
                INSERT INTO shrna_generated_artifacts (
                    project_id, artifact_json, generated_at
                ) VALUES (?, ?, ?)
                """,
                tuple(
                    (
                        project_id,
                        dumps_record(artifact),
                        artifact.generated_at.isoformat(),
                    )
                    for artifact in artifacts
                ),
            )

    def list_artifacts(self, project_id: str) -> tuple[GeneratedArtifact, ...]:
        with self._connection() as connection:
            rows = connection.execute(
                """
                SELECT artifact_json
                FROM shrna_generated_artifacts
                WHERE project_id = ?
                ORDER BY artifact_id
                """,
                (project_id,),
            ).fetchall()
        artifacts = tuple(loads_record(row["artifact_json"]) for row in rows)
        if not all(isinstance(item, GeneratedArtifact) for item in artifacts):
            raise TypeError("数据库中的 shRNA 输出记录类型无效")
        return artifacts
