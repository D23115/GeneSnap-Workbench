"""Application orchestration shared by the desktop UI and acceptance tests."""

from __future__ import annotations

import hashlib
import shutil
from dataclasses import dataclass, replace
from datetime import date, datetime
from pathlib import Path
import re
from uuid import uuid4

from genesnap_workbench.domain.syn import (
    SYNAuditEvent,
    SYNDesignInput,
    SYNDesignVersion,
    SYNManualOverrideRecord,
    SYNProjectSnapshot,
)
from genesnap_workbench.domain.shrna import (
    ShRNACloneResultRecord,
    ShRNAAuditEvent,
    ShRNACandidate,
    ShRNADesignInput,
    ShRNAProjectSnapshot,
)
from genesnap_workbench.domain.expression import (
    ExpressionAuditEvent,
    ExpressionCloneResultRecord,
    ExpressionDesignInput,
    ExpressionProjectSnapshot,
)
from genesnap_workbench.domain.reporter import (
    ReporterAuditEvent,
    ReporterCloneResultRecord,
    ReporterDesignInput,
    ReporterProjectSnapshot,
)
from genesnap_workbench.domain.sequencing_submission import SequencingSubmissionRecord
from genesnap_workbench.domain.project_tracking import ProjectInterruptionRecord
from genesnap_workbench.sequencing.expression import (
    ExpressionCloneJudgmentStatus,
    judge_expression_read,
)
from genesnap_workbench.sequencing.shrna import (
    CloneJudgmentStatus,
    judge_shrna_read,
    match_shrna_sequence_files,
    read_sequence_file,
)
from genesnap_workbench.project_workflow.project_folders import (
    ProjectWorkspace,
    create_project_folder,
)
from genesnap_workbench.sequence_core.syn_design import (
    confirm_syn_design_warnings,
    create_syn_design,
)
from genesnap_workbench.sequence_core.shrna import (
    create_shrna_design,
    select_initial_candidates,
)
from genesnap_workbench.sequence_core.dna import normalize_dna
from genesnap_workbench.sequence_core.expression import (
    confirm_expression_design,
    create_expression_design,
)
from genesnap_workbench.sequence_core.reporter import (
    confirm_reporter_design,
    create_reporter_design,
)
from genesnap_workbench.storage.syn_repository import (
    DuplicateProjectError,
    SQLiteSYNProjectRepository,
    StoredSYNProject,
    SYNProjectSummary,
)
from genesnap_workbench.storage.shrna_repository import (
    SQLiteShRNAProjectRepository,
    StoredShRNAProject,
)
from genesnap_workbench.storage.expression_repository import (
    SQLiteExpressionProjectRepository,
    StoredExpressionProject,
)
from genesnap_workbench.storage.reporter_repository import (
    SQLiteReporterProjectRepository,
    StoredReporterProject,
)
from genesnap_workbench.storage.visibility import (
    LocalProjectVisibilityStore,
    ProjectVisibilityEvent,
)
from genesnap_workbench.storage.application_settings import LocalApplicationSettingsStore
from genesnap_workbench.template_engine.shrna_exports import (
    export_shrna_analysis_report,
    export_shrna_bundle,
)
from genesnap_workbench.template_engine.syn_exports import GeneratedArtifact, export_syn_bundle
from genesnap_workbench.template_engine.sequencing_forms import export_sequencing_form
from genesnap_workbench.template_engine.expression_exports import (
    export_expression_analysis_report,
    export_expression_bundle,
)
from genesnap_workbench.template_engine.reporter_exports import (
    export_reporter_analysis_report,
    export_reporter_bundle,
)
from genesnap_workbench.template_engine.workbook_templates import (
    ContactProfile,
    LocalContactProfileStore,
    LocalWorkbookTemplateStore,
    WorkbookTemplateInspection,
    WorkbookTemplateProfile,
)
from genesnap_workbench.vector_library.expression import (
    apply_expression_protocol,
    expression_rules_from_protocol,
)
from genesnap_workbench.vector_library.expression_profiles import (
    ExpressionProtocolProfileSummary,
    LocalExpressionProtocolStore,
)
from genesnap_workbench.vector_library.reporter_profiles import (
    LocalReporterProtocolStore,
    ReporterProtocolProfileSummary,
)
from genesnap_workbench.vector_library.models import (
    ExpressionVectorProtocol,
    ReporterVectorProtocol,
    SYNVectorProtocol,
    VectorRecord,
)
from genesnap_workbench.vector_library.reporter import apply_reporter_protocol
from genesnap_workbench.vector_library.starters import load_public_puc57_starter
from genesnap_workbench.vector_library.starters import load_public_plko1_puro_starter
from genesnap_workbench.project_workflow.business_calendar import ChinaBusinessCalendar


class DesignConfirmationRequired(ValueError):
    pass


def _latest_tracking_values(snapshot) -> tuple[str, str]:
    submissions = getattr(snapshot, "sequencing_submissions", ())
    if not submissions:
        return "", ""
    latest = submissions[-1]
    return latest.internal_submission_no, latest.vendor_order_no


def _molecular_tracking_values(snapshot) -> tuple[str, str, str]:
    return (
        getattr(snapshot, "internal_project_no", ""),
        getattr(snapshot, "primer_submission_no", ""),
        getattr(snapshot, "primer_vendor_order_no", ""),
    )


def _mark_latest_submission_analyzed(snapshot):
    submissions = getattr(snapshot, "sequencing_submissions", ())
    if not submissions:
        return snapshot
    latest = replace(submissions[-1], status="analyzed")
    return replace(snapshot, sequencing_submissions=submissions[:-1] + (latest,))


def _file_sha256(path: Path) -> str:
    digest = hashlib.sha256()
    with Path(path).open("rb") as handle:
        for chunk in iter(lambda: handle.read(1024 * 1024), b""):
            digest.update(chunk)
    return digest.hexdigest()


def _effective_due_date(snapshot, stored_due_date: date) -> date:
    return getattr(snapshot, "effective_due_date", None) or stored_due_date


def _validate_project_dates(received_date: date, due_date: date) -> None:
    if due_date < received_date:
        raise ValueError("标准完工日期不能早于接收日期")


@dataclass(frozen=True, slots=True)
class NewSYNProjectCommand:
    project_id: str
    target_name: str
    raw_sequence: str
    input_format: str
    linearization_site: str
    received_date: date
    due_date: date
    actor: str
    vector_sequence_confirmed: bool


@dataclass(frozen=True, slots=True)
class PreparedSYNProject:
    vector: VectorRecord
    protocol: SYNVectorProtocol
    design: SYNDesignVersion


@dataclass(frozen=True, slots=True)
class NewShRNAProjectCommand:
    project_id: str
    gene_symbol: str
    species: str
    cds_sequence: str
    candidates: tuple[ShRNACandidate, ...]
    target_count: int
    clones_per_target: int
    received_date: date
    due_date: date
    actor: str
    vector_sequence_confirmed: bool
    transcript_accession: str | None = None
    gene_id: str | None = None
    ccds_id: str | None = None
    primer_vendor_name: str = "标准"
    sequencing_vendor_name: str = "标准"
    primer_template_id: str | None = None
    sequencing_template_id: str | None = None


@dataclass(frozen=True, slots=True)
class NewExpressionProjectCommand:
    project_id: str
    gene_symbol: str
    species: str
    source_cds: str
    construct_lines: tuple[str, ...]
    received_date: date
    due_date: date
    actor: str
    vector: VectorRecord
    protocol: ExpressionVectorProtocol
    clones_per_construct: int = 5
    transcript_accession: str | None = None
    primer_vendor_name: str = "标准"
    sequencing_vendor_name: str = "标准"
    sequencing_method: str = "Nanopore"
    design_confirmation_reason: str | None = None
    primer_template_id: str | None = None
    sequencing_template_id: str | None = None
    gene_id: str | None = None


@dataclass(frozen=True, slots=True)
class NewReporterProjectCommand:
    project_id: str
    gene_symbol: str
    species: str
    promoter_sequence: str
    construct_lines: tuple[str, ...]
    mutation_definitions: tuple[str, ...]
    received_date: date
    due_date: date
    actor: str
    vector: VectorRecord
    protocol: ReporterVectorProtocol
    clones_per_construct: int = 5
    transcript_accession: str | None = None
    primer_vendor_name: str = "标准"
    sequencing_vendor_name: str = "标准"
    sequencing_method: str = "Nanopore"
    design_confirmation_reason: str | None = None
    primer_template_id: str | None = None
    sequencing_template_id: str | None = None
    gene_id: str | None = None


@dataclass(frozen=True, slots=True)
class UnifiedProjectSummary:
    project_id: str
    target_name: str
    project_category: str
    workflow_type: str
    status: str
    received_date: date
    due_date: date
    project_folder: Path
    folder_suffix: str
    revision: int
    design_summary: str
    usable_clone_names: tuple[str, ...] = ()
    latest_internal_submission_no: str = ""
    latest_vendor_order_no: str = ""
    interruption_type: str | None = None
    frozen_remaining_workdays: int | None = None
    is_manually_hidden: bool = False
    internal_project_no: str = ""
    primer_submission_no: str = ""
    primer_vendor_order_no: str = ""


@dataclass(frozen=True, slots=True)
class ShRNAAnalysisOutcome:
    project: StoredShRNAProject
    analysis_report: Path
    unmatched_files: tuple[Path, ...]
    ambiguous_files: tuple[Path, ...]


@dataclass(frozen=True, slots=True)
class AddOnSequencingPreview:
    project_id: str
    workflow_type: str
    round_no: int
    affected_owner_ids: tuple[str, ...]
    affected_owner_labels: tuple[str, ...]
    clones_per_owner: int
    sample_names: tuple[str, ...]
    reason: str


@dataclass(frozen=True, slots=True)
class ExpressionAnalysisOutcome:
    project: StoredExpressionProject
    analysis_report: Path
    unmatched_files: tuple[Path, ...]
    ambiguous_files: tuple[Path, ...]


@dataclass(frozen=True, slots=True)
class ReporterAnalysisOutcome:
    project: StoredReporterProject
    analysis_report: Path
    unmatched_files: tuple[Path, ...]
    ambiguous_files: tuple[Path, ...]


class GeneSnapApplicationService:
    MOLECULAR_TRANSITIONS = {
        "mark_primers_ordered": ("design_completed", "primers_ordered"),
        "mark_primers_arrived": ("primers_ordered", "primers_arrived"),
        "start_cloning": ("primers_arrived", "cloning_in_progress"),
        "mark_sent_for_sequencing": (
            "cloning_in_progress",
            "sequencing_in_progress",
        ),
        "start_plasmid_prep": (
            "analysis_completed",
            "plasmid_prep_in_progress",
        ),
        "complete_plasmid_prep": (
            "plasmid_prep_in_progress",
            "plasmid_prep_completed",
        ),
        "complete_project": ("plasmid_prep_completed", "project_completed"),
    }

    def __init__(self, data_root: Path) -> None:
        self.data_root = Path(data_root)
        self.repository = SQLiteSYNProjectRepository(
            self.data_root / "genesnap_workbench.db",
        )
        self.repository.initialize()
        self.shrna_repository = SQLiteShRNAProjectRepository(
            self.data_root / "genesnap_workbench.db",
        )
        self.shrna_repository.initialize()
        self.expression_repository = SQLiteExpressionProjectRepository(
            self.data_root / "genesnap_workbench.db",
        )
        self.expression_repository.initialize()
        self.reporter_repository = SQLiteReporterProjectRepository(
            self.data_root / "genesnap_workbench.db",
        )
        self.reporter_repository.initialize()
        self.expression_profile_store = LocalExpressionProtocolStore(
            self.data_root / "vector_protocols" / "expression",
        )
        self.reporter_profile_store = LocalReporterProtocolStore(
            self.data_root / "vector_protocols" / "reporter",
        )
        self.workbook_template_store = LocalWorkbookTemplateStore(
            self.data_root / "workbook_templates",
        )
        self.contact_profile_store = LocalContactProfileStore(
            self.data_root / "contact_profile.json",
        )
        self.visibility_store = LocalProjectVisibilityStore(
            self.data_root / "project_visibility.json",
        )
        self.application_settings_store = LocalApplicationSettingsStore(
            self.data_root / "application_settings.json",
        )

    @property
    def projects_root(self) -> Path:
        configured = self.application_settings_store.load().get("projects_root")
        if isinstance(configured, str) and configured.strip():
            return Path(configured)
        return self.data_root / "projects"

    @property
    def has_custom_projects_root(self) -> bool:
        configured = self.application_settings_store.load().get("projects_root")
        return isinstance(configured, str) and bool(configured.strip())

    def set_projects_root(self, path: Path) -> Path:
        chosen = Path(path).expanduser()
        if not str(chosen).strip():
            raise ValueError("项目保存位置不能为空")
        chosen.mkdir(parents=True, exist_ok=True)
        resolved = chosen.resolve()
        self.application_settings_store.update(projects_root=str(resolved))
        return resolved

    def has_accepted_broad_terms(self, version: str = "2025-12-08") -> bool:
        return self.application_settings_store.load().get("broad_terms_version") == version

    def accept_broad_terms(self, version: str = "2025-12-08") -> None:
        self.application_settings_store.update(broad_terms_version=version)

    def list_projects(self) -> tuple[SYNProjectSummary, ...]:
        return self.repository.list_projects()

    def save_expression_profile(
        self,
        vector: VectorRecord,
        protocol: ExpressionVectorProtocol,
    ) -> ExpressionProtocolProfileSummary:
        return self.expression_profile_store.save_profile(vector, protocol)

    def list_expression_profiles(self) -> tuple[ExpressionProtocolProfileSummary, ...]:
        return self.expression_profile_store.list_profiles()

    def load_expression_profile(
        self,
        profile_id: str,
    ) -> tuple[VectorRecord, ExpressionVectorProtocol]:
        return self.expression_profile_store.load_profile(profile_id)

    def save_reporter_profile(
        self,
        vector: VectorRecord,
        protocol: ReporterVectorProtocol,
    ) -> ReporterProtocolProfileSummary:
        return self.reporter_profile_store.save_profile(vector, protocol)

    def list_reporter_profiles(self) -> tuple[ReporterProtocolProfileSummary, ...]:
        return self.reporter_profile_store.list_profiles()

    def load_reporter_profile(
        self,
        profile_id: str,
    ) -> tuple[VectorRecord, ReporterVectorProtocol]:
        return self.reporter_profile_store.load_profile(profile_id)

    def save_workbook_template(
        self,
        source_path: Path,
        *,
        display_name: str,
        inspected: WorkbookTemplateInspection,
    ) -> WorkbookTemplateProfile:
        return self.workbook_template_store.save_import(
            source_path,
            display_name=display_name,
            inspected=inspected,
        )

    def list_workbook_templates(
        self,
        kind: str | None = None,
    ) -> tuple[WorkbookTemplateProfile, ...]:
        return self.workbook_template_store.list_profiles(kind)

    def load_contact_profile(self) -> ContactProfile:
        return self.contact_profile_store.load()

    def save_contact_profile(self, profile: ContactProfile) -> None:
        self.contact_profile_store.save(profile)

    def set_project_hidden(
        self,
        project_id: str,
        *,
        hidden: bool,
        reason: str,
        actor: str,
        occurred_at: datetime,
    ) -> None:
        if project_id not in {item.project_id for item in self.list_all_projects()}:
            raise KeyError(f"未找到项目：{project_id}")
        self.visibility_store.append(
            ProjectVisibilityEvent(
                project_id=project_id,
                hidden=hidden,
                reason=reason.strip(),
                actor=actor,
                occurred_at=occurred_at,
            ),
        )

    def list_all_projects(self) -> tuple[UnifiedProjectSummary, ...]:
        summaries: list[tuple[datetime, UnifiedProjectSummary]] = []
        for item in self.repository.list_projects():
            stored = self.repository.load_project(item.project_id)
            latest_confirmations = {}
            for record in stored.snapshot.sequencing_confirmations:
                latest_confirmations[record.clone_id] = record
            usable = tuple(
                record.clone_id
                for record in latest_confirmations.values()
                if record.result.value == "correct"
            )
            route = (
                "模块化"
                if stored.design.module_plan.route.value == "modular"
                else "单池"
            )
            summaries.append(
                (
                    stored.created_at,
                    UnifiedProjectSummary(
                        project_id=item.project_id,
                        target_name=item.target_name,
                        project_category="合成/组装类",
                        workflow_type="de_novo_gene_synthesis",
                        status=item.status,
                        received_date=item.received_date,
                        due_date=_effective_due_date(stored.snapshot, item.due_date),
                        project_folder=item.project_folder,
                        folder_suffix=item.folder_suffix,
                        revision=item.revision,
                        design_summary=(
                            f"SYN：{len(stored.design.final_sequence)} bp / {route} / "
                            f"{len(stored.design.oligos)} oligos"
                        ),
                        usable_clone_names=usable,
                        latest_internal_submission_no=_latest_tracking_values(stored.snapshot)[0],
                        latest_vendor_order_no=_latest_tracking_values(stored.snapshot)[1],
                        interruption_type=getattr(stored.snapshot, "interruption_type", None),
                        frozen_remaining_workdays=getattr(
                            stored.snapshot,
                            "frozen_remaining_workdays",
                            None,
                        ),
                        is_manually_hidden=self.visibility_store.is_hidden(item.project_id),
                    ),
                ),
            )
        for item in self.shrna_repository.list_projects():
            stored = self.shrna_repository.load_project(item.project_id)
            latest_by_clone = {}
            for record in stored.snapshot.clone_results:
                latest_by_clone[record.clone_name] = record
            usable = tuple(
                name
                for name, record in latest_by_clone.items()
                if record.status == "pass" or record.manually_confirmed_usable
            )
            summaries.append(
                (
                    stored.created_at,
                    UnifiedProjectSummary(
                        project_id=item.project_id,
                        target_name=item.gene_symbol,
                        project_category="沉默/敲低类",
                        workflow_type="shrna_knockdown",
                        status=item.status,
                        received_date=item.received_date,
                        due_date=_effective_due_date(stored.snapshot, item.due_date),
                        project_folder=item.project_folder,
                        folder_suffix=item.folder_suffix,
                        revision=item.revision,
                        design_summary=(
                            f"shRNA：{stored.design.target_count} targets / "
                            f"{stored.design.target_count * stored.design.clones_per_target} clones"
                        ),
                        usable_clone_names=usable,
                        latest_internal_submission_no=_latest_tracking_values(stored.snapshot)[0],
                        latest_vendor_order_no=_latest_tracking_values(stored.snapshot)[1],
                        internal_project_no=_molecular_tracking_values(stored.snapshot)[0],
                        primer_submission_no=_molecular_tracking_values(stored.snapshot)[1],
                        primer_vendor_order_no=_molecular_tracking_values(stored.snapshot)[2],
                        interruption_type=getattr(stored.snapshot, "interruption_type", None),
                        frozen_remaining_workdays=getattr(
                            stored.snapshot,
                            "frozen_remaining_workdays",
                            None,
                        ),
                        is_manually_hidden=self.visibility_store.is_hidden(item.project_id),
                    ),
                ),
            )
        for item in self.expression_repository.list_projects():
            stored = self.expression_repository.load_project(item.project_id)
            prefix = f"{stored.design.gene_symbol}-"
            construct_names = tuple(
                construct.construct_name.removeprefix(prefix)
                for construct in stored.design.constructs
            )
            latest_by_clone = {}
            for record in stored.snapshot.clone_results:
                latest_by_clone[record.clone_name] = record
            usable = tuple(
                name
                for name, record in latest_by_clone.items()
                if record.status == "pass" or record.manually_confirmed_usable
            )
            summaries.append(
                (
                    stored.created_at,
                    UnifiedProjectSummary(
                        project_id=item.project_id,
                        target_name=item.gene_symbol,
                        project_category="表达类",
                        workflow_type="expression",
                        status=item.status,
                        received_date=item.received_date,
                        due_date=_effective_due_date(stored.snapshot, item.due_date),
                        project_folder=item.project_folder,
                        folder_suffix=item.folder_suffix,
                        revision=item.revision,
                        design_summary=(
                            f"表达类：{' + '.join(construct_names)}，"
                            f"共 {len(construct_names)} 个构建"
                        ),
                        usable_clone_names=usable,
                        latest_internal_submission_no=_latest_tracking_values(stored.snapshot)[0],
                        latest_vendor_order_no=_latest_tracking_values(stored.snapshot)[1],
                        internal_project_no=_molecular_tracking_values(stored.snapshot)[0],
                        primer_submission_no=_molecular_tracking_values(stored.snapshot)[1],
                        primer_vendor_order_no=_molecular_tracking_values(stored.snapshot)[2],
                        interruption_type=getattr(stored.snapshot, "interruption_type", None),
                        frozen_remaining_workdays=getattr(
                            stored.snapshot,
                            "frozen_remaining_workdays",
                            None,
                        ),
                        is_manually_hidden=self.visibility_store.is_hidden(item.project_id),
                    ),
                ),
            )
        for item in self.reporter_repository.list_projects():
            stored = self.reporter_repository.load_project(item.project_id)
            full_length = max(
                construct.retained_promoter_length
                for construct in stored.design.constructs
            )
            construct_labels = []
            for construct in stored.design.constructs:
                if construct.mutation_names:
                    prefix = (
                        ""
                        if construct.retained_promoter_length == full_length
                        else f"P{construct.retained_promoter_length}+"
                    )
                    construct_labels.append(prefix + "+".join(construct.mutation_names))
                elif construct.retained_promoter_length == full_length:
                    construct_labels.append("WT")
                else:
                    construct_labels.append(f"P{construct.retained_promoter_length}")
            latest = {}
            for record in stored.snapshot.clone_results:
                latest[record.clone_name] = record
            usable = tuple(
                name
                for name, record in latest.items()
                if record.status == "pass" or record.manually_confirmed_usable
            )
            summaries.append(
                (
                    stored.created_at,
                    UnifiedProjectSummary(
                        project_id=item.project_id,
                        target_name=item.gene_symbol,
                        project_category="报告/检测类",
                        workflow_type="promoter_luciferase_reporter",
                        status=item.status,
                        received_date=item.received_date,
                        due_date=_effective_due_date(stored.snapshot, item.due_date),
                        project_folder=item.project_folder,
                        folder_suffix=item.folder_suffix,
                        revision=item.revision,
                        design_summary=(
                            f"报告类：{' + '.join(construct_labels)}，"
                            f"共 {len(construct_labels)} 个构建"
                        ),
                        usable_clone_names=usable,
                        latest_internal_submission_no=_latest_tracking_values(stored.snapshot)[0],
                        latest_vendor_order_no=_latest_tracking_values(stored.snapshot)[1],
                        internal_project_no=_molecular_tracking_values(stored.snapshot)[0],
                        primer_submission_no=_molecular_tracking_values(stored.snapshot)[1],
                        primer_vendor_order_no=_molecular_tracking_values(stored.snapshot)[2],
                        interruption_type=getattr(stored.snapshot, "interruption_type", None),
                        frozen_remaining_workdays=getattr(
                            stored.snapshot,
                            "frozen_remaining_workdays",
                            None,
                        ),
                        is_manually_hidden=self.visibility_store.is_hidden(item.project_id),
                    ),
                ),
            )
        summaries.sort(key=lambda pair: (pair[0], pair[1].project_id), reverse=True)
        return tuple(item for _, item in summaries)

    def load_any_project(
        self,
        project_id: str,
        workflow_type: str,
    ) -> (
        StoredSYNProject
        | StoredShRNAProject
        | StoredExpressionProject
        | StoredReporterProject
    ):
        if workflow_type == "promoter_luciferase_reporter":
            return self.load_reporter_project(project_id)
        if workflow_type == "expression":
            return self.load_expression_project(project_id)
        if workflow_type == "shrna_knockdown":
            return self.load_shrna_project(project_id)
        if workflow_type == "de_novo_gene_synthesis":
            return self.load_project(project_id)
        raise KeyError(f"未知工作流：{workflow_type}")

    def load_project(self, project_id: str) -> StoredSYNProject:
        return self.repository.load_project(project_id)

    def load_shrna_project(self, project_id: str) -> StoredShRNAProject:
        return self.shrna_repository.load_project(project_id)

    def load_expression_project(self, project_id: str) -> StoredExpressionProject:
        return self.expression_repository.load_project(project_id)

    def load_reporter_project(self, project_id: str) -> StoredReporterProject:
        return self.reporter_repository.load_project(project_id)

    def transition_molecular_project(
        self,
        project_id: str,
        *,
        workflow_type: str,
        action: str,
        actor: str,
        occurred_at: datetime,
        selected_clone_names: tuple[str, ...] = (),
        note: str | None = None,
        internal_project_no: str = "",
        primer_submission_no: str = "",
        primer_vendor_order_no: str = "",
        internal_submission_no: str = "",
        vendor_order_no: str = "",
    ) -> StoredExpressionProject | StoredShRNAProject | StoredReporterProject:
        try:
            expected_status, next_status = self.MOLECULAR_TRANSITIONS[action]
        except KeyError as error:
            raise ValueError(f"未知项目操作：{action}") from error
        if workflow_type == "expression":
            repository = self.expression_repository
            stored = repository.load_project(project_id)
            event_type = ExpressionAuditEvent
        elif workflow_type == "shrna_knockdown":
            repository = self.shrna_repository
            stored = repository.load_project(project_id)
            event_type = ShRNAAuditEvent
        elif workflow_type == "promoter_luciferase_reporter":
            repository = self.reporter_repository
            stored = repository.load_project(project_id)
            event_type = ReporterAuditEvent
        else:
            raise ValueError(f"该操作不支持工作流：{workflow_type}")
        if stored.snapshot.status != expected_status:
            raise ValueError(
                f"当前状态 {stored.snapshot.status} 不能执行 {action}；"
                f"需要 {expected_status}",
            )

        snapshot = stored.snapshot
        if action == "mark_primers_ordered":
            snapshot = replace(
                snapshot,
                internal_project_no=internal_project_no.strip(),
                primer_submission_no=primer_submission_no.strip(),
                primer_vendor_order_no=primer_vendor_order_no.strip(),
            )
        elif action == "mark_sent_for_sequencing":
            if workflow_type in {"expression", "promoter_luciferase_reporter"}:
                sample_names = tuple(
                    f"{construct.construct_name}-{clone_no}"
                    for construct in stored.design.constructs
                    for clone_no in range(1, snapshot.clones_per_construct + 1)
                )
            else:
                sample_names = tuple(
                    clone_name
                    for target in stored.design.targets
                    for clone_name in target.clone_names
                )
            submission = SequencingSubmissionRecord(
                submission_id=f"submission-{uuid4()}",
                round_no=len(snapshot.sequencing_submissions) + 1,
                submission_kind=("initial" if not snapshot.sequencing_submissions else "post_rework"),
                created_at=occurred_at,
                sent_at=occurred_at,
                sample_names=sample_names,
                internal_submission_no=internal_submission_no.strip(),
                vendor_order_no=vendor_order_no.strip(),
                note=note,
            )
            snapshot = replace(
                snapshot,
                sequencing_submissions=snapshot.sequencing_submissions + (submission,),
            )
        elif action == "start_plasmid_prep":
            selected = tuple(dict.fromkeys(selected_clone_names))
            if not selected:
                raise ValueError("开始质粒抽提前必须选择可用克隆")
            latest = {}
            for record in snapshot.clone_results:
                latest[record.clone_name] = record
            usable = {
                name: record
                for name, record in latest.items()
                if record.status == "pass" or record.manually_confirmed_usable
            }
            invalid = tuple(name for name in selected if name not in usable)
            if invalid:
                raise ValueError(f"选择中包含不可用克隆：{', '.join(invalid)}")
            if workflow_type in {"expression", "promoter_luciferase_reporter"}:
                required_owners = {
                    construct.construct_id for construct in stored.design.constructs
                }
                selected_owners = {usable[name].construct_id for name in selected}
            else:
                required_owners = {target.target_id for target in stored.design.targets}
                selected_owners = {usable[name].target_id for name in selected}
            if selected_owners != required_owners:
                raise ValueError("每个构建或 target 至少需要选择 1 个可用克隆")
            snapshot = replace(snapshot, selected_prep_clone_names=selected)
        elif action == "complete_plasmid_prep":
            snapshot = replace(snapshot, plasmid_prep_completed_at=occurred_at)
        elif action == "complete_project":
            snapshot = replace(snapshot, actual_completed_at=occurred_at)

        snapshot = replace(snapshot, status=next_status)
        snapshot = snapshot.append_status_event(
            event_type(
                event_id=f"{action}-{uuid4()}",
                event_type=action,
                occurred_at=occurred_at,
                actor=actor,
                from_status=expected_status,
                to_status=next_status,
                note=note,
                source="user",
            ),
        )
        repository.save_snapshot(
            project_id,
            snapshot,
            expected_revision=stored.snapshot.revision,
            updated_at=occurred_at,
        )
        return repository.load_project(project_id)

    def update_molecular_tracking_numbers(
        self,
        project_id: str,
        *,
        workflow_type: str,
        internal_project_no: str,
        primer_submission_no: str,
        primer_vendor_order_no: str,
        actor: str,
        occurred_at: datetime,
        note: str,
    ) -> StoredExpressionProject | StoredShRNAProject | StoredReporterProject:
        clean_note = note.strip()
        if not clean_note:
            raise ValueError("修改项目/引物编号时必须填写修改说明")
        if workflow_type == "expression":
            repository = self.expression_repository
            stored = repository.load_project(project_id)
            event_type = ExpressionAuditEvent
        elif workflow_type == "shrna_knockdown":
            repository = self.shrna_repository
            stored = repository.load_project(project_id)
            event_type = ShRNAAuditEvent
        elif workflow_type == "promoter_luciferase_reporter":
            repository = self.reporter_repository
            stored = repository.load_project(project_id)
            event_type = ReporterAuditEvent
        else:
            raise ValueError(f"该操作不支持工作流：{workflow_type}")
        old_internal_project_no = stored.snapshot.internal_project_no
        old_primer_submission_no = stored.snapshot.primer_submission_no
        old_primer_vendor_order_no = stored.snapshot.primer_vendor_order_no
        new_internal_project_no = internal_project_no.strip()
        new_primer_submission_no = primer_submission_no.strip()
        new_primer_vendor_order_no = primer_vendor_order_no.strip()
        snapshot = replace(
            stored.snapshot,
            internal_project_no=new_internal_project_no,
            primer_submission_no=new_primer_submission_no,
            primer_vendor_order_no=new_primer_vendor_order_no,
        )
        audit_note = (
            f"{clean_note}；"
            f"内部编号：{old_internal_project_no or '-'} -> {new_internal_project_no or '-'}；"
            f"引物送单号：{old_primer_submission_no or '-'} -> "
            f"{new_primer_submission_no or '-'}；"
            f"引物订单号：{old_primer_vendor_order_no or '-'} -> "
            f"{new_primer_vendor_order_no or '-'}"
        )
        snapshot = snapshot.append_status_event(
            event_type(
                event_id=f"update-molecular-tracking-numbers-{uuid4()}",
                event_type="update_molecular_tracking_numbers",
                occurred_at=occurred_at,
                actor=actor,
                from_status=stored.snapshot.status,
                to_status=stored.snapshot.status,
                note=audit_note,
                source="user",
            ),
        )
        repository.save_snapshot(
            project_id,
            snapshot,
            expected_revision=stored.snapshot.revision,
            updated_at=occurred_at,
        )
        return repository.load_project(project_id)

    def update_latest_sequencing_tracking(
        self,
        project_id: str,
        *,
        workflow_type: str,
        internal_submission_no: str,
        vendor_order_no: str,
        actor: str,
        occurred_at: datetime,
        note: str,
    ) -> StoredExpressionProject | StoredShRNAProject | StoredReporterProject:
        if not note.strip():
            raise ValueError("修改送测编号或订单号时必须填写说明")
        if workflow_type == "expression":
            repository = self.expression_repository
            stored = repository.load_project(project_id)
            event_type = ExpressionAuditEvent
        elif workflow_type == "shrna_knockdown":
            repository = self.shrna_repository
            stored = repository.load_project(project_id)
            event_type = ShRNAAuditEvent
        elif workflow_type == "promoter_luciferase_reporter":
            repository = self.reporter_repository
            stored = repository.load_project(project_id)
            event_type = ReporterAuditEvent
        else:
            raise ValueError(f"该操作不支持工作流：{workflow_type}")
        if not stored.snapshot.sequencing_submissions:
            raise ValueError("项目还没有送测记录")
        latest = stored.snapshot.sequencing_submissions[-1]
        corrected = replace(
            latest,
            internal_submission_no=internal_submission_no.strip(),
            vendor_order_no=vendor_order_no.strip(),
        )
        snapshot = replace(
            stored.snapshot,
            sequencing_submissions=stored.snapshot.sequencing_submissions[:-1] + (corrected,),
        )
        snapshot = snapshot.append_status_event(
            event_type(
                event_id=f"update-sequencing-tracking-{uuid4()}",
                event_type="update_sequencing_tracking",
                occurred_at=occurred_at,
                actor=actor,
                from_status=stored.snapshot.status,
                to_status=stored.snapshot.status,
                note=note.strip(),
                source="user",
            ),
        )
        repository.save_snapshot(
            project_id,
            snapshot,
            expected_revision=stored.snapshot.revision,
            updated_at=occurred_at,
        )
        return repository.load_project(project_id)

    def adjust_molecular_due_date(
        self,
        project_id: str,
        *,
        workflow_type: str,
        new_due_date: date,
        note: str,
        actor: str,
        occurred_at: datetime,
    ) -> StoredExpressionProject | StoredShRNAProject | StoredReporterProject:
        clean_note = note.strip()
        if not clean_note:
            raise ValueError("修正标准完工日期时必须填写原因")
        if workflow_type == "expression":
            repository = self.expression_repository
            stored = repository.load_project(project_id)
            event_type = ExpressionAuditEvent
        elif workflow_type == "shrna_knockdown":
            repository = self.shrna_repository
            stored = repository.load_project(project_id)
            event_type = ShRNAAuditEvent
        elif workflow_type == "promoter_luciferase_reporter":
            repository = self.reporter_repository
            stored = repository.load_project(project_id)
            event_type = ReporterAuditEvent
        else:
            raise ValueError(f"该操作不支持工作流：{workflow_type}")
        if stored.snapshot.status == "abnormal_or_paused":
            raise ValueError("暂停或异常期间不能修正标准完工日期")
        if stored.snapshot.status == "project_completed":
            raise ValueError("已完成项目需要先进入完成后修正流程")
        _validate_project_dates(stored.received_date, new_due_date)
        previous_due_date = _effective_due_date(stored.snapshot, stored.due_date)
        if new_due_date == previous_due_date:
            raise ValueError("新标准完工日期与当前日期相同")

        snapshot = replace(
            stored.snapshot,
            effective_due_date=new_due_date,
        )
        snapshot = snapshot.append_status_event(
            event_type(
                event_id=f"adjust-due-date-{uuid4()}",
                event_type="adjust_due_date",
                occurred_at=occurred_at,
                actor=actor,
                from_status=stored.snapshot.status,
                to_status=stored.snapshot.status,
                note=(
                    f"标准完工日期：{previous_due_date.isoformat()} -> "
                    f"{new_due_date.isoformat()}；{clean_note}"
                ),
                source="user",
            ),
        )
        repository.save_snapshot(
            project_id,
            snapshot,
            expected_revision=stored.snapshot.revision,
            updated_at=occurred_at,
        )
        return repository.load_project(project_id)

    def preview_addon_sequencing(
        self,
        project_id: str,
        *,
        workflow_type: str,
        clones_per_owner: int | None = None,
    ) -> AddOnSequencingPreview:
        if workflow_type == "expression":
            stored = self.expression_repository.load_project(project_id)
        elif workflow_type == "shrna_knockdown":
            stored = self.shrna_repository.load_project(project_id)
        elif workflow_type == "promoter_luciferase_reporter":
            stored = self.reporter_repository.load_project(project_id)
        else:
            raise ValueError(f"该操作不支持工作流：{workflow_type}")
        snapshot = stored.snapshot
        if snapshot.status != "analysis_completed":
            raise ValueError("必须先完成本轮全部测序分析和人工复核，才能生成加测")
        if not snapshot.sequencing_submissions:
            raise ValueError("项目没有可追溯的送测记录")
        if snapshot.sequencing_submissions[-1].status != "analyzed":
            raise ValueError("当前送测轮次尚未完成分析")

        latest_by_clone = {}
        for record in snapshot.clone_results:
            latest_by_clone[record.clone_name] = record
        unresolved = tuple(
            record.clone_name
            for record in latest_by_clone.values()
            if record.status == "warning"
            and not record.manually_confirmed_usable
            and getattr(record, "manual_review_status", None) is None
        )
        if unresolved:
            raise ValueError("仍有 WARNING 克隆未完成人工复核，不能生成加测")

        previous_sample_names = tuple(
            sample_name
            for submission in snapshot.sequencing_submissions
            for sample_name in submission.sample_names
        )
        affected_ids: list[str] = []
        affected_labels: list[str] = []
        generated_names: list[str] = []
        default_count = 10 if workflow_type == "shrna_knockdown" else 5
        count = default_count if clones_per_owner is None else clones_per_owner
        if not 1 <= count <= 96:
            raise ValueError("每个 target/构建的加测克隆数必须在 1 到 96 之间")

        if workflow_type == "shrna_knockdown":
            owners = tuple(
                (
                    target.target_id,
                    f"shRNA-{target.target_no}",
                    target.clone_names[0].rsplit("-", 1)[0],
                )
                for target in stored.design.targets
            )
            owner_field = "target_id"
        else:
            owners = tuple(
                (
                    construct.construct_id,
                    construct.construct_name,
                    construct.construct_name,
                )
                for construct in stored.design.constructs
            )
            owner_field = "construct_id"

        for owner_id, owner_label, sample_prefix in owners:
            owner_records = tuple(
                record
                for record in latest_by_clone.values()
                if getattr(record, owner_field) == owner_id
            )
            if any(
                record.status == "pass" or record.manually_confirmed_usable
                for record in owner_records
            ):
                continue
            pattern = re.compile(rf"^{re.escape(sample_prefix)}-(\d+)$")
            existing_numbers = tuple(
                int(match.group(1))
                for sample_name in previous_sample_names
                if (match := pattern.fullmatch(sample_name))
            )
            start_no = max(existing_numbers, default=0) + 1
            affected_ids.append(owner_id)
            affected_labels.append(owner_label)
            generated_names.extend(
                f"{sample_prefix}-{clone_no}"
                for clone_no in range(start_no, start_no + count)
            )
        if not affected_ids:
            raise ValueError("每个 target/构建都已有可用克隆，不需要加测")
        return AddOnSequencingPreview(
            project_id=project_id,
            workflow_type=workflow_type,
            round_no=len(snapshot.sequencing_submissions) + 1,
            affected_owner_ids=tuple(affected_ids),
            affected_owner_labels=tuple(affected_labels),
            clones_per_owner=count,
            sample_names=tuple(generated_names),
            reason="当前轮次分析完成后，部分 target/构建仍无可用克隆",
        )

    def confirm_addon_sequencing(
        self,
        preview: AddOnSequencingPreview,
        *,
        actor: str,
        occurred_at: datetime,
        sequencing_vendor_name: str = "标准",
        sequencing_template_id: str | None = None,
        internal_submission_no: str = "",
        vendor_order_no: str = "",
    ) -> StoredExpressionProject | StoredShRNAProject | StoredReporterProject:
        fresh = self.preview_addon_sequencing(
            preview.project_id,
            workflow_type=preview.workflow_type,
            clones_per_owner=preview.clones_per_owner,
        )
        if fresh != preview:
            raise ValueError("项目数据已变化，请重新预览加测方案")
        if preview.workflow_type == "expression":
            repository = self.expression_repository
            stored = repository.load_project(preview.project_id)
            event_type = ExpressionAuditEvent
            sequencing_method = "Nanopore"
            primer_name = ""
        elif preview.workflow_type == "shrna_knockdown":
            repository = self.shrna_repository
            stored = repository.load_project(preview.project_id)
            event_type = ShRNAAuditEvent
            sequencing_method = "Sanger"
            primer_name = "U6"
        elif preview.workflow_type == "promoter_luciferase_reporter":
            repository = self.reporter_repository
            stored = repository.load_project(preview.project_id)
            event_type = ReporterAuditEvent
            sequencing_method = "Nanopore"
            primer_name = ""
        else:
            raise ValueError(f"该操作不支持工作流：{preview.workflow_type}")

        safe_vendor = re.sub(r'[<>:"/\\|?*]+', "_", sequencing_vendor_name).strip() or "标准"
        safe_gene = re.sub(r'[<>:"/\\|?*]+', "_", stored.design.gene_symbol).strip() or "gene"
        output_path = (
            stored.project_folder
            / "02_orders"
            / (
                f"{occurred_at:%Y%m%d}-{safe_vendor}测序表-"
                f"{safe_gene}_加测{preview.round_no}.xlsx"
            )
        )
        if output_path.exists():
            output_path = output_path.with_name(
                f"{output_path.stem}-{occurred_at:%H%M%S}{output_path.suffix}",
            )
        records = tuple(
            {
                "sample_name": sample_name,
                "primer_name": primer_name,
                "gene_symbol": stored.design.gene_symbol,
                "clone_no": int(sample_name.rsplit("-", 1)[-1]),
                "method": sequencing_method,
                "note": "加测",
            }
            for sample_name in preview.sample_names
        )
        export_sequencing_form(
            records,
            output_path,
            template_store=self.workbook_template_store,
            contact_profile=self.contact_profile_store.load(),
            template_id=sequencing_template_id,
        )
        submission = SequencingSubmissionRecord(
            submission_id=f"submission-{uuid4()}",
            round_no=preview.round_no,
            submission_kind="add_on",
            created_at=occurred_at,
            sent_at=occurred_at,
            sample_names=preview.sample_names,
            internal_submission_no=internal_submission_no.strip(),
            vendor_order_no=vendor_order_no.strip(),
            template_id=sequencing_template_id,
            form_path=str(output_path),
            status="sent",
            note=preview.reason,
        )
        snapshot = replace(
            stored.snapshot,
            status="add_on_in_progress",
            sequencing_submissions=stored.snapshot.sequencing_submissions + (submission,),
        )
        snapshot = snapshot.append_status_event(
            event_type(
                event_id=f"confirm-addon-sequencing-{uuid4()}",
                event_type="confirm_addon_sequencing",
                occurred_at=occurred_at,
                actor=actor,
                from_status=stored.snapshot.status,
                to_status="add_on_in_progress",
                note=(
                    f"{preview.reason}；{', '.join(preview.affected_owner_labels)}；"
                    f"新增 {len(preview.sample_names)} 个克隆"
                ),
                source="user",
            ),
        )
        repository.save_snapshot(
            preview.project_id,
            snapshot,
            expected_revision=stored.snapshot.revision,
            updated_at=occurred_at,
        )
        repository.append_artifacts(
            preview.project_id,
            (
                GeneratedArtifact(
                    artifact_type="addon_sequencing_order_xlsx",
                    design_version_id=stored.design.design_version_id,
                    generated_at=occurred_at,
                    path=output_path,
                    content_sha256=_file_sha256(output_path),
                ),
            ),
        )
        return repository.load_project(preview.project_id)

    def start_molecular_rework(
        self,
        project_id: str,
        *,
        workflow_type: str,
        actor: str,
        occurred_at: datetime,
        note: str,
    ) -> StoredExpressionProject | StoredShRNAProject | StoredReporterProject:
        if not note.strip():
            raise ValueError("重新连接/转化必须填写原因")
        failed = self.preview_addon_sequencing(
            project_id,
            workflow_type=workflow_type,
        )
        if workflow_type == "expression":
            repository = self.expression_repository
            stored = repository.load_project(project_id)
            event_type = ExpressionAuditEvent
        elif workflow_type == "shrna_knockdown":
            repository = self.shrna_repository
            stored = repository.load_project(project_id)
            event_type = ShRNAAuditEvent
        elif workflow_type == "promoter_luciferase_reporter":
            repository = self.reporter_repository
            stored = repository.load_project(project_id)
            event_type = ReporterAuditEvent
        else:
            raise ValueError(f"该操作不支持工作流：{workflow_type}")
        snapshot = replace(
            stored.snapshot,
            status="rework_in_progress",
            experiment_attempt_no=stored.snapshot.experiment_attempt_no + 1,
            rework_owner_ids=failed.affected_owner_ids,
        )
        snapshot = snapshot.append_status_event(
            event_type(
                event_id=f"start-rework-{uuid4()}",
                event_type="restart_cloning",
                occurred_at=occurred_at,
                actor=actor,
                from_status=stored.snapshot.status,
                to_status="rework_in_progress",
                note=f"{note.strip()}；对象：{', '.join(failed.affected_owner_labels)}",
                source="user",
            ),
        )
        repository.save_snapshot(
            project_id,
            snapshot,
            expected_revision=stored.snapshot.revision,
            updated_at=occurred_at,
        )
        return repository.load_project(project_id)

    def preview_rework_submission(
        self,
        project_id: str,
        *,
        workflow_type: str,
        clones_per_owner: int = 5,
    ) -> AddOnSequencingPreview:
        if not 1 <= clones_per_owner <= 96:
            raise ValueError("每个重做对象的送测克隆数必须在 1 到 96 之间")
        if workflow_type == "expression":
            stored = self.expression_repository.load_project(project_id)
        elif workflow_type == "shrna_knockdown":
            stored = self.shrna_repository.load_project(project_id)
        elif workflow_type == "promoter_luciferase_reporter":
            stored = self.reporter_repository.load_project(project_id)
        else:
            raise ValueError(f"该操作不支持工作流：{workflow_type}")
        snapshot = stored.snapshot
        if snapshot.status != "rework_in_progress" or not snapshot.rework_owner_ids:
            raise ValueError("项目当前不在重做中，或没有记录重做对象")
        affected_labels: list[str] = []
        sample_names: list[str] = []
        if workflow_type == "shrna_knockdown":
            suffix_no = snapshot.experiment_attempt_no - 1
            rework_suffix = "n" if suffix_no == 1 else f"n{suffix_no}"
            for target in stored.design.targets:
                if target.target_id not in snapshot.rework_owner_ids:
                    continue
                affected_labels.append(f"shRNA-{target.target_no}")
                sample_names.extend(
                    f"{stored.design.gene_symbol}_{target.target_no}{rework_suffix}_{clone_no}"
                    for clone_no in range(1, clones_per_owner + 1)
                )
        else:
            suffix_no = snapshot.experiment_attempt_no - 1
            rework_suffix = "n" if suffix_no == 1 else f"n{suffix_no}"
            for construct in stored.design.constructs:
                if construct.construct_id not in snapshot.rework_owner_ids:
                    continue
                affected_labels.append(construct.construct_name)
                sample_names.extend(
                    f"{construct.construct_name}-{rework_suffix}-{clone_no}"
                    for clone_no in range(1, clones_per_owner + 1)
                )
        return AddOnSequencingPreview(
            project_id=project_id,
            workflow_type=workflow_type,
            round_no=len(snapshot.sequencing_submissions) + 1,
            affected_owner_ids=snapshot.rework_owner_ids,
            affected_owner_labels=tuple(affected_labels),
            clones_per_owner=clones_per_owner,
            sample_names=tuple(sample_names),
            reason="重新连接/转化后的新一轮送测",
        )

    def confirm_rework_submission(
        self,
        preview: AddOnSequencingPreview,
        *,
        actor: str,
        occurred_at: datetime,
        sequencing_vendor_name: str = "标准",
        sequencing_template_id: str | None = None,
        internal_submission_no: str = "",
        vendor_order_no: str = "",
    ) -> StoredExpressionProject | StoredShRNAProject | StoredReporterProject:
        fresh = self.preview_rework_submission(
            preview.project_id,
            workflow_type=preview.workflow_type,
            clones_per_owner=preview.clones_per_owner,
        )
        if fresh != preview:
            raise ValueError("项目数据已变化，请重新预览重做送测方案")
        if preview.workflow_type == "expression":
            repository = self.expression_repository
            stored = repository.load_project(preview.project_id)
            event_type = ExpressionAuditEvent
            sequencing_method, primer_name = "Nanopore", ""
        elif preview.workflow_type == "shrna_knockdown":
            repository = self.shrna_repository
            stored = repository.load_project(preview.project_id)
            event_type = ShRNAAuditEvent
            sequencing_method, primer_name = "Sanger", "U6"
        elif preview.workflow_type == "promoter_luciferase_reporter":
            repository = self.reporter_repository
            stored = repository.load_project(preview.project_id)
            event_type = ReporterAuditEvent
            sequencing_method, primer_name = "Nanopore", ""
        else:
            raise ValueError(f"该操作不支持工作流：{preview.workflow_type}")
        safe_vendor = re.sub(r'[<>:"/\\|?*]+', "_", sequencing_vendor_name).strip() or "标准"
        safe_gene = re.sub(r'[<>:"/\\|?*]+', "_", stored.design.gene_symbol).strip() or "gene"
        output_path = (
            stored.project_folder
            / "02_orders"
            / f"{occurred_at:%Y%m%d}-{safe_vendor}测序表-{safe_gene}_加测{preview.round_no}.xlsx"
        )
        if output_path.exists():
            output_path = output_path.with_name(
                f"{output_path.stem}-{occurred_at:%H%M%S}{output_path.suffix}",
            )
        records = tuple(
            {
                "sample_name": sample_name,
                "primer_name": primer_name,
                "gene_symbol": stored.design.gene_symbol,
                "clone_no": int(re.search(r"(\d+)$", sample_name).group(1)),
                "method": sequencing_method,
                "note": "重做后送测",
            }
            for sample_name in preview.sample_names
        )
        export_sequencing_form(
            records,
            output_path,
            template_store=self.workbook_template_store,
            contact_profile=self.contact_profile_store.load(),
            template_id=sequencing_template_id,
        )
        submission = SequencingSubmissionRecord(
            submission_id=f"submission-{uuid4()}",
            round_no=preview.round_no,
            submission_kind="post_rework",
            created_at=occurred_at,
            sent_at=occurred_at,
            sample_names=preview.sample_names,
            internal_submission_no=internal_submission_no.strip(),
            vendor_order_no=vendor_order_no.strip(),
            template_id=sequencing_template_id,
            form_path=str(output_path),
            experiment_attempt_no=stored.snapshot.experiment_attempt_no,
            status="sent",
            note=preview.reason,
        )
        snapshot = replace(
            stored.snapshot,
            status="sequencing_in_progress",
            sequencing_submissions=stored.snapshot.sequencing_submissions + (submission,),
        )
        snapshot = snapshot.append_status_event(
            event_type(
                event_id=f"confirm-rework-submission-{uuid4()}",
                event_type="post_rework_submission",
                occurred_at=occurred_at,
                actor=actor,
                from_status=stored.snapshot.status,
                to_status="sequencing_in_progress",
                note=f"{preview.reason}；{', '.join(preview.affected_owner_labels)}",
                source="user",
            ),
        )
        repository.save_snapshot(
            preview.project_id,
            snapshot,
            expected_revision=stored.snapshot.revision,
            updated_at=occurred_at,
        )
        repository.append_artifacts(
            preview.project_id,
            (
                GeneratedArtifact(
                    artifact_type="post_rework_sequencing_order_xlsx",
                    design_version_id=stored.design.design_version_id,
                    generated_at=occurred_at,
                    path=output_path,
                    content_sha256=_file_sha256(output_path),
                ),
            ),
        )
        return repository.load_project(preview.project_id)

    def mark_molecular_interrupted(
        self,
        project_id: str,
        *,
        workflow_type: str,
        interruption_type: str,
        note: str,
        actor: str,
        occurred_at: datetime,
    ) -> StoredExpressionProject | StoredShRNAProject | StoredReporterProject:
        if interruption_type not in {"pause", "abnormal"}:
            raise ValueError("中断类型必须选择 pause 或 abnormal")
        if not note.strip():
            raise ValueError("标记暂停/异常时必须填写原因")
        if workflow_type == "expression":
            repository = self.expression_repository
            stored = repository.load_project(project_id)
            event_type = ExpressionAuditEvent
        elif workflow_type == "shrna_knockdown":
            repository = self.shrna_repository
            stored = repository.load_project(project_id)
            event_type = ShRNAAuditEvent
        elif workflow_type == "promoter_luciferase_reporter":
            repository = self.reporter_repository
            stored = repository.load_project(project_id)
            event_type = ReporterAuditEvent
        else:
            raise ValueError(f"该操作不支持工作流：{workflow_type}")
        if stored.snapshot.status in {"abnormal_or_paused", "project_completed"}:
            raise ValueError("当前项目状态不能再次标记暂停/异常")
        calendar = ChinaBusinessCalendar.for_2026()
        effective_due = _effective_due_date(stored.snapshot, stored.due_date)
        frozen = calendar.remaining_workdays(occurred_at.date(), effective_due)
        snapshot = replace(
            stored.snapshot,
            status="abnormal_or_paused",
            effective_due_date=effective_due,
            interruption_type=interruption_type,
            interrupted_at=occurred_at,
            interrupted_previous_status=stored.snapshot.status,
            frozen_remaining_workdays=frozen,
        )
        snapshot = snapshot.append_status_event(
            event_type(
                event_id=f"mark-interrupted-{uuid4()}",
                event_type="mark_abnormal_or_paused",
                occurred_at=occurred_at,
                actor=actor,
                from_status=stored.snapshot.status,
                to_status="abnormal_or_paused",
                note=f"{interruption_type}：{note.strip()}",
                source="user",
            ),
        )
        repository.save_snapshot(
            project_id,
            snapshot,
            expected_revision=stored.snapshot.revision,
            updated_at=occurred_at,
        )
        return repository.load_project(project_id)

    def resume_molecular_project(
        self,
        project_id: str,
        *,
        workflow_type: str,
        actor: str,
        occurred_at: datetime,
        note: str,
        resume_status: str | None = None,
    ) -> StoredExpressionProject | StoredShRNAProject | StoredReporterProject:
        if not note.strip():
            raise ValueError("恢复项目时必须填写原因")
        if workflow_type == "expression":
            repository = self.expression_repository
            stored = repository.load_project(project_id)
            event_type = ExpressionAuditEvent
        elif workflow_type == "shrna_knockdown":
            repository = self.shrna_repository
            stored = repository.load_project(project_id)
            event_type = ShRNAAuditEvent
        elif workflow_type == "promoter_luciferase_reporter":
            repository = self.reporter_repository
            stored = repository.load_project(project_id)
            event_type = ReporterAuditEvent
        else:
            raise ValueError(f"该操作不支持工作流：{workflow_type}")
        snapshot = stored.snapshot
        if (
            snapshot.status != "abnormal_or_paused"
            or snapshot.interrupted_at is None
            or snapshot.interrupted_previous_status is None
            or snapshot.interruption_type is None
            or snapshot.frozen_remaining_workdays is None
        ):
            raise ValueError("项目当前没有可恢复的暂停/异常记录")
        target_status = resume_status or snapshot.interrupted_previous_status
        if target_status in {"abnormal_or_paused", "project_completed"}:
            raise ValueError("恢复目标状态无效")
        calendar = ChinaBusinessCalendar.for_2026()
        paused_workdays = calendar.workdays_in_half_open_interval(
            snapshot.interrupted_at.date(),
            occurred_at.date(),
        )
        effective_due = _effective_due_date(snapshot, stored.due_date)
        extended_due = calendar.add_workdays(effective_due, paused_workdays)
        interruption = ProjectInterruptionRecord(
            interruption_id=f"interruption-{uuid4()}",
            interruption_type=snapshot.interruption_type,
            started_at=snapshot.interrupted_at,
            resumed_at=occurred_at,
            previous_status=snapshot.interrupted_previous_status,
            resume_status=target_status,
            frozen_remaining_workdays=snapshot.frozen_remaining_workdays,
            paused_workdays=paused_workdays,
            start_note=next(
                (
                    event.note or ""
                    for event in reversed(snapshot.status_history)
                    if event.event_type == "mark_abnormal_or_paused"
                ),
                "",
            ),
            resume_note=note.strip(),
        )
        updated = replace(
            snapshot,
            status=target_status,
            effective_due_date=extended_due,
            interruption_type=None,
            interrupted_at=None,
            interrupted_previous_status=None,
            frozen_remaining_workdays=None,
            accumulated_paused_workdays=(
                snapshot.accumulated_paused_workdays + paused_workdays
            ),
            interruption_history=snapshot.interruption_history + (interruption,),
        )
        updated = updated.append_status_event(
            event_type(
                event_id=f"resume-interrupted-{uuid4()}",
                event_type="resume_from_abnormal_or_paused",
                occurred_at=occurred_at,
                actor=actor,
                from_status="abnormal_or_paused",
                to_status=target_status,
                note=(
                    f"{note.strip()}；暂停 {paused_workdays} 个工作日；"
                    f"完工日顺延至 {extended_due.isoformat()}"
                ),
                source="user",
            ),
        )
        repository.save_snapshot(
            project_id,
            updated,
            expected_revision=snapshot.revision,
            updated_at=occurred_at,
        )
        return repository.load_project(project_id)

    def _assert_unique_project_id(self, project_id: str) -> None:
        for repository in (
            self.repository,
            self.shrna_repository,
            self.expression_repository,
            self.reporter_repository,
        ):
            try:
                repository.load_project(project_id)
            except KeyError:
                continue
            raise DuplicateProjectError(f"项目号已存在：{project_id}")

    def create_expression_project(
        self,
        command: NewExpressionProjectCommand,
        *,
        created_at: datetime,
    ) -> StoredExpressionProject:
        _validate_project_dates(command.received_date, command.due_date)
        self._assert_unique_project_id(command.project_id)
        design = create_expression_design(
            ExpressionDesignInput(
                project_id=command.project_id,
                gene_symbol=command.gene_symbol,
                species=command.species,
                source_cds=command.source_cds,
                construct_lines=command.construct_lines,
                transcript_accession=command.transcript_accession,
            ),
            expression_rules_from_protocol(command.protocol),
            design_version_id=f"{command.project_id}-v1",
            created_at=created_at,
        )
        design = replace(design, gene_id=command.gene_id)
        if design.requires_confirmation:
            if not (command.design_confirmation_reason or "").strip():
                raise DesignConfirmationRequired(
                    "表达设计含点突变或需复核项目，请确认最终序列并填写原因",
                )
            design = confirm_expression_design(
                design,
                confirmation_id=f"expression-confirm-{uuid4()}",
                reason=command.design_confirmation_reason or "",
                actor=command.actor,
                occurred_at=created_at,
            )
        vector_design = apply_expression_protocol(
            design,
            command.vector,
            command.protocol,
        )
        snapshot = ExpressionProjectSnapshot(
            project_id=command.project_id,
            revision=1,
            status="design_completed",
            active_design_version_id=design.design_version_id,
            clones_per_construct=command.clones_per_construct,
            status_history=(
                ExpressionAuditEvent(
                    event_id=f"complete-design-{uuid4()}",
                    event_type="complete_design",
                    occurred_at=created_at,
                    actor=command.actor,
                    from_status="recorded",
                    to_status="design_completed",
                    source="application",
                ),
            ),
        )
        workspace: ProjectWorkspace | None = None
        database_created = False
        try:
            workspace = create_project_folder(
                self.projects_root,
                project_id=command.project_id,
                target_name=command.gene_symbol,
                folder_suffix="OE",
            )
            calendar = ChinaBusinessCalendar.for_2026()
            bundle = export_expression_bundle(
                design,
                vector_design,
                workspace,
                generated_at=created_at,
                primer_order_date=created_at.date(),
                sequencing_order_date=calendar.add_workdays(created_at.date(), 2),
                clones_per_construct=command.clones_per_construct,
                primer_vendor_name=command.primer_vendor_name,
                sequencing_vendor_name=command.sequencing_vendor_name,
                sequencing_method=command.sequencing_method,
                workbook_template_store=self.workbook_template_store,
                contact_profile=self.contact_profile_store.load(),
                primer_template_id=command.primer_template_id,
                sequencing_template_id=command.sequencing_template_id,
            )
            self.expression_repository.create_project(
                project_id=command.project_id,
                gene_symbol=command.gene_symbol,
                species=command.species,
                vector_name=command.protocol.display_name,
                received_date=command.received_date,
                due_date=command.due_date,
                project_folder=workspace.root,
                design=design,
                vector_design=vector_design,
                snapshot=snapshot,
                created_at=created_at,
            )
            database_created = True
            self.expression_repository.append_artifacts(
                command.project_id,
                bundle.artifacts,
            )
        except Exception:
            if workspace is not None and not database_created:
                self._remove_new_workspace(workspace)
            raise
        return self.expression_repository.load_project(command.project_id)

    def create_reporter_project(
        self,
        command: NewReporterProjectCommand,
        *,
        created_at: datetime,
    ) -> StoredReporterProject:
        _validate_project_dates(command.received_date, command.due_date)
        self._assert_unique_project_id(command.project_id)
        design = create_reporter_design(
            ReporterDesignInput(
                project_id=command.project_id,
                gene_symbol=command.gene_symbol,
                species=command.species,
                promoter_sequence=command.promoter_sequence,
                construct_lines=command.construct_lines,
                mutation_definitions=command.mutation_definitions,
                transcript_accession=command.transcript_accession,
            ),
            protocol_version_id=command.protocol.protocol_version_id,
            design_version_id=f"{command.project_id}-v1",
            created_at=created_at,
        )
        design = replace(design, gene_id=command.gene_id)
        if design.requires_confirmation:
            if not (command.design_confirmation_reason or "").strip():
                raise DesignConfirmationRequired(
                    "reporter 突变设计需要确认最终替换序列并填写原因",
                )
            design = confirm_reporter_design(
                design,
                confirmation_id=f"reporter-confirm-{uuid4()}",
                reason=command.design_confirmation_reason or "",
                actor=command.actor,
                occurred_at=created_at,
            )
        vector_design = apply_reporter_protocol(
            design,
            command.vector,
            command.protocol,
        )
        snapshot = ReporterProjectSnapshot(
            project_id=command.project_id,
            revision=1,
            status="design_completed",
            active_design_version_id=design.design_version_id,
            clones_per_construct=command.clones_per_construct,
            status_history=(
                ReporterAuditEvent(
                    event_id=f"complete-design-{uuid4()}",
                    event_type="complete_design",
                    occurred_at=created_at,
                    actor=command.actor,
                    from_status="recorded",
                    to_status="design_completed",
                    source="application",
                ),
            ),
        )
        workspace: ProjectWorkspace | None = None
        database_created = False
        try:
            workspace = create_project_folder(
                self.projects_root,
                project_id=command.project_id,
                target_name=command.gene_symbol,
                folder_suffix="RPT",
            )
            calendar = ChinaBusinessCalendar.for_2026()
            bundle = export_reporter_bundle(
                design,
                vector_design,
                workspace,
                generated_at=created_at,
                primer_order_date=created_at.date(),
                sequencing_order_date=calendar.add_workdays(created_at.date(), 2),
                clones_per_construct=command.clones_per_construct,
                primer_vendor_name=command.primer_vendor_name,
                sequencing_vendor_name=command.sequencing_vendor_name,
                sequencing_method=command.sequencing_method,
                workbook_template_store=self.workbook_template_store,
                contact_profile=self.contact_profile_store.load(),
                primer_template_id=command.primer_template_id,
                sequencing_template_id=command.sequencing_template_id,
            )
            self.reporter_repository.create_project(
                project_id=command.project_id,
                gene_symbol=command.gene_symbol,
                species=command.species,
                vector_name=command.protocol.display_name,
                received_date=command.received_date,
                due_date=command.due_date,
                project_folder=workspace.root,
                design=design,
                vector_design=vector_design,
                snapshot=snapshot,
                created_at=created_at,
            )
            database_created = True
            self.reporter_repository.append_artifacts(
                command.project_id,
                bundle.artifacts,
            )
        except Exception:
            if workspace is not None and not database_created:
                self._remove_new_workspace(workspace)
            raise
        return self.reporter_repository.load_project(command.project_id)

    def analyze_reporter_sequencing(
        self,
        project_id: str,
        *,
        actor: str,
        analyzed_at: datetime,
    ) -> ReporterAnalysisOutcome:
        stored = self.reporter_repository.load_project(project_id)
        expected_by_clone = {
            f"{construct.construct_name}-{clone_no}": construct
            for construct in stored.design.constructs
            for clone_no in range(1, stored.snapshot.clones_per_construct + 1)
        }
        for submission in stored.snapshot.sequencing_submissions:
            for sample_name in submission.sample_names:
                construct = next(
                    (
                        item
                        for item in stored.design.constructs
                        if sample_name.startswith(f"{item.construct_name}-")
                    ),
                    None,
                )
                if construct is not None:
                    expected_by_clone[sample_name] = construct
        match_plan = match_shrna_sequence_files(
            stored.project_folder / "03_sequencing",
            tuple(expected_by_clone),
        )
        records = []
        for clone_name, construct in expected_by_clone.items():
            paths = match_plan.files_for(clone_name)
            if not paths:
                records.append(
                    ReporterCloneResultRecord(
                        result_id=f"reporter-seq-{uuid4()}",
                        clone_name=clone_name,
                        construct_id=construct.construct_id,
                        status="warning",
                        reason="未找到测序文件",
                        analyzed_at=analyzed_at,
                        source_files=(),
                    ),
                )
                continue
            judgments = []
            errors = []
            for path in paths:
                try:
                    read = read_sequence_file(path)
                    judgments.append(
                        judge_expression_read(
                            clone_name=clone_name,
                            construct_id=construct.construct_id,
                            read_sequence=read.sequence,
                            expected_insert_sequence=construct.insert_sequence,
                            expected_coding_sequence=construct.insert_sequence,
                        ),
                    )
                except Exception as error:
                    errors.append(f"{path.name}: {error}")
            passed = next(
                (
                    item
                    for item in judgments
                    if item.status is ExpressionCloneJudgmentStatus.PASS
                ),
                None,
            )
            if passed is not None:
                selected = passed
                status = "pass"
                reason = passed.reason
            elif judgments:
                selected = max(judgments, key=lambda item: (item.coverage, item.identity))
                status = "warning"
                reason = selected.reason
                if errors:
                    reason += f"；另有 {len(errors)} 个文件无法读取"
            else:
                selected = None
                status = "warning"
                reason = "；".join(errors) or "没有可判读的测序文件"
            records.append(
                ReporterCloneResultRecord(
                    result_id=f"reporter-seq-{uuid4()}",
                    clone_name=clone_name,
                    construct_id=construct.construct_id,
                    status=status,
                    reason=reason,
                    analyzed_at=analyzed_at,
                    source_files=tuple(str(path) for path in paths),
                    coverage=selected.coverage if selected is not None else 0.0,
                    identity=selected.identity if selected is not None else 0.0,
                    substitution_count=(
                        selected.substitution_count if selected is not None else 0
                    ),
                    insertion_count=selected.insertion_count if selected is not None else 0,
                    deletion_count=selected.deletion_count if selected is not None else 0,
                    frameshift=selected.frameshift if selected is not None else False,
                ),
            )
        updated = stored.snapshot.append_clone_results(tuple(records))
        updated = updated.append_status_event(
            ReporterAuditEvent(
                event_id=f"analyze-reporter-sequencing-{uuid4()}",
                event_type="analyze_sequencing",
                occurred_at=analyzed_at,
                actor=actor,
                from_status=stored.snapshot.status,
                to_status="analysis_completed",
                source="application",
            ),
        )
        updated = _mark_latest_submission_analyzed(updated)
        updated = replace(updated, status="analysis_completed")
        artifact = export_reporter_analysis_report(
            stored.design,
            stored.vector_design,
            tuple(records),
            stored.project_folder / "04_reports",
            analyzed_at=analyzed_at,
            unmatched_files=match_plan.unmatched_files,
            ambiguous_files=match_plan.ambiguous_files,
        )
        self.reporter_repository.save_snapshot(
            project_id,
            updated,
            expected_revision=stored.snapshot.revision,
            updated_at=analyzed_at,
        )
        self.reporter_repository.append_artifacts(project_id, (artifact,))
        return ReporterAnalysisOutcome(
            project=self.reporter_repository.load_project(project_id),
            analysis_report=artifact.path,
            unmatched_files=match_plan.unmatched_files,
            ambiguous_files=match_plan.ambiguous_files,
        )

    def analyze_expression_sequencing(
        self,
        project_id: str,
        *,
        actor: str,
        analyzed_at: datetime,
    ) -> ExpressionAnalysisOutcome:
        stored = self.expression_repository.load_project(project_id)
        expected_by_clone = {}
        for construct in stored.design.constructs:
            for clone_no in range(1, stored.snapshot.clones_per_construct + 1):
                expected_by_clone[f"{construct.construct_name}-{clone_no}"] = construct
        for submission in stored.snapshot.sequencing_submissions:
            for sample_name in submission.sample_names:
                construct = next(
                    (
                        item
                        for item in stored.design.constructs
                        if sample_name.startswith(f"{item.construct_name}-")
                    ),
                    None,
                )
                if construct is not None:
                    expected_by_clone[sample_name] = construct
        match_plan = match_shrna_sequence_files(
            stored.project_folder / "03_sequencing",
            tuple(expected_by_clone),
        )
        records: list[ExpressionCloneResultRecord] = []
        for clone_name, construct in expected_by_clone.items():
            paths = match_plan.files_for(clone_name)
            if not paths:
                records.append(
                    ExpressionCloneResultRecord(
                        result_id=f"expression-seq-{uuid4()}",
                        clone_name=clone_name,
                        construct_id=construct.construct_id,
                        status="warning",
                        reason="未找到测序文件",
                        analyzed_at=analyzed_at,
                        source_files=(),
                    ),
                )
                continue
            judgments = []
            errors = []
            for path in paths:
                try:
                    read = read_sequence_file(path)
                    judgments.append(
                        judge_expression_read(
                            clone_name=clone_name,
                            construct_id=construct.construct_id,
                            read_sequence=read.sequence,
                            expected_insert_sequence=construct.insert_sequence,
                            expected_coding_sequence=construct.coding_sequence,
                        ),
                    )
                except Exception as error:
                    errors.append(f"{path.name}: {error}")
            passed = next(
                (
                    item
                    for item in judgments
                    if item.status is ExpressionCloneJudgmentStatus.PASS
                ),
                None,
            )
            if passed is not None:
                selected = passed
                status = "pass"
                reason = (
                    f"{len(paths)} 个文件中至少 1 个完整匹配：{selected.reason}"
                    if len(paths) > 1
                    else selected.reason
                )
            elif judgments:
                selected = max(
                    judgments,
                    key=lambda item: (item.coverage, item.identity),
                )
                status = "warning"
                reason = selected.reason
                if errors:
                    reason += f"；另有 {len(errors)} 个文件无法读取"
            else:
                selected = None
                status = "warning"
                reason = "；".join(errors) or "没有可判读的测序文件"
            records.append(
                ExpressionCloneResultRecord(
                    result_id=f"expression-seq-{uuid4()}",
                    clone_name=clone_name,
                    construct_id=construct.construct_id,
                    status=status,
                    reason=reason,
                    analyzed_at=analyzed_at,
                    source_files=tuple(str(path) for path in paths),
                    coverage=selected.coverage if selected is not None else 0.0,
                    identity=selected.identity if selected is not None else 0.0,
                    substitution_count=(
                        selected.substitution_count if selected is not None else 0
                    ),
                    insertion_count=(
                        selected.insertion_count if selected is not None else 0
                    ),
                    deletion_count=(
                        selected.deletion_count if selected is not None else 0
                    ),
                    frameshift=selected.frameshift if selected is not None else False,
                    premature_stop=(
                        selected.premature_stop if selected is not None else False
                    ),
                ),
            )
        updated = stored.snapshot.append_clone_results(tuple(records))
        updated = updated.append_status_event(
            ExpressionAuditEvent(
                event_id=f"analyze-expression-sequencing-{uuid4()}",
                event_type="analyze_sequencing",
                occurred_at=analyzed_at,
                actor=actor,
                from_status=stored.snapshot.status,
                to_status="analysis_completed",
                source="application",
            ),
        )
        updated = _mark_latest_submission_analyzed(updated)
        updated = replace(updated, status="analysis_completed")
        artifact = export_expression_analysis_report(
            stored.design,
            stored.vector_design,
            tuple(records),
            stored.project_folder / "04_reports",
            analyzed_at=analyzed_at,
            unmatched_files=match_plan.unmatched_files,
            ambiguous_files=match_plan.ambiguous_files,
        )
        self.expression_repository.save_snapshot(
            project_id,
            updated,
            expected_revision=stored.snapshot.revision,
            updated_at=analyzed_at,
        )
        self.expression_repository.append_artifacts(project_id, (artifact,))
        return ExpressionAnalysisOutcome(
            project=self.expression_repository.load_project(project_id),
            analysis_report=artifact.path,
            unmatched_files=match_plan.unmatched_files,
            ambiguous_files=match_plan.ambiguous_files,
        )

    def confirm_expression_clone_review(
        self,
        project_id: str,
        *,
        clone_name: str,
        usable: bool,
        note: str,
        actor: str,
        occurred_at: datetime,
    ) -> StoredExpressionProject:
        clean_note = note.strip()
        if not clean_note:
            raise ValueError("人工复核必须填写说明")
        stored = self.expression_repository.load_project(project_id)
        latest = None
        for record in stored.snapshot.clone_results:
            if record.clone_name == clone_name:
                latest = record
        if latest is None:
            raise KeyError(f"未找到克隆判读记录：{clone_name}")
        if latest.status != "warning":
            raise ValueError("只有 WARNING 克隆需要人工复核")
        reviewed = replace(
            latest,
            result_id=f"expression-review-{uuid4()}",
            analyzed_at=occurred_at,
            manually_confirmed_usable=usable,
            manual_review_status="usable" if usable else "unusable",
            manual_note=clean_note,
            reason=f"{latest.reason}；人工复核：{clean_note}",
        )
        updated = stored.snapshot.append_clone_results((reviewed,))
        updated = updated.append_status_event(
            ExpressionAuditEvent(
                event_id=f"review-expression-clone-{uuid4()}",
                event_type=(
                    "confirm_clone_usable" if usable else "confirm_clone_unusable"
                ),
                occurred_at=occurred_at,
                actor=actor,
                from_status=stored.snapshot.status,
                to_status=stored.snapshot.status,
                note=f"{clone_name}：{clean_note}",
                source="user",
            ),
        )
        self.expression_repository.save_snapshot(
            project_id,
            updated,
            expected_revision=stored.snapshot.revision,
            updated_at=occurred_at,
        )
        return self.expression_repository.load_project(project_id)

    def confirm_shrna_clone_review(
        self,
        project_id: str,
        *,
        clone_name: str,
        usable: bool,
        note: str,
        actor: str,
        reviewed_at: datetime,
    ) -> StoredShRNAProject:
        clean_note = note.strip()
        if not clean_note:
            raise ValueError("人工复核必须填写说明")
        stored = self.shrna_repository.load_project(project_id)
        latest = None
        for record in stored.snapshot.clone_results:
            if record.clone_name == clone_name:
                latest = record
        if latest is None:
            raise KeyError(f"未找到克隆判读记录：{clone_name}")
        if latest.status != "warning":
            raise ValueError("只有 WARNING 克隆需要人工复核")
        reviewed = replace(
            latest,
            result_id=f"shrna-review-{uuid4()}",
            analyzed_at=reviewed_at,
            manually_confirmed_usable=usable,
            manual_review_status="usable" if usable else "unusable",
            manual_note=clean_note,
            reason=f"{latest.reason}；人工复核：{clean_note}",
        )
        updated = stored.snapshot.append_clone_results((reviewed,))
        updated = updated.append_status_event(
            ShRNAAuditEvent(
                event_id=f"review-shrna-clone-{uuid4()}",
                event_type="confirm_clone_usable" if usable else "confirm_clone_unusable",
                occurred_at=reviewed_at,
                actor=actor,
                from_status=stored.snapshot.status,
                to_status=stored.snapshot.status,
                note=f"{clone_name}：{clean_note}",
                source="user",
            ),
        )
        self.shrna_repository.save_snapshot(
            project_id,
            updated,
            expected_revision=stored.snapshot.revision,
            updated_at=reviewed_at,
        )
        return self.shrna_repository.load_project(project_id)

    def create_shrna_project(
        self,
        command: NewShRNAProjectCommand,
        *,
        created_at: datetime,
    ) -> StoredShRNAProject:
        _validate_project_dates(command.received_date, command.due_date)
        self._assert_unique_project_id(command.project_id)
        vector, protocol = load_public_plko1_puro_starter(
            user_confirmed=command.vector_sequence_confirmed,
        )
        selection = select_initial_candidates(
            command.candidates,
            target_count=command.target_count,
        )
        design_input = ShRNADesignInput(
            project_id=command.project_id,
            gene_symbol=command.gene_symbol,
            species=command.species,
            cds_sequence=normalize_dna(command.cds_sequence),
            vector_protocol_version_id=protocol.protocol_version_id,
            transcript_accession=command.transcript_accession,
            gene_id=command.gene_id,
            ccds_id=command.ccds_id,
            target_count=command.target_count,
            clones_per_target=command.clones_per_target,
        )
        design = create_shrna_design(
            design_input,
            selection.selected,
            vector,
            protocol,
            design_version_id=f"{command.project_id}-v1",
            created_at=created_at,
        )
        snapshot = ShRNAProjectSnapshot(
            project_id=command.project_id,
            revision=1,
            status="design_completed",
            active_design_version_id=design.design_version_id,
            clone_results=(),
            status_history=(
                ShRNAAuditEvent(
                    event_id=f"complete-design-{uuid4()}",
                    event_type="complete_design",
                    occurred_at=created_at,
                    actor=command.actor,
                    from_status="recorded",
                    to_status="design_completed",
                    source="application",
                ),
            ),
        )
        workspace: ProjectWorkspace | None = None
        database_created = False
        try:
            workspace = create_project_folder(
                self.projects_root,
                project_id=command.project_id,
                target_name=command.gene_symbol,
                folder_suffix="KD",
            )
            calendar = ChinaBusinessCalendar.for_2026()
            bundle = export_shrna_bundle(
                design,
                workspace,
                generated_at=created_at,
                primer_order_date=created_at.date(),
                sequencing_order_date=calendar.add_workdays(created_at.date(), 2),
                primer_vendor_name=command.primer_vendor_name,
                sequencing_vendor_name=command.sequencing_vendor_name,
                workbook_template_store=self.workbook_template_store,
                contact_profile=self.contact_profile_store.load(),
                primer_template_id=command.primer_template_id,
                sequencing_template_id=command.sequencing_template_id,
            )
            self.shrna_repository.create_project(
                project_id=command.project_id,
                gene_symbol=command.gene_symbol,
                species=command.species,
                received_date=command.received_date,
                due_date=command.due_date,
                project_folder=workspace.root,
                design=design,
                snapshot=snapshot,
                created_at=created_at,
            )
            database_created = True
            self.shrna_repository.append_artifacts(command.project_id, bundle.artifacts)
        except Exception:
            if workspace is not None and not database_created:
                self._remove_new_workspace(workspace)
            raise
        return self.shrna_repository.load_project(command.project_id)

    def analyze_shrna_sequencing(
        self,
        project_id: str,
        *,
        actor: str,
        analyzed_at: datetime,
    ) -> ShRNAAnalysisOutcome:
        stored = self.shrna_repository.load_project(project_id)
        target_by_clone = {
            clone_name: target
            for target in stored.design.targets
            for clone_name in target.clone_names
        }
        for submission in stored.snapshot.sequencing_submissions:
            for sample_name in submission.sample_names:
                target = next(
                    (
                        item
                        for item in stored.design.targets
                        if sample_name.startswith(f"{item.clone_names[0].rsplit('-', 1)[0]}-")
                        or re.fullmatch(
                            rf"{re.escape(stored.design.gene_symbol)}_"
                            rf"{item.target_no}n(?:\d+)?_\d+",
                            sample_name,
                        )
                    ),
                    None,
                )
                if target is not None:
                    target_by_clone[sample_name] = target
        clone_names = tuple(target_by_clone)
        match_plan = match_shrna_sequence_files(
            stored.project_folder / "03_sequencing",
            clone_names,
        )
        records: list[ShRNACloneResultRecord] = []
        for clone_name in clone_names:
            target = target_by_clone[clone_name]
            paths = match_plan.files_for(clone_name)
            if not paths:
                records.append(
                    ShRNACloneResultRecord(
                        result_id=f"seq-result-{uuid4()}",
                        clone_name=clone_name,
                        target_id=target.target_id,
                        status="warning",
                        reason="未找到测序文件",
                        analyzed_at=analyzed_at,
                        source_files=(),
                    ),
                )
                continue
            judgments = []
            read_errors = []
            for path in paths:
                try:
                    sequencing_read = read_sequence_file(path)
                    judgments.append(
                        judge_shrna_read(
                            clone_name=clone_name,
                            target_id=target.target_id,
                            read_sequence=sequencing_read.sequence,
                            expected_forward_oligo=target.oligos.forward_sequence,
                            source_files=(path,),
                        ),
                    )
                except Exception as error:
                    read_errors.append(f"{path.name}: {error}")
            statuses = {item.status.value for item in judgments}
            if read_errors or not judgments:
                status = CloneJudgmentStatus.WARNING.value
                reason = "；".join(read_errors) or "没有可判读的测序文件"
                match_start = None
            elif len(statuses) > 1:
                status = CloneJudgmentStatus.WARNING.value
                reason = "同一克隆的多个测序文件判读结论不一致，需要人工确认"
                match_start = None
            else:
                status = judgments[0].status.value
                reason = judgments[0].reason
                if len(judgments) > 1:
                    reason = f"{len(judgments)} 个测序文件判读一致：{reason}"
                match_start = next(
                    (
                        item.match_start
                        for item in judgments
                        if item.match_start is not None
                    ),
                    None,
                )
            records.append(
                ShRNACloneResultRecord(
                    result_id=f"seq-result-{uuid4()}",
                    clone_name=clone_name,
                    target_id=target.target_id,
                    status=status,
                    reason=reason,
                    analyzed_at=analyzed_at,
                    source_files=tuple(str(path) for path in paths),
                    match_start=match_start,
                ),
            )

        updated = stored.snapshot.append_clone_results(tuple(records))
        updated = updated.append_status_event(
            ShRNAAuditEvent(
                event_id=f"analyze-sequencing-{uuid4()}",
                event_type="analyze_sequencing",
                occurred_at=analyzed_at,
                actor=actor,
                from_status=stored.snapshot.status,
                to_status="analysis_completed",
                source="application",
            ),
        )
        updated = _mark_latest_submission_analyzed(updated)
        updated = replace(updated, status="analysis_completed")
        artifact = export_shrna_analysis_report(
            stored.design,
            tuple(records),
            stored.project_folder / "04_reports",
            analyzed_at=analyzed_at,
            unmatched_files=match_plan.unmatched_files,
            ambiguous_files=match_plan.ambiguous_files,
        )
        self.shrna_repository.save_snapshot(
            project_id,
            updated,
            expected_revision=stored.snapshot.revision,
            updated_at=analyzed_at,
        )
        self.shrna_repository.append_artifacts(project_id, (artifact,))
        return ShRNAAnalysisOutcome(
            project=self.shrna_repository.load_project(project_id),
            analysis_report=artifact.path,
            unmatched_files=match_plan.unmatched_files,
            ambiguous_files=match_plan.ambiguous_files,
        )

    def prepare_syn_project(
        self,
        command: NewSYNProjectCommand,
        *,
        created_at: datetime,
    ) -> PreparedSYNProject:
        _validate_project_dates(command.received_date, command.due_date)
        try:
            self.repository.load_project(command.project_id)
        except KeyError:
            pass
        else:
            raise DuplicateProjectError(f"项目号已存在：{command.project_id}")
        vector, protocol = load_public_puc57_starter(
            command.linearization_site,
            user_confirmed=command.vector_sequence_confirmed,
        )
        design_input = SYNDesignInput(
            project_id=command.project_id,
            target_name=command.target_name,
            raw_sequence=command.raw_sequence,
            input_format=command.input_format,
            vector_protocol_version_id=protocol.protocol_version_id,
        )
        design = create_syn_design(
            design_input,
            vector,
            protocol,
            design_version_id=f"{command.project_id}-v1",
            created_at=created_at,
        )
        return PreparedSYNProject(
            vector=vector,
            protocol=protocol,
            design=design,
        )

    def save_prepared_syn_project(
        self,
        command: NewSYNProjectCommand,
        prepared: PreparedSYNProject,
        *,
        design_confirmation_reason: str | None,
        created_at: datetime,
    ) -> StoredSYNProject:
        _validate_project_dates(command.received_date, command.due_date)
        design = prepared.design
        if design.project_id != command.project_id:
            raise ValueError("预计算设计与待保存项目号不一致")
        if design.requires_confirmation:
            if not (design_confirmation_reason or "").strip():
                raise DesignConfirmationRequired(
                    "该设计含需人工复核项，请填写确认原因后保存",
                )
            design = confirm_syn_design_warnings(
                design,
                override_id=f"design-confirm-{uuid4()}",
                reason=design_confirmation_reason,
                actor=command.actor,
                occurred_at=created_at,
            )
        vector_confirmation = SYNManualOverrideRecord(
            override_id=f"vector-confirm-{uuid4()}",
            field_path="vector_reference_confirmation",
            old_value="unconfirmed",
            new_value="confirmed",
            reason="用户确认实际使用载体与内置 SnapGene pUC57 公开参考序列一致",
            occurred_at=created_at,
            actor=command.actor,
        )
        design = replace(
            design,
            manual_overrides=design.manual_overrides + (vector_confirmation,),
        )
        snapshot = SYNProjectSnapshot(
            project_id=command.project_id,
            revision=1,
            status="design_completed",
            resuspension_data_status="missing",
            syn_assembly_round_no=0,
            syn_assembly_substep=None,
            active_design_version_id=design.design_version_id,
            attempts=(),
            colonies=(),
            prep_records=(),
            sequencing_confirmations=(),
            status_history=(
                SYNAuditEvent(
                    event_id=f"complete-design-{uuid4()}",
                    event_type="complete_design",
                    occurred_at=created_at,
                    actor=command.actor,
                    from_status="recorded",
                    to_status="design_completed",
                    source="application",
                ),
            ),
            manual_override_history=(),
        )
        workspace: ProjectWorkspace | None = None
        database_created = False
        try:
            workspace = create_project_folder(
                self.projects_root,
                project_id=command.project_id,
                target_name=command.target_name,
            )
            bundle = export_syn_bundle(
                design,
                workspace.folder("01_design"),
                target_name=command.target_name,
                generated_at=created_at,
            )
            self.repository.create_project(
                project_id=command.project_id,
                target_name=command.target_name,
                received_date=command.received_date,
                due_date=command.due_date,
                project_folder=workspace.root,
                design=design,
                snapshot=snapshot,
                created_at=created_at,
            )
            database_created = True
            self.repository.append_artifacts(command.project_id, bundle.artifacts)
        except Exception:
            if workspace is not None and not database_created:
                self._remove_new_workspace(workspace)
            raise
        return self.repository.load_project(command.project_id)

    def _remove_new_workspace(self, workspace: ProjectWorkspace) -> None:
        projects_root = self.projects_root.resolve()
        workspace_root = workspace.root.resolve()
        if workspace_root.parent != projects_root:
            raise RuntimeError("拒绝清理项目根目录之外的路径")
        shutil.rmtree(workspace_root)
