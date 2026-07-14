"""Promoter-luciferase 报告载体的领域记录。"""

from __future__ import annotations

from dataclasses import dataclass, replace
from datetime import date, datetime

from .sequencing_submission import SequencingSubmissionRecord
from .project_tracking import ProjectInterruptionRecord


def _strict_dna(field_name: str, value: str) -> str:
    lines = value.lstrip("\ufeff").splitlines()
    normalized = "".join(
        line.strip() for line in lines if not line.lstrip().startswith(">")
    ).upper()
    if not normalized:
        raise ValueError(f"{field_name} 不能为空")
    invalid = sorted(set(normalized) - set("ACGT"))
    if invalid:
        raise ValueError(f"{field_name} 含不支持碱基：{', '.join(invalid)}")
    return normalized


@dataclass(frozen=True, slots=True)
class ReporterDesignInput:
    project_id: str
    gene_symbol: str
    species: str
    promoter_sequence: str
    construct_lines: tuple[str, ...]
    mutation_definitions: tuple[str, ...] = ()
    transcript_accession: str | None = None
    promoter_length_requested: int | None = None

    def __post_init__(self) -> None:
        if not self.project_id.strip() or not self.gene_symbol.strip():
            raise ValueError("project_id 和 gene_symbol 不能为空")
        species = self.species.strip().lower()
        if species not in {"human", "mouse", "rat"}:
            raise ValueError("species 必须是 human、mouse 或 rat")
        object.__setattr__(self, "species", species)
        object.__setattr__(self, "gene_symbol", self.gene_symbol.strip())
        object.__setattr__(
            self,
            "promoter_sequence",
            _strict_dna("promoter_sequence", self.promoter_sequence),
        )
        object.__setattr__(
            self,
            "construct_lines",
            tuple(line.strip() for line in self.construct_lines if line.strip()),
        )
        object.__setattr__(
            self,
            "mutation_definitions",
            tuple(line.strip() for line in self.mutation_definitions if line.strip()),
        )
        if not self.construct_lines:
            raise ValueError("construct_lines 不能为空")


@dataclass(frozen=True, slots=True)
class PromoterMutationDefinition:
    name: str
    start: int
    end: int
    replacement_sequence: str
    original_sequence: str


@dataclass(frozen=True, slots=True)
class ReporterConstruct:
    construct_id: str
    construct_name: str
    request_line: str
    retained_promoter_length: int
    retained_source_start: int
    mutation_names: tuple[str, ...]
    insert_sequence: str


@dataclass(frozen=True, slots=True)
class ReporterDesignVersion:
    design_version_id: str
    project_id: str
    gene_symbol: str
    species: str
    transcript_accession: str | None
    promoter_source_checksum: str
    protocol_version_id: str
    created_at: datetime
    mutation_definitions: tuple[PromoterMutationDefinition, ...]
    constructs: tuple[ReporterConstruct, ...]
    design_warnings: tuple[str, ...]
    requires_confirmation: bool
    confirmation_history: tuple[ReporterDesignConfirmation, ...] = ()
    gene_id: str | None = None


@dataclass(frozen=True, slots=True)
class ReporterDesignConfirmation:
    confirmation_id: str
    reason: str
    actor: str
    occurred_at: datetime


@dataclass(frozen=True, slots=True)
class ReporterPrimerPlan:
    primer_id: str
    name: str
    sequence: str
    direction: str
    anneal_length: int


@dataclass(frozen=True, slots=True)
class ReporterConstructVectorPlan:
    construct_id: str
    construct_name: str
    forward_primer: ReporterPrimerPlan
    reverse_primer: ReporterPrimerPlan
    expected_plasmid_sequence: str
    expected_plasmid_checksum: str


@dataclass(frozen=True, slots=True)
class ReporterVectorDesignResult:
    design_version_id: str
    vector_record_id: str
    vector_checksum: str
    protocol_version_id: str
    construct_plans: tuple[ReporterConstructVectorPlan, ...]


@dataclass(frozen=True, slots=True)
class ReporterAuditEvent:
    event_id: str
    event_type: str
    occurred_at: datetime
    actor: str
    from_status: str | None = None
    to_status: str | None = None
    note: str | None = None
    source: str = "user"


@dataclass(frozen=True, slots=True)
class ReporterCloneResultRecord:
    result_id: str
    clone_name: str
    construct_id: str
    status: str
    reason: str
    analyzed_at: datetime
    source_files: tuple[str, ...]
    coverage: float = 0.0
    identity: float = 0.0
    substitution_count: int = 0
    insertion_count: int = 0
    deletion_count: int = 0
    frameshift: bool = False
    manually_confirmed_usable: bool = False
    manual_review_status: str | None = None
    manual_note: str | None = None

    def __post_init__(self) -> None:
        if self.status not in {"pass", "fail", "warning"}:
            raise ValueError("clone result status must be pass, fail, or warning")


@dataclass(frozen=True, slots=True)
class ReporterProjectSnapshot:
    project_id: str
    revision: int
    status: str
    active_design_version_id: str
    status_history: tuple[ReporterAuditEvent, ...]
    clones_per_construct: int = 5
    clone_results: tuple[ReporterCloneResultRecord, ...] = ()
    selected_prep_clone_names: tuple[str, ...] = ()
    plasmid_prep_completed_at: datetime | None = None
    actual_completed_at: datetime | None = None
    sequencing_submissions: tuple[SequencingSubmissionRecord, ...] = ()
    experiment_attempt_no: int = 1
    rework_owner_ids: tuple[str, ...] = ()
    effective_due_date: date | None = None
    interruption_type: str | None = None
    interrupted_at: datetime | None = None
    interrupted_previous_status: str | None = None
    frozen_remaining_workdays: int | None = None
    accumulated_paused_workdays: int = 0
    interruption_history: tuple[ProjectInterruptionRecord, ...] = ()
    internal_project_no: str = ""
    primer_submission_no: str = ""
    primer_vendor_order_no: str = ""

    def __post_init__(self) -> None:
        if not 1 <= self.clones_per_construct <= 96:
            raise ValueError("clones_per_construct must be between 1 and 96")

    def append_status_event(self, event: ReporterAuditEvent) -> ReporterProjectSnapshot:
        if any(item.event_id == event.event_id for item in self.status_history):
            raise ValueError("status event_id must be unique")
        return replace(
            self,
            revision=self.revision + 1,
            status_history=self.status_history + (event,),
        )

    def append_clone_results(
        self,
        records: tuple[ReporterCloneResultRecord, ...],
    ) -> ReporterProjectSnapshot:
        known = {item.result_id for item in self.clone_results}
        incoming = [item.result_id for item in records]
        if len(incoming) != len(set(incoming)) or any(item in known for item in incoming):
            raise ValueError("clone result_id must be append-only and unique")
        return replace(
            self,
            revision=self.revision + 1,
            clone_results=self.clone_results + records,
        )
