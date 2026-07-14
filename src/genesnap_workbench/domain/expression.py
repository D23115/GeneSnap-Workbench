"""表达类项目的通用构建记录。"""

from __future__ import annotations

from dataclasses import dataclass, replace
from datetime import date, datetime
from enum import Enum

from .sequencing_submission import SequencingSubmissionRecord
from .project_tracking import ProjectInterruptionRecord


class ExpressionConstructKind(str, Enum):
    FULL_LENGTH = "full_length"
    TRUNCATION = "truncation"
    DELETION = "deletion"
    MUTATION = "mutation"


def _strict_dna(field_name: str, value: str) -> str:
    normalized = "".join(value.split()).upper()
    if normalized.startswith(">"):
        lines = value.lstrip("\ufeff").splitlines()
        normalized = "".join(
            line.strip() for line in lines if not line.lstrip().startswith(">")
        ).upper()
    if not normalized:
        raise ValueError(f"{field_name} must not be blank")
    invalid = sorted(set(normalized) - set("ACGT"))
    if invalid:
        raise ValueError(f"{field_name} contains unsupported bases: {', '.join(invalid)}")
    return normalized


@dataclass(frozen=True, slots=True)
class ExpressionDesignInput:
    project_id: str
    gene_symbol: str
    species: str
    source_cds: str
    construct_lines: tuple[str, ...]
    transcript_accession: str | None = None

    def __post_init__(self) -> None:
        if not self.project_id.strip() or not self.gene_symbol.strip():
            raise ValueError("project_id and gene_symbol must not be blank")
        object.__setattr__(self, "gene_symbol", self.gene_symbol.strip())
        object.__setattr__(self, "species", self.species.strip().lower())
        object.__setattr__(self, "source_cds", _strict_dna("source_cds", self.source_cds))
        object.__setattr__(
            self,
            "construct_lines",
            tuple(line.strip() for line in self.construct_lines if line.strip()),
        )
        if self.species not in {"human", "mouse", "rat"}:
            raise ValueError("species must be human, mouse, or rat")
        if not self.construct_lines:
            raise ValueError("construct_lines must not be empty")


@dataclass(frozen=True, slots=True)
class ExpressionVectorRules:
    protocol_version_id: str
    kozak_sequence: str
    stop_codon_rule: str
    c_terminal_fusion_name: str | None
    single_fragment_max_bp: int = 7000

    def __post_init__(self) -> None:
        object.__setattr__(
            self,
            "kozak_sequence",
            _strict_dna("kozak_sequence", self.kozak_sequence),
        )
        if self.stop_codon_rule not in {
            "preserve",
            "remove_for_c_terminal_fusion",
        }:
            raise ValueError("unsupported stop_codon_rule")
        if self.single_fragment_max_bp <= 0:
            raise ValueError("single_fragment_max_bp must be positive")


@dataclass(frozen=True, slots=True)
class ExpressionPCRFragment:
    fragment_no: int
    start: int
    end: int
    sequence: str

    def __post_init__(self) -> None:
        if self.fragment_no <= 0 or self.start < 0 or self.end <= self.start:
            raise ValueError("invalid expression PCR fragment")
        if len(self.sequence) != self.end - self.start:
            raise ValueError("fragment sequence does not match its interval")


@dataclass(frozen=True, slots=True)
class ExpressionMutation:
    notation: str
    amino_acid_position: int
    original_amino_acid: str
    new_amino_acid: str
    original_codon: str
    new_codon: str


@dataclass(frozen=True, slots=True)
class ExpressionConstruct:
    construct_id: str
    construct_name: str
    request_line: str
    kind: ExpressionConstructKind
    coding_sequence: str
    insert_sequence: str
    start_codon_reintroduced: bool
    terminal_stop_present: bool
    c_terminal_fusion_name: str | None
    mutations: tuple[ExpressionMutation, ...]
    fragments: tuple[ExpressionPCRFragment, ...]


@dataclass(frozen=True, slots=True)
class ExpressionDesignVersion:
    design_version_id: str
    project_id: str
    gene_symbol: str
    species: str
    transcript_accession: str | None
    source_cds_checksum: str
    protocol_version_id: str
    created_at: datetime
    constructs: tuple[ExpressionConstruct, ...]
    unparsed_lines: tuple[str, ...]
    design_warnings: tuple[str, ...]
    requires_confirmation: bool
    confirmation_history: tuple[ExpressionDesignConfirmation, ...] = ()
    gene_id: str | None = None


@dataclass(frozen=True, slots=True)
class ExpressionDesignConfirmation:
    confirmation_id: str
    reason: str
    actor: str
    occurred_at: datetime


@dataclass(frozen=True, slots=True)
class ExpressionPrimerPlan:
    name: str
    sequence: str
    direction: str
    fragment_no: int
    anneal_length: int
    overlap_length: int


@dataclass(frozen=True, slots=True)
class ExpressionConstructVectorPlan:
    construct_id: str
    construct_name: str
    forward_primer: str
    reverse_primer: str
    forward_anneal_length: int
    reverse_anneal_length: int
    expected_plasmid_sequence: str
    expected_plasmid_checksum: str
    left_junction: str
    right_junction: str
    primers: tuple[ExpressionPrimerPlan, ...]


@dataclass(frozen=True, slots=True)
class ExpressionVectorDesignResult:
    design_version_id: str
    vector_record_id: str
    vector_checksum: str
    protocol_version_id: str
    construct_plans: tuple[ExpressionConstructVectorPlan, ...]


@dataclass(frozen=True, slots=True)
class ExpressionAuditEvent:
    event_id: str
    event_type: str
    occurred_at: datetime
    actor: str
    from_status: str | None = None
    to_status: str | None = None
    note: str | None = None
    source: str = "user"


@dataclass(frozen=True, slots=True)
class ExpressionCloneResultRecord:
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
    premature_stop: bool = False
    manually_confirmed_usable: bool = False
    manual_review_status: str | None = None
    manual_note: str | None = None

    def __post_init__(self) -> None:
        if self.status not in {"pass", "fail", "warning"}:
            raise ValueError("clone result status must be pass, fail, or warning")
        if self.manual_review_status not in {None, "usable", "unusable"}:
            raise ValueError("manual_review_status must be usable, unusable, or None")


@dataclass(frozen=True, slots=True)
class ExpressionProjectSnapshot:
    project_id: str
    revision: int
    status: str
    active_design_version_id: str
    status_history: tuple[ExpressionAuditEvent, ...]
    clones_per_construct: int = 5
    clone_results: tuple[ExpressionCloneResultRecord, ...] = ()
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

    def append_status_event(
        self,
        event: ExpressionAuditEvent,
    ) -> ExpressionProjectSnapshot:
        if any(item.event_id == event.event_id for item in self.status_history):
            raise ValueError("status event_id must be unique")
        return replace(
            self,
            revision=self.revision + 1,
            status_history=self.status_history + (event,),
        )

    def append_clone_results(
        self,
        records: tuple[ExpressionCloneResultRecord, ...],
    ) -> ExpressionProjectSnapshot:
        known = {item.result_id for item in self.clone_results}
        incoming = [item.result_id for item in records]
        if len(incoming) != len(set(incoming)) or any(
            result_id in known for result_id in incoming
        ):
            raise ValueError("clone result_id must be append-only and unique")
        return replace(
            self,
            revision=self.revision + 1,
            clone_results=self.clone_results + records,
        )
