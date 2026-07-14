"""pLKO/shRNA 工作流的不可变领域记录。"""

from __future__ import annotations

from dataclasses import dataclass, replace
from datetime import date, datetime
from decimal import Decimal
from enum import Enum

from .sequencing_submission import SequencingSubmissionRecord
from .project_tracking import ProjectInterruptionRecord

class BlastScreenStatus(str, Enum):
    PENDING = "pending"
    PASS = "pass"
    FAIL = "fail"
    UNAVAILABLE = "unavailable"
    MANUALLY_ACCEPTED = "manually_accepted"


def _require_nonblank(field_name: str, value: str) -> None:
    if not value.strip():
        raise ValueError(f"{field_name} must not be blank")


def _normalize_strict_dna(field_name: str, value: str) -> str:
    normalized = "".join(value.split()).upper()
    if not normalized:
        raise ValueError(f"{field_name} must not be blank")
    invalid = sorted(set(normalized) - set("ACGT"))
    if invalid:
        raise ValueError(
            f"{field_name} contains unsupported DNA characters: {', '.join(invalid)}",
        )
    return normalized


@dataclass(frozen=True, slots=True)
class ShRNADesignInput:
    project_id: str
    gene_symbol: str
    species: str
    cds_sequence: str
    vector_protocol_version_id: str
    transcript_accession: str | None = None
    gene_id: str | None = None
    ccds_id: str | None = None
    target_count: int = 3
    clones_per_target: int = 5

    def __post_init__(self) -> None:
        for field_name in (
            "project_id",
            "gene_symbol",
            "species",
            "cds_sequence",
            "vector_protocol_version_id",
        ):
            _require_nonblank(field_name, getattr(self, field_name))
        object.__setattr__(
            self,
            "cds_sequence",
            _normalize_strict_dna("cds_sequence", self.cds_sequence),
        )
        object.__setattr__(self, "gene_symbol", self.gene_symbol.strip())
        object.__setattr__(self, "species", self.species.strip().lower())
        if self.species not in {"human", "mouse", "rat"}:
            raise ValueError("species must be human, mouse, or rat")
        if not 1 <= self.target_count <= 3:
            raise ValueError("target_count must be between 1 and 3")
        if not 1 <= self.clones_per_target <= 96:
            raise ValueError("clones_per_target must be between 1 and 96")


@dataclass(frozen=True, slots=True)
class ShRNACandidate:
    candidate_id: str
    target_sequence: str
    start_position: int | None
    intrinsic_score: Decimal
    source_rank: int
    blast_status: BlastScreenStatus = BlastScreenStatus.PENDING
    first_offtarget_gene: str | None = None
    first_offtarget_mismatches: int | None = None
    blast_note: str | None = None
    forward_oligo_sequence: str | None = None
    reverse_oligo_sequence: str | None = None
    oligo_source: str | None = None
    oligo_comparison_note: str | None = None

    @property
    def gc_percent(self) -> Decimal:
        gc_count = self.target_sequence.count("G") + self.target_sequence.count("C")
        return (Decimal(gc_count) * Decimal("100") / Decimal(len(self.target_sequence))).quantize(
            Decimal("0.01"),
        )

    def __post_init__(self) -> None:
        _require_nonblank("candidate_id", self.candidate_id)
        normalized = _normalize_strict_dna("target_sequence", self.target_sequence)
        if not 15 <= len(normalized) <= 30:
            raise ValueError("target_sequence must contain 15-30 nt")
        object.__setattr__(self, "target_sequence", normalized)
        if self.start_position is not None and self.start_position < 0:
            raise ValueError("start_position must not be negative")
        if self.source_rank <= 0:
            raise ValueError("source_rank must be positive")
        if self.first_offtarget_mismatches is not None and self.first_offtarget_mismatches < 0:
            raise ValueError("first_offtarget_mismatches must not be negative")
        for field_name in ("forward_oligo_sequence", "reverse_oligo_sequence"):
            value = getattr(self, field_name)
            if value is not None:
                object.__setattr__(self, field_name, _normalize_strict_dna(field_name, value))
        if (self.forward_oligo_sequence is None) != (self.reverse_oligo_sequence is None):
            raise ValueError("Broad oligo sequences must be provided as a pair")


@dataclass(frozen=True, slots=True)
class ShRNAOligoPair:
    target_id: str
    forward_name: str
    forward_sequence: str
    reverse_name: str
    reverse_sequence: str

    def __post_init__(self) -> None:
        for field_name in ("target_id", "forward_name", "reverse_name"):
            _require_nonblank(field_name, getattr(self, field_name))
        object.__setattr__(
            self,
            "forward_sequence",
            _normalize_strict_dna("forward_sequence", self.forward_sequence),
        )
        object.__setattr__(
            self,
            "reverse_sequence",
            _normalize_strict_dna("reverse_sequence", self.reverse_sequence),
        )


@dataclass(frozen=True, slots=True)
class ShRNATargetDesign:
    target_id: str
    target_no: int
    candidate: ShRNACandidate
    oligos: ShRNAOligoPair
    clone_names: tuple[str, ...]

    def __post_init__(self) -> None:
        _require_nonblank("target_id", self.target_id)
        if self.target_no <= 0:
            raise ValueError("target_no must be positive")
        if self.oligos.target_id != self.target_id:
            raise ValueError("oligo pair belongs to a different target")
        if not self.clone_names:
            raise ValueError("target must contain at least one clone name")
        if len(self.clone_names) != len(set(self.clone_names)):
            raise ValueError("clone names must be unique")


@dataclass(frozen=True, slots=True)
class ShRNATargetSelectionResult:
    selected: tuple[ShRNACandidate, ...]
    remaining: tuple[ShRNACandidate, ...]
    spacing_relaxed: bool
    minimum_spacing_bp: int


@dataclass(frozen=True, slots=True)
class ShRNABlastResolution:
    selected: tuple[ShRNACandidate, ...]
    needs_screening: tuple[ShRNACandidate, ...]
    completed: bool
    requires_confirmation: bool
    note: str


@dataclass(frozen=True, slots=True)
class ShRNADesignVersion:
    design_version_id: str
    project_id: str
    gene_symbol: str
    species: str
    vector_record_id: str
    vector_checksum: str
    vector_protocol_version_id: str
    created_at: datetime
    cds_checksum: str
    transcript_accession: str | None
    gene_id: str | None
    ccds_id: str | None
    target_count: int
    clones_per_target: int
    targets: tuple[ShRNATargetDesign, ...]
    plasmid_simulations: tuple[ShRNAPlasmidSimulation, ...]
    rule_versions: tuple[str, ...]
    design_warnings: tuple[str, ...] = ()
    requires_confirmation: bool = False

    def __post_init__(self) -> None:
        if len(self.targets) != self.target_count:
            raise ValueError("target count does not match the design input")
        for expected_no, target in enumerate(self.targets, start=1):
            if target.target_no != expected_no:
                raise ValueError("targets must use contiguous target numbers")
            if len(target.clone_names) != self.clones_per_target:
                raise ValueError("clone count does not match the design input")
        if len(self.plasmid_simulations) != self.target_count:
            raise ValueError("each shRNA target must have one plasmid simulation")
        target_ids = {target.target_id for target in self.targets}
        simulation_ids = {item.target_id for item in self.plasmid_simulations}
        if target_ids != simulation_ids:
            raise ValueError("plasmid simulations do not match shRNA targets")
        if any(
            item.vector_record_id != self.vector_record_id
            or item.vector_checksum != self.vector_checksum
            or item.protocol_version_id != self.vector_protocol_version_id
            for item in self.plasmid_simulations
        ):
            raise ValueError("plasmid simulations are not bound to the design vector")


@dataclass(frozen=True, slots=True)
class ShRNAPlasmidSimulation:
    target_id: str
    vector_record_id: str
    vector_checksum: str
    protocol_version_id: str
    expected_plasmid_sequence: str
    expected_plasmid_checksum: str
    left_cut_position: int
    right_cut_position: int
    left_site_count: int
    right_site_count: int

    def __post_init__(self) -> None:
        _require_nonblank("target_id", self.target_id)
        if self.left_cut_position < 0 or self.right_cut_position <= self.left_cut_position:
            raise ValueError("invalid shRNA insertion cut positions")
        object.__setattr__(
            self,
            "expected_plasmid_sequence",
            _normalize_strict_dna(
                "expected_plasmid_sequence",
                self.expected_plasmid_sequence,
            ),
        )


@dataclass(frozen=True, slots=True)
class ShRNAAuditEvent:
    event_id: str
    event_type: str
    occurred_at: datetime
    actor: str
    from_status: str | None = None
    to_status: str | None = None
    note: str | None = None
    source: str = "user"


@dataclass(frozen=True, slots=True)
class ShRNACloneResultRecord:
    result_id: str
    clone_name: str
    target_id: str
    status: str
    reason: str
    analyzed_at: datetime
    source_files: tuple[str, ...]
    match_start: int | None = None
    manually_confirmed_usable: bool = False
    manual_review_status: str | None = None
    manual_note: str | None = None

    def __post_init__(self) -> None:
        if self.status not in {"pass", "fail", "warning"}:
            raise ValueError("clone result status must be pass, fail, or warning")
        if self.manual_review_status not in {None, "usable", "unusable"}:
            raise ValueError("manual_review_status must be usable, unusable, or None")


@dataclass(frozen=True, slots=True)
class ShRNAProjectSnapshot:
    project_id: str
    revision: int
    status: str
    active_design_version_id: str
    clone_results: tuple[ShRNACloneResultRecord, ...]
    status_history: tuple[ShRNAAuditEvent, ...]
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

    def append_clone_results(
        self,
        records: tuple[ShRNACloneResultRecord, ...],
    ) -> ShRNAProjectSnapshot:
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

    def append_status_event(self, event: ShRNAAuditEvent) -> ShRNAProjectSnapshot:
        if any(item.event_id == event.event_id for item in self.status_history):
            raise ValueError("status event_id must be unique")
        return replace(
            self,
            revision=self.revision + 1,
            status_history=self.status_history + (event,),
        )
