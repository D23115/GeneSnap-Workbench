"""Immutable domain records for the de novo gene synthesis workflow."""

from __future__ import annotations

from dataclasses import dataclass, replace
from datetime import datetime
from decimal import Decimal
from enum import Enum


class SYNRoute(str, Enum):
    SINGLE_POOL = "single_pool"
    MODULAR = "modular"


class ResuspensionStatus(str, Enum):
    MISSING = "missing"
    PARTIAL = "partial"
    COMPLETE = "complete"


class SYNAssemblyAttemptResult(str, Enum):
    PENDING = "pending"
    SUCCESS = "success"
    FAILED = "failed"
    PARTIAL = "partial"


class SYNColonyPCRResult(str, Enum):
    PENDING = "pending"
    POSITIVE = "positive"
    NEGATIVE = "negative"
    UNCERTAIN = "uncertain"


class PlasmidPrepStatus(str, Enum):
    PENDING = "pending"
    IN_PROGRESS = "in_progress"
    COMPLETED = "completed"


class SYNSequencingResult(str, Enum):
    CORRECT = "correct"
    INCORRECT = "incorrect"
    UNCERTAIN = "uncertain"


def _require_nonblank(field_name: str, value: str) -> None:
    if not value.strip():
        raise ValueError(f"{field_name} must not be blank")


def _validate_interval(start: int, end: int) -> None:
    if start < 0 or end <= start:
        raise ValueError(f"Invalid half-open interval: [{start}, {end})")


@dataclass(frozen=True, slots=True)
class SYNDesignInput:
    project_id: str
    target_name: str
    raw_sequence: str
    input_format: str
    vector_protocol_version_id: str
    gene_symbol: str | None = None
    species: str | None = None
    site_retention_rule: str = "not_required"
    reagent_profile_ids: tuple[str, ...] = ()

    def __post_init__(self) -> None:
        for field_name in (
            "project_id",
            "target_name",
            "raw_sequence",
            "input_format",
            "vector_protocol_version_id",
        ):
            _require_nonblank(field_name, getattr(self, field_name))


@dataclass(frozen=True, slots=True)
class SYNQCRisk:
    rule_key: str
    severity: str
    start: int
    end: int
    observed_value: str
    message: str
    requires_confirmation: bool

    def __post_init__(self) -> None:
        _validate_interval(self.start, self.end)


@dataclass(frozen=True, slots=True)
class SYNSequenceQCResult:
    design_version_id: str
    rules_version: str
    sequence_checksum: str
    sequence_length: int
    overall_gc_percent: Decimal
    risks: tuple[SYNQCRisk, ...]
    blocked_reasons: tuple[str, ...]
    confirmable_warnings: tuple[str, ...]


@dataclass(frozen=True, slots=True)
class SYNThermodynamicMetadata:
    analyzer_name: str
    analyzer_version: str
    tm_celsius: Decimal
    parameters: tuple[tuple[str, str], ...] = ()


@dataclass(frozen=True, slots=True)
class SYNAssemblyOligo:
    oligo_id: str
    design_version_id: str
    name: str
    sequence: str
    strand: str
    start: int
    end: int
    pool_id: str
    module_id: str
    overlap_left: tuple[int, int] | None
    overlap_right: tuple[int, int] | None
    tm_metadata: SYNThermodynamicMetadata

    def __post_init__(self) -> None:
        if len(self.sequence) > 65:
            raise ValueError("Assembly oligo must not exceed 65 nt")
        _validate_interval(self.start, self.end)


@dataclass(frozen=True, slots=True)
class SYNModule:
    design_version_id: str
    module_id: str
    ordinal: int
    start: int
    end: int
    sequence_checksum: str
    oligo_ids: tuple[str, ...]
    boundary_reason: str
    left_overlap: tuple[int, int] | None = None
    right_overlap: tuple[int, int] | None = None

    def __post_init__(self) -> None:
        _validate_interval(self.start, self.end)
        for overlap in (self.left_overlap, self.right_overlap):
            if overlap is not None:
                _validate_interval(*overlap)


@dataclass(frozen=True, slots=True)
class SYNModulePlan:
    design_version_id: str
    route: SYNRoute
    modules: tuple[SYNModule, ...]
    requires_confirmation: bool
    routing_reason: str

    def __post_init__(self) -> None:
        if not self.modules:
            raise ValueError("Module plan must contain at least one module")
        previous_end = 0
        for expected_ordinal, module in enumerate(self.modules, start=1):
            if module.design_version_id != self.design_version_id:
                raise ValueError("Module design_version_id does not match module plan")
            if module.ordinal != expected_ordinal or module.start != previous_end:
                raise ValueError("Modules must be ordered and contiguous")
            previous_end = module.end


@dataclass(frozen=True, slots=True)
class SYNPlasmidFeature:
    label: str
    feature_type: str
    start: int
    end: int
    strand: int = 1

    def __post_init__(self) -> None:
        _validate_interval(self.start, self.end)


@dataclass(frozen=True, slots=True)
class SYNPlasmidSimulation:
    design_version_id: str
    vector_record_id: str
    vector_checksum: str
    protocol_version_id: str
    linearization_sites: tuple[str, ...]
    site_retention_rule: str
    homology_arms: tuple[str, str]
    junctions: tuple[str, ...]
    expected_plasmid_sequence: str
    expected_plasmid_checksum: str
    features: tuple[SYNPlasmidFeature, ...] = ()
    expected_digest_fragments_bp: tuple[int, ...] = ()


@dataclass(frozen=True, slots=True)
class SYNAuditEvent:
    event_id: str
    event_type: str
    occurred_at: datetime
    actor: str
    note: str | None = None
    from_status: str | None = None
    to_status: str | None = None
    source: str = "user"
    related_entity_id: str | None = None


@dataclass(frozen=True, slots=True)
class SYNManualOverrideRecord:
    override_id: str
    field_path: str
    old_value: str
    new_value: str
    reason: str
    occurred_at: datetime
    actor: str


@dataclass(frozen=True, slots=True)
class SYNDesignVersion:
    design_version_id: str
    project_id: str
    version_no: int
    created_at: datetime
    raw_sequence_checksum: str
    normalized_sequence: str
    normalized_checksum: str
    final_sequence: str
    final_checksum: str
    qc_result: SYNSequenceQCResult
    module_plan: SYNModulePlan
    oligos: tuple[SYNAssemblyOligo, ...]
    plasmid_simulation: SYNPlasmidSimulation
    rule_versions: tuple[str, ...]
    manual_overrides: tuple[SYNManualOverrideRecord, ...]
    design_warnings: tuple[str, ...] = ()
    requires_confirmation: bool = False

    def __post_init__(self) -> None:
        nested_version_ids = (
            self.qc_result.design_version_id,
            self.module_plan.design_version_id,
            self.plasmid_simulation.design_version_id,
            *(oligo.design_version_id for oligo in self.oligos),
        )
        if any(item != self.design_version_id for item in nested_version_ids):
            raise ValueError("Nested design_version_id does not match design version")
        if self.qc_result.sequence_checksum != self.final_checksum:
            raise ValueError("QC sequence checksum must match final_checksum")
        if self.module_plan.modules[-1].end != len(self.final_sequence):
            raise ValueError("Module plan must cover the complete final sequence")

        known_oligo_ids = {oligo.oligo_id for oligo in self.oligos}
        for module in self.module_plan.modules:
            missing = set(module.oligo_ids) - known_oligo_ids
            if missing:
                names = ", ".join(sorted(missing))
                raise ValueError(f"Module references missing oligo: {names}")


@dataclass(frozen=True, slots=True)
class OligoResuspensionItem:
    oligo_id: str
    planned_amount_nmol: Decimal | None
    actual_amount_nmol: Decimal | None
    target_stock_concentration_uM: Decimal
    water_volume_ul: Decimal | None
    actual_stock_concentration_uM: Decimal | None


@dataclass(frozen=True, slots=True)
class OligoResuspensionPlan:
    design_version_id: str
    items: tuple[OligoResuspensionItem, ...]
    status: ResuspensionStatus
    generated_at: datetime


@dataclass(frozen=True, slots=True)
class OligoMixItem:
    pool_id: str
    oligo_id: str
    reference_concentration_uM: Decimal
    actual_concentration_uM: Decimal | None
    sample_volume_ul: Decimal | None


@dataclass(frozen=True, slots=True)
class OligoMixPlan:
    design_version_id: str
    standard_volume_per_oligo_ul: Decimal | None
    items: tuple[OligoMixItem, ...]
    generated_at: datetime

    @property
    def is_formal_export_ready(self) -> bool:
        return self.standard_volume_per_oligo_ul is not None and all(
            item.sample_volume_ul is not None for item in self.items
        )


@dataclass(frozen=True, slots=True)
class SYNAssemblyAttemptRecord:
    attempt_id: str
    project_id: str
    design_version_id: str
    syn_assembly_round_no: int
    restart_from_substep: str
    result: SYNAssemblyAttemptResult
    started_at: datetime
    completed_at: datetime | None = None
    note: str | None = None


@dataclass(frozen=True, slots=True)
class SYNAssemblyStepRecord:
    record_id: str
    attempt_id: str
    substep: str
    step_attempt_no: int
    result: SYNAssemblyAttemptResult
    recorded_at: datetime
    actor: str
    note: str | None = None

    def __post_init__(self) -> None:
        if self.step_attempt_no <= 0:
            raise ValueError("step_attempt_no must be positive")


@dataclass(frozen=True, slots=True)
class SYNColonyPCRRecord:
    clone_id: str
    attempt_id: str
    clone_no: int
    display_name: str
    result: SYNColonyPCRResult
    expected_band_bp: int | None = None
    observed_note: str | None = None
    recorded_at: datetime | None = None


@dataclass(frozen=True, slots=True)
class PlasmidPrepRecord:
    clone_id: str
    selected_for_prep: bool
    status: PlasmidPrepStatus
    completed_at: datetime | None = None
    note: str | None = None


@dataclass(frozen=True, slots=True)
class SYNSequencingConfirmation:
    confirmation_id: str
    clone_id: str
    result: SYNSequencingResult
    confirmed_at: datetime
    actor: str
    note: str | None = None
    supersedes_confirmation_id: str | None = None


@dataclass(frozen=True, slots=True)
class SYNProjectSnapshot:
    project_id: str
    revision: int
    status: str
    resuspension_data_status: str
    syn_assembly_round_no: int
    syn_assembly_substep: str | None
    active_design_version_id: str
    attempts: tuple[SYNAssemblyAttemptRecord, ...]
    colonies: tuple[SYNColonyPCRRecord, ...]
    prep_records: tuple[PlasmidPrepRecord, ...]
    sequencing_confirmations: tuple[SYNSequencingConfirmation, ...]
    status_history: tuple[SYNAuditEvent, ...]
    manual_override_history: tuple[SYNManualOverrideRecord, ...]
    step_records: tuple[SYNAssemblyStepRecord, ...] = ()
    actual_completed_at: datetime | None = None

    def __post_init__(self) -> None:
        for attempt in self.attempts:
            self._validate_attempt_owner(attempt)
        attempt_ids = {attempt.attempt_id for attempt in self.attempts}
        for record in self.step_records:
            if record.attempt_id not in attempt_ids:
                raise ValueError("Assembly step record references missing attempt")
        record_ids = [record.record_id for record in self.step_records]
        if len(record_ids) != len(set(record_ids)):
            raise ValueError("Assembly step record_id must be unique")

    def _validate_attempt_owner(self, attempt: SYNAssemblyAttemptRecord) -> None:
        if attempt.project_id != self.project_id:
            raise ValueError("Assembly attempt project_id does not match snapshot")
        if attempt.design_version_id != self.active_design_version_id:
            raise ValueError(
                "Assembly attempt design_version_id does not match active design",
            )

    def append_attempt(self, attempt: SYNAssemblyAttemptRecord) -> SYNProjectSnapshot:
        """Return a new snapshot with one additional assembly attempt."""
        self._validate_attempt_owner(attempt)
        if any(item.attempt_id == attempt.attempt_id for item in self.attempts):
            raise ValueError(f"Assembly attempt already exists: {attempt.attempt_id}")
        if attempt.syn_assembly_round_no <= self.syn_assembly_round_no:
            raise ValueError("New assembly attempt must advance the assembly round")
        return replace(
            self,
            revision=self.revision + 1,
            syn_assembly_round_no=max(
                self.syn_assembly_round_no,
                attempt.syn_assembly_round_no,
            ),
            attempts=self.attempts + (attempt,),
        )

    def append_status_event(self, event: SYNAuditEvent) -> SYNProjectSnapshot:
        """Return a new snapshot with one additional status event."""
        if any(item.event_id == event.event_id for item in self.status_history):
            raise ValueError(f"Status event already exists: {event.event_id}")
        return replace(
            self,
            revision=self.revision + 1,
            status_history=self.status_history + (event,),
        )
