"""Stable vector and protocol records used by construction workflows."""

from __future__ import annotations

from dataclasses import dataclass
from enum import Enum
import hashlib

from genesnap_workbench.sequence_core.dna import normalize_dna

from .comparison import canonical_circular_sequence


def normalized_circular_checksum(sequence: str) -> str:
    canonical = canonical_circular_sequence(sequence)
    return hashlib.sha256(canonical.encode("ascii")).hexdigest()


class SiteRetentionRule(str, Enum):
    NOT_REQUIRED = "not_required"
    REBUILD_FLANKING_SITES = "rebuild_flanking_sites"


@dataclass(frozen=True, slots=True)
class RestrictionSite:
    name: str
    sequence: str
    cut_offset: int

    def __post_init__(self) -> None:
        normalized = normalize_dna(self.sequence)
        object.__setattr__(self, "sequence", normalized)
        if not 0 < self.cut_offset < len(normalized):
            raise ValueError("Restriction-site cut offset must be inside recognition site")


@dataclass(frozen=True, slots=True)
class VectorRecord:
    vector_record_id: str
    structural_display_name: str
    sequence: str
    topology: str
    normalized_circular_sha256: str
    local_aliases: tuple[str, ...] = ()
    backbone_family: str = "unknown"
    public_source_urls: tuple[str, ...] = ()
    public_equivalence_status: str = "unknown"

    @classmethod
    def from_sequence(
        cls,
        *,
        vector_record_id: str,
        structural_display_name: str,
        sequence: str,
        topology: str = "circular",
        local_aliases: tuple[str, ...] = (),
        backbone_family: str = "unknown",
        public_source_urls: tuple[str, ...] = (),
        public_equivalence_status: str = "unknown",
    ) -> VectorRecord:
        normalized = normalize_dna(sequence)
        return cls(
            vector_record_id=vector_record_id,
            structural_display_name=structural_display_name,
            sequence=normalized,
            topology=topology,
            normalized_circular_sha256=normalized_circular_checksum(normalized),
            local_aliases=local_aliases,
            backbone_family=backbone_family,
            public_source_urls=public_source_urls,
            public_equivalence_status=public_equivalence_status,
        )


@dataclass(frozen=True, slots=True)
class SYNVectorProtocol:
    protocol_id: str
    protocol_version_id: str
    display_name: str
    status: str
    experimental_validation_status: str
    vector_record_id: str
    vector_checksum: str
    workflow_type: str
    insertion_mode: str
    linearization_site: RestrictionSite
    site_retention_rule: SiteRetentionRule
    release_site: RestrictionSite | None
    homology_arm_length: int


@dataclass(frozen=True, slots=True)
class ShRNAVectorProtocol:
    protocol_id: str
    protocol_version_id: str
    display_name: str
    status: str
    experimental_validation_status: str
    vector_record_id: str
    vector_checksum: str
    workflow_type: str
    insertion_mode: str
    left_site: RestrictionSite
    right_site: RestrictionSite
    sequencing_primer_name: str
    default_target_count: int = 3
    default_clones_per_target: int = 5

    def __post_init__(self) -> None:
        if not 1 <= self.default_target_count <= 3:
            raise ValueError("default_target_count must be between 1 and 3")
        if not 1 <= self.default_clones_per_target <= 96:
            raise ValueError("default_clones_per_target must be between 1 and 96")


@dataclass(frozen=True, slots=True)
class ExpressionVectorProtocol:
    protocol_id: str
    protocol_version_id: str
    display_name: str
    status: str
    experimental_validation_status: str
    vector_record_id: str
    vector_checksum: str
    workflow_type: str
    insertion_mode: str
    left_boundary: int
    right_boundary: int
    left_primer_homology: str
    right_primer_homology: str
    kozak_sequence: str
    stop_codon_rule: str
    c_terminal_fusion_name: str | None
    single_fragment_max_bp: int = 7000
    anneal_min_bp: int = 18
    anneal_max_bp: int = 25

    def __post_init__(self) -> None:
        for field_name in (
            "left_primer_homology",
            "right_primer_homology",
            "kozak_sequence",
        ):
            object.__setattr__(self, field_name, normalize_dna(getattr(self, field_name)))
        if self.left_boundary < 0 or self.right_boundary < self.left_boundary:
            raise ValueError("表达载体插入边界无效：右边界不能小于左边界")
        if not 1 <= self.anneal_min_bp <= self.anneal_max_bp:
            raise ValueError("表达载体退火长度范围无效：最小长度必须不大于最大长度")


@dataclass(frozen=True, slots=True)
class ReporterVectorProtocol:
    protocol_id: str
    protocol_version_id: str
    display_name: str
    status: str
    experimental_validation_status: str
    vector_record_id: str
    vector_checksum: str
    workflow_type: str
    insertion_mode: str
    left_boundary: int
    right_boundary: int
    left_primer_homology: str
    right_primer_homology: str
    default_sequencing_method: str = "Nanopore"
    anneal_min_bp: int = 18
    anneal_max_bp: int = 25

    def __post_init__(self) -> None:
        object.__setattr__(
            self,
            "left_primer_homology",
            normalize_dna(self.left_primer_homology),
        )
        object.__setattr__(
            self,
            "right_primer_homology",
            normalize_dna(self.right_primer_homology),
        )
        if self.left_boundary < 0 or self.right_boundary < self.left_boundary:
            raise ValueError("invalid reporter insertion boundaries")
        if not 1 <= self.anneal_min_bp <= self.anneal_max_bp:
            raise ValueError("invalid reporter annealing length range")


@dataclass(frozen=True, slots=True)
class ProtocolValidationIssue:
    code: str
    severity: str
    message: str


@dataclass(frozen=True, slots=True)
class SYNVectorProtocolValidationResult:
    errors: tuple[ProtocolValidationIssue, ...]
    warnings: tuple[ProtocolValidationIssue, ...]

    @property
    def is_valid(self) -> bool:
        return not self.errors

    @property
    def error_codes(self) -> tuple[str, ...]:
        return tuple(issue.code for issue in self.errors)


@dataclass(frozen=True, slots=True)
class ShRNAVectorProtocolValidationResult:
    errors: tuple[ProtocolValidationIssue, ...]
    warnings: tuple[ProtocolValidationIssue, ...]

    @property
    def is_valid(self) -> bool:
        return not self.errors

    @property
    def error_codes(self) -> tuple[str, ...]:
        return tuple(issue.code for issue in self.errors)


@dataclass(frozen=True, slots=True)
class ExpressionVectorProtocolValidationResult:
    errors: tuple[ProtocolValidationIssue, ...]
    warnings: tuple[ProtocolValidationIssue, ...]

    @property
    def is_valid(self) -> bool:
        return not self.errors

    @property
    def error_codes(self) -> tuple[str, ...]:
        return tuple(issue.code for issue in self.errors)


@dataclass(frozen=True, slots=True)
class ReporterVectorProtocolValidationResult:
    errors: tuple[ProtocolValidationIssue, ...]
    warnings: tuple[ProtocolValidationIssue, ...]

    @property
    def is_valid(self) -> bool:
        return not self.errors

    @property
    def error_codes(self) -> tuple[str, ...]:
        return tuple(issue.code for issue in self.errors)
