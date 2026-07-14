"""Validation and final-plasmid simulation for SYN vector protocols."""

from __future__ import annotations

from genesnap_workbench.domain.syn import SYNPlasmidFeature, SYNPlasmidSimulation
from genesnap_workbench.sequence_core.dna import normalize_dna, reverse_complement

from .models import (
    ProtocolValidationIssue,
    SiteRetentionRule,
    SYNVectorProtocol,
    SYNVectorProtocolValidationResult,
    VectorRecord,
    normalized_circular_checksum,
)


def circular_sequence_checksum(sequence: str) -> str:
    """Return a stable checksum independent of circular origin and strand."""
    return normalized_circular_checksum(sequence)


def _circular_site_positions(sequence: str, site: str) -> tuple[int, ...]:
    normalized = normalize_dna(sequence)
    motif = normalize_dna(site)
    circular = normalized + normalized[: len(motif) - 1]
    return tuple(
        start
        for start in range(len(normalized))
        if circular.startswith(motif, start)
    )


def _contains_site(sequence: str, site: str) -> bool:
    normalized = normalize_dna(sequence)
    motif = normalize_dna(site)
    reverse = reverse_complement(motif)
    return motif in normalized or reverse in normalized


def validate_syn_vector_protocol(
    protocol: SYNVectorProtocol,
    vector: VectorRecord,
    insert_sequence: str,
) -> SYNVectorProtocolValidationResult:
    """Validate exact vector binding, cut uniqueness, and release-site rules."""
    insert = normalize_dna(insert_sequence)
    errors: list[ProtocolValidationIssue] = []
    warnings: list[ProtocolValidationIssue] = []

    actual_vector_checksum = circular_sequence_checksum(vector.sequence)
    if actual_vector_checksum != vector.normalized_circular_sha256:
        errors.append(
            ProtocolValidationIssue(
                code="VECTOR_RECORD_HASH_INVALID",
                severity="error",
                message="载体记录中的序列与其自身校验值不一致",
            ),
        )
    if protocol.vector_record_id != vector.vector_record_id:
        errors.append(
            ProtocolValidationIssue(
                code="VECTOR_ID_MISMATCH",
                severity="error",
                message="protocol 绑定的载体记录与当前载体不一致",
            ),
        )
    if protocol.status != "enabled":
        errors.append(
            ProtocolValidationIssue(
                code="PROTOCOL_NOT_ENABLED",
                severity="error",
                message="只有已启用的 protocol 可以用于正式设计",
            ),
        )
    if protocol.workflow_type != "de_novo_gene_synthesis":
        errors.append(
            ProtocolValidationIssue(
                code="WORKFLOW_MISMATCH",
                severity="error",
                message="protocol 不适用于 SYN 工作流",
            ),
        )
    if protocol.vector_checksum != vector.normalized_circular_sha256:
        errors.append(
            ProtocolValidationIssue(
                code="VECTOR_HASH_MISMATCH",
                severity="error",
                message="protocol 绑定的载体序列校验值与当前载体不一致",
            ),
        )

    site_positions = _circular_site_positions(
        vector.sequence,
        protocol.linearization_site.sequence,
    )
    if not site_positions:
        errors.append(
            ProtocolValidationIssue(
                code="SITE_NOT_FOUND",
                severity="error",
                message="载体中未找到指定线性化位点",
            ),
        )
    elif len(site_positions) > 1:
        errors.append(
            ProtocolValidationIssue(
                code="SITE_NOT_UNIQUE",
                severity="error",
                message=f"指定线性化位点在载体中出现 {len(site_positions)} 次",
            ),
        )

    if protocol.site_retention_rule == SiteRetentionRule.REBUILD_FLANKING_SITES:
        if protocol.release_site is None:
            errors.append(
                ProtocolValidationIssue(
                    code="RELEASE_SITE_REQUIRED",
                    severity="error",
                    message="重建两侧位点时必须指定释放酶切位点",
                ),
            )
        elif _contains_site(insert, protocol.release_site.sequence):
            errors.append(
                ProtocolValidationIssue(
                    code="RELEASE_SITE_INTERNAL",
                    severity="error",
                    message="insert 内部含有同名释放位点，不能完整酶切释放",
                ),
            )
        elif _circular_site_positions(
            vector.sequence,
            protocol.release_site.sequence,
        ):
            errors.append(
                ProtocolValidationIssue(
                    code="RELEASE_SITE_IN_BACKBONE",
                    severity="error",
                    message="载体 backbone 内部含有同名释放位点，不能按两片段酶切解释",
                ),
            )

    return SYNVectorProtocolValidationResult(
        errors=tuple(errors),
        warnings=tuple(warnings),
    )


class SYNVectorProtocolError(ValueError):
    def __init__(self, validation: SYNVectorProtocolValidationResult) -> None:
        self.validation = validation
        super().__init__("；".join(issue.message for issue in validation.errors))


def simulate_syn_plasmid(
    vector: VectorRecord,
    protocol: SYNVectorProtocol,
    insert_sequence: str,
    *,
    design_version_id: str,
) -> SYNPlasmidSimulation:
    """Simulate one circular product from a validated SYN vector protocol."""
    insert = normalize_dna(insert_sequence)
    validation = validate_syn_vector_protocol(protocol, vector, insert)
    if not validation.is_valid:
        raise SYNVectorProtocolError(validation)

    site_start = _circular_site_positions(
        vector.sequence,
        protocol.linearization_site.sequence,
    )[0]
    cut_position = site_start + protocol.linearization_site.cut_offset
    linear_backbone = vector.sequence[cut_position:] + vector.sequence[:cut_position]

    cassette = insert
    digest_fragments: tuple[int, ...] = ()
    if protocol.site_retention_rule == SiteRetentionRule.REBUILD_FLANKING_SITES:
        assert protocol.release_site is not None
        release_sequence = protocol.release_site.sequence
        cassette = release_sequence + insert + release_sequence
        insert_fragment_length = len(insert) + len(release_sequence)
        digest_fragments = (
            insert_fragment_length,
            len(linear_backbone) + len(release_sequence),
        )

    arm_length = protocol.homology_arm_length
    left_arm = linear_backbone[-arm_length:]
    right_arm = linear_backbone[:arm_length]
    expected_sequence = linear_backbone + cassette
    release_length = (
        len(protocol.release_site.sequence)
        if protocol.site_retention_rule == SiteRetentionRule.REBUILD_FLANKING_SITES
        and protocol.release_site is not None
        else 0
    )
    insert_start = len(linear_backbone) + release_length
    insert_end = insert_start + len(insert)
    features = [
        SYNPlasmidFeature(
            label=vector.structural_display_name,
            feature_type="vector_backbone",
            start=0,
            end=len(linear_backbone),
        ),
        SYNPlasmidFeature(
            label="SYN insert",
            feature_type="misc_feature",
            start=insert_start,
            end=insert_end,
        ),
    ]
    if release_length:
        features.extend(
            (
                SYNPlasmidFeature(
                    label=f"{protocol.release_site.name} release site (left)",
                    feature_type="restriction_site",
                    start=len(linear_backbone),
                    end=insert_start,
                ),
                SYNPlasmidFeature(
                    label=f"{protocol.release_site.name} release site (right)",
                    feature_type="restriction_site",
                    start=insert_end,
                    end=insert_end + release_length,
                ),
            ),
        )
    junction_window = min(10, arm_length, len(cassette))
    junctions = (
        linear_backbone[-junction_window:] + cassette[:junction_window],
        cassette[-junction_window:] + linear_backbone[:junction_window],
    )
    return SYNPlasmidSimulation(
        design_version_id=design_version_id,
        vector_record_id=vector.vector_record_id,
        vector_checksum=vector.normalized_circular_sha256,
        protocol_version_id=protocol.protocol_version_id,
        linearization_sites=(protocol.linearization_site.name,),
        site_retention_rule=protocol.site_retention_rule.value,
        homology_arms=(left_arm, right_arm),
        junctions=junctions,
        expected_plasmid_sequence=expected_sequence,
        expected_plasmid_checksum=circular_sequence_checksum(expected_sequence),
        features=tuple(features),
        expected_digest_fragments_bp=digest_fragments,
    )
