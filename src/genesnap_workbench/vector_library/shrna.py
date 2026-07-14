"""pLKO/shRNA protocol 校验与预期质粒模拟。"""

from __future__ import annotations

from genesnap_workbench.domain.shrna import ShRNAPlasmidSimulation
from genesnap_workbench.sequence_core.dna import normalize_dna

from .models import (
    ProtocolValidationIssue,
    ShRNAVectorProtocol,
    ShRNAVectorProtocolValidationResult,
    VectorRecord,
    normalized_circular_checksum,
)


class ShRNAVectorProtocolError(ValueError):
    def __init__(self, validation: ShRNAVectorProtocolValidationResult) -> None:
        self.validation = validation
        message = "; ".join(issue.message for issue in validation.errors)
        super().__init__(message or "shRNA vector protocol is not ready")


def _site_positions(sequence: str, motif: str) -> tuple[int, ...]:
    positions: list[int] = []
    start = 0
    while True:
        position = sequence.find(motif, start)
        if position < 0:
            return tuple(positions)
        positions.append(position)
        start = position + 1


def validate_shrna_protocol(
    vector: VectorRecord,
    protocol: ShRNAVectorProtocol,
) -> ShRNAVectorProtocolValidationResult:
    errors: list[ProtocolValidationIssue] = []
    warnings: list[ProtocolValidationIssue] = []
    if protocol.vector_record_id != vector.vector_record_id:
        errors.append(
            ProtocolValidationIssue(
                "vector_id_mismatch",
                "error",
                "protocol 绑定的载体 ID 与当前载体不一致",
            ),
        )
    if protocol.vector_checksum != vector.normalized_circular_sha256:
        errors.append(
            ProtocolValidationIssue(
                "vector_checksum_mismatch",
                "error",
                "protocol 绑定的载体序列校验值与当前载体不一致",
            ),
        )
    if protocol.workflow_type != "shrna_knockdown":
        errors.append(
            ProtocolValidationIssue(
                "workflow_mismatch",
                "error",
                "protocol 不适用于 shRNA 工作流",
            ),
        )

    left_positions = _site_positions(vector.sequence, protocol.left_site.sequence)
    right_positions = _site_positions(vector.sequence, protocol.right_site.sequence)
    if len(left_positions) != 1:
        errors.append(
            ProtocolValidationIssue(
                "left_site_not_unique",
                "error",
                f"{protocol.left_site.name} 位点数量为 {len(left_positions)}，无法唯一定位",
            ),
        )
    if len(right_positions) != 1:
        errors.append(
            ProtocolValidationIssue(
                "right_site_not_unique",
                "error",
                f"{protocol.right_site.name} 位点数量为 {len(right_positions)}，无法唯一定位",
            ),
        )
    if len(left_positions) == 1 and len(right_positions) == 1:
        left_cut = left_positions[0] + protocol.left_site.cut_offset
        right_cut = right_positions[0] + protocol.right_site.cut_offset
        if left_cut >= right_cut:
            errors.append(
                ProtocolValidationIssue(
                    "site_order_invalid",
                    "error",
                    "AgeI/EcoRI 切口顺序无效",
                ),
            )
    if protocol.experimental_validation_status != "verified":
        warnings.append(
            ProtocolValidationIssue(
                "protocol_unverified",
                "warning",
                "该 protocol 尚未绑定湿实验验证记录",
            ),
        )
    return ShRNAVectorProtocolValidationResult(tuple(errors), tuple(warnings))


def simulate_shrna_plasmid(
    vector: VectorRecord,
    protocol: ShRNAVectorProtocol,
    *,
    target_id: str,
    forward_oligo: str,
) -> ShRNAPlasmidSimulation:
    validation = validate_shrna_protocol(vector, protocol)
    if protocol.status != "enabled":
        validation = ShRNAVectorProtocolValidationResult(
            errors=validation.errors
            + (
                ProtocolValidationIssue(
                    "protocol_not_enabled",
                    "error",
                    "载体序列尚未确认，不能生成正式预期质粒",
                ),
            ),
            warnings=validation.warnings,
        )
    if not validation.is_valid:
        raise ShRNAVectorProtocolError(validation)

    insert = normalize_dna(forward_oligo)
    if not (
        insert.startswith("CCGG")
        and "CTCGAG" in insert
        and insert.endswith("TTTTTG")
    ):
        invalid = ProtocolValidationIssue(
            "invalid_plko_forward_oligo",
            "error",
            "正向 oligo 不符合已确认的 pLKO hairpin 结构",
        )
        raise ShRNAVectorProtocolError(
            ShRNAVectorProtocolValidationResult((invalid,), validation.warnings),
        )

    left_start = vector.sequence.index(protocol.left_site.sequence)
    right_start = vector.sequence.index(protocol.right_site.sequence)
    left_cut = left_start + protocol.left_site.cut_offset
    right_cut = right_start + protocol.right_site.cut_offset
    expected = vector.sequence[:left_cut] + insert + vector.sequence[right_cut:]
    return ShRNAPlasmidSimulation(
        target_id=target_id,
        vector_record_id=vector.vector_record_id,
        vector_checksum=vector.normalized_circular_sha256,
        protocol_version_id=protocol.protocol_version_id,
        expected_plasmid_sequence=expected,
        expected_plasmid_checksum=normalized_circular_checksum(expected),
        left_cut_position=left_cut,
        right_cut_position=right_cut,
        left_site_count=expected.count(protocol.left_site.sequence),
        right_site_count=expected.count(protocol.right_site.sequence),
    )
