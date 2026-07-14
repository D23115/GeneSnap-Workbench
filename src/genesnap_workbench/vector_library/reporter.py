"""Promoter-luciferase 设计与数据化载体 protocol 的连接层。"""

from __future__ import annotations

from genesnap_workbench.domain.reporter import (
    ReporterConstructVectorPlan,
    ReporterDesignVersion,
    ReporterPrimerPlan,
    ReporterVectorDesignResult,
)
from genesnap_workbench.sequence_core.dna import reverse_complement

from .models import (
    ProtocolValidationIssue,
    ReporterVectorProtocol,
    ReporterVectorProtocolValidationResult,
    VectorRecord,
    normalized_circular_checksum,
)


class ReporterVectorProtocolError(ValueError):
    def __init__(self, validation: ReporterVectorProtocolValidationResult) -> None:
        self.validation = validation
        super().__init__("；".join(item.message for item in validation.errors))


def validate_reporter_protocol(
    vector: VectorRecord,
    protocol: ReporterVectorProtocol,
) -> ReporterVectorProtocolValidationResult:
    errors = []
    warnings = []
    if protocol.vector_record_id != vector.vector_record_id:
        errors.append(ProtocolValidationIssue("vector_id_mismatch", "error", "载体 ID 不一致"))
    if protocol.vector_checksum != vector.normalized_circular_sha256:
        errors.append(
            ProtocolValidationIssue(
                "vector_checksum_mismatch",
                "error",
                "载体序列校验值不一致",
            ),
        )
    if protocol.workflow_type != "promoter_luciferase_reporter":
        errors.append(ProtocolValidationIssue("workflow_mismatch", "error", "工作流类型不一致"))
    if protocol.status != "enabled":
        errors.append(ProtocolValidationIssue("protocol_not_enabled", "error", "protocol 尚未启用"))
    if protocol.right_boundary > len(vector.sequence):
        errors.append(ProtocolValidationIssue("boundary_out_of_range", "error", "插入边界超出载体"))
    else:
        left_start = protocol.left_boundary - len(protocol.left_primer_homology)
        if (
            left_start < 0
            or vector.sequence[left_start : protocol.left_boundary]
            != protocol.left_primer_homology
        ):
            errors.append(
                ProtocolValidationIssue(
                    "left_homology_mismatch",
                    "error",
                    "F 引物载体同源序列与插入边界不一致",
                ),
            )
        right_top = reverse_complement(protocol.right_primer_homology)
        if (
            vector.sequence[
                protocol.right_boundary : protocol.right_boundary + len(right_top)
            ]
            != right_top
        ):
            errors.append(
                ProtocolValidationIssue(
                    "right_homology_mismatch",
                    "error",
                    "R 引物载体同源序列与插入边界不一致",
                ),
            )
    if protocol.experimental_validation_status != "verified":
        warnings.append(
            ProtocolValidationIssue(
                "protocol_unverified",
                "warning",
                "该 reporter protocol 尚未绑定湿实验验证项目",
            ),
        )
    return ReporterVectorProtocolValidationResult(tuple(errors), tuple(warnings))


def _anneal_score(sequence: str) -> tuple[float, int, int]:
    gc = (sequence.count("G") + sequence.count("C")) / len(sequence) * 100
    gc_penalty = max(0.0, 40.0 - gc, gc - 60.0)
    tm = 2 * (sequence.count("A") + sequence.count("T")) + 4 * (
        sequence.count("G") + sequence.count("C")
    )
    return gc_penalty, abs(tm - 65), abs(len(sequence) - 22)


def _choose_forward(sequence: str, protocol: ReporterVectorProtocol) -> str:
    candidates = tuple(
        sequence[:length]
        for length in range(
            protocol.anneal_min_bp,
            min(protocol.anneal_max_bp, len(sequence)) + 1,
        )
    )
    if not candidates:
        raise ValueError("promoter insert 太短，无法生成 F 引物")
    return min(candidates, key=_anneal_score)


def _choose_reverse(sequence: str, protocol: ReporterVectorProtocol) -> str:
    candidates = tuple(
        reverse_complement(sequence[-length:])
        for length in range(
            protocol.anneal_min_bp,
            min(protocol.anneal_max_bp, len(sequence)) + 1,
        )
    )
    if not candidates:
        raise ValueError("promoter insert 太短，无法生成 R 引物")
    return min(candidates, key=_anneal_score)


def apply_reporter_protocol(
    design: ReporterDesignVersion,
    vector: VectorRecord,
    protocol: ReporterVectorProtocol,
) -> ReporterVectorDesignResult:
    validation = validate_reporter_protocol(vector, protocol)
    if not validation.is_valid:
        raise ReporterVectorProtocolError(validation)
    if design.protocol_version_id != protocol.protocol_version_id:
        raise ValueError("reporter 设计与载体 protocol 版本不一致")
    plans = []
    for construct in design.constructs:
        forward_anneal = _choose_forward(construct.insert_sequence, protocol)
        reverse_anneal = _choose_reverse(construct.insert_sequence, protocol)
        forward = ReporterPrimerPlan(
            primer_id=f"{construct.construct_id}-primer-F",
            name=f"{design.gene_symbol}-P{construct.retained_promoter_length}-F",
            sequence=protocol.left_primer_homology + forward_anneal,
            direction="F",
            anneal_length=len(forward_anneal),
        )
        reverse = ReporterPrimerPlan(
            primer_id=f"{construct.construct_id}-primer-R",
            name=f"{design.gene_symbol}-P-R",
            sequence=protocol.right_primer_homology + reverse_anneal,
            direction="R",
            anneal_length=len(reverse_anneal),
        )
        expected = (
            vector.sequence[: protocol.left_boundary]
            + construct.insert_sequence
            + vector.sequence[protocol.right_boundary :]
        )
        plans.append(
            ReporterConstructVectorPlan(
                construct_id=construct.construct_id,
                construct_name=construct.construct_name,
                forward_primer=forward,
                reverse_primer=reverse,
                expected_plasmid_sequence=expected,
                expected_plasmid_checksum=normalized_circular_checksum(expected),
            ),
        )
    return ReporterVectorDesignResult(
        design_version_id=design.design_version_id,
        vector_record_id=vector.vector_record_id,
        vector_checksum=vector.normalized_circular_sha256,
        protocol_version_id=protocol.protocol_version_id,
        construct_plans=tuple(plans),
    )
