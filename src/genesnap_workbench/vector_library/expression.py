"""通用表达设计与数据化载体 protocol 的连接层。"""

from __future__ import annotations

from genesnap_workbench.domain.expression import (
    ExpressionPrimerPlan,
    ExpressionConstructVectorPlan,
    ExpressionDesignVersion,
    ExpressionVectorDesignResult,
    ExpressionVectorRules,
)
from genesnap_workbench.sequence_core.dna import reverse_complement

from .models import (
    ExpressionVectorProtocol,
    ExpressionVectorProtocolValidationResult,
    ProtocolValidationIssue,
    VectorRecord,
    normalized_circular_checksum,
)


class ExpressionVectorProtocolError(ValueError):
    def __init__(self, validation: ExpressionVectorProtocolValidationResult) -> None:
        self.validation = validation
        super().__init__("; ".join(item.message for item in validation.errors))


def expression_rules_from_protocol(
    protocol: ExpressionVectorProtocol,
) -> ExpressionVectorRules:
    return ExpressionVectorRules(
        protocol_version_id=protocol.protocol_version_id,
        kozak_sequence=protocol.kozak_sequence,
        stop_codon_rule=protocol.stop_codon_rule,
        c_terminal_fusion_name=protocol.c_terminal_fusion_name,
        single_fragment_max_bp=protocol.single_fragment_max_bp,
    )


def validate_expression_protocol(
    vector: VectorRecord,
    protocol: ExpressionVectorProtocol,
) -> ExpressionVectorProtocolValidationResult:
    errors: list[ProtocolValidationIssue] = []
    warnings: list[ProtocolValidationIssue] = []
    if protocol.vector_record_id != vector.vector_record_id:
        errors.append(
            ProtocolValidationIssue("vector_id_mismatch", "error", "protocol 载体 ID 不一致"),
        )
    if protocol.vector_checksum != vector.normalized_circular_sha256:
        errors.append(
            ProtocolValidationIssue(
                "vector_checksum_mismatch",
                "error",
                "protocol 载体序列校验值不一致",
            ),
        )
    if protocol.workflow_type != "expression":
        errors.append(
            ProtocolValidationIssue("workflow_mismatch", "error", "protocol 不属于表达类"),
        )
    if protocol.status != "enabled":
        errors.append(
            ProtocolValidationIssue("protocol_not_enabled", "error", "protocol 尚未启用"),
        )
    if protocol.right_boundary > len(vector.sequence):
        errors.append(
            ProtocolValidationIssue("boundary_out_of_range", "error", "插入边界超出载体长度"),
        )
    else:
        left_start = protocol.left_boundary - len(protocol.left_primer_homology)
        left_context = vector.sequence[max(0, left_start) : protocol.left_boundary]
        if left_start < 0 or left_context != protocol.left_primer_homology:
            errors.append(
                ProtocolValidationIssue(
                    "left_homology_mismatch",
                    "error",
                    "左侧引物同源序列与载体边界不一致",
                ),
            )
        right_top = reverse_complement(protocol.right_primer_homology)
        right_context = vector.sequence[
            protocol.right_boundary : protocol.right_boundary + len(right_top)
        ]
        if right_context != right_top:
            errors.append(
                ProtocolValidationIssue(
                    "right_homology_mismatch",
                    "error",
                    "右侧引物同源序列与载体边界不一致",
                ),
            )
    if protocol.experimental_validation_status != "verified":
        warnings.append(
            ProtocolValidationIssue(
                "protocol_unverified",
                "warning",
                "该表达 protocol 尚未绑定湿实验验证项目",
            ),
        )
    return ExpressionVectorProtocolValidationResult(tuple(errors), tuple(warnings))


def _gc_percent(sequence: str) -> float:
    return (sequence.count("G") + sequence.count("C")) / len(sequence) * 100


def _wallace_tm(sequence: str) -> int:
    return 2 * (sequence.count("A") + sequence.count("T")) + 4 * (
        sequence.count("G") + sequence.count("C")
    )


def _anneal_score(sequence: str) -> tuple[float, float, int]:
    gc = _gc_percent(sequence)
    gc_penalty = max(0.0, 40.0 - gc, gc - 60.0)
    return gc_penalty, abs(_wallace_tm(sequence) - 65), abs(len(sequence) - 22)


def _choose_forward_anneal(
    sequence: str,
    protocol: ExpressionVectorProtocol,
) -> str:
    if len(sequence) < protocol.anneal_min_bp:
        raise ValueError("insert 太短，无法生成满足最短退火长度的引物")
    candidates = tuple(
        sequence[:length]
        for length in range(
            protocol.anneal_min_bp,
            min(protocol.anneal_max_bp, len(sequence)) + 1,
        )
    )
    return min(candidates, key=_anneal_score)


def _choose_reverse_anneal(
    sequence: str,
    protocol: ExpressionVectorProtocol,
) -> str:
    candidates = tuple(
        reverse_complement(sequence[-length:])
        for length in range(
            protocol.anneal_min_bp,
            min(protocol.anneal_max_bp, len(sequence)) + 1,
        )
    )
    return min(candidates, key=_anneal_score)


def apply_expression_protocol(
    design: ExpressionDesignVersion,
    vector: VectorRecord,
    protocol: ExpressionVectorProtocol,
) -> ExpressionVectorDesignResult:
    validation = validate_expression_protocol(vector, protocol)
    if not validation.is_valid:
        raise ExpressionVectorProtocolError(validation)
    if design.protocol_version_id != protocol.protocol_version_id:
        raise ValueError("表达设计与载体 protocol 版本不一致")
    plans = []
    for construct in design.constructs:
        forward_anneal = _choose_forward_anneal(construct.coding_sequence, protocol)
        reverse_anneal = _choose_reverse_anneal(construct.coding_sequence, protocol)
        forward_primer = (
            protocol.left_primer_homology
            + protocol.kozak_sequence
            + forward_anneal
        )
        reverse_primer = protocol.right_primer_homology + reverse_anneal
        if len(construct.fragments) == 1:
            primers = (
                ExpressionPrimerPlan(
                    name=f"{construct.construct_name}-F",
                    sequence=forward_primer,
                    direction="F",
                    fragment_no=1,
                    anneal_length=len(forward_anneal),
                    overlap_length=len(protocol.left_primer_homology) + len(protocol.kozak_sequence),
                ),
                ExpressionPrimerPlan(
                    name=f"{construct.construct_name}-R",
                    sequence=reverse_primer,
                    direction="R",
                    fragment_no=1,
                    anneal_length=len(reverse_anneal),
                    overlap_length=len(protocol.right_primer_homology),
                ),
            )
        elif len(construct.fragments) == 2:
            split = (len(construct.coding_sequence) // 2 // 3) * 3
            left_template = construct.coding_sequence[:split]
            right_template = construct.coding_sequence[split:]
            left_reverse_anneal = _choose_reverse_anneal(left_template, protocol)
            right_forward_anneal = _choose_forward_anneal(right_template, protocol)
            internal_overlap = right_template[: len(right_forward_anneal)]
            primers = (
                ExpressionPrimerPlan(
                    name=f"{construct.construct_name}-P1-F",
                    sequence=forward_primer,
                    direction="F",
                    fragment_no=1,
                    anneal_length=len(forward_anneal),
                    overlap_length=len(protocol.left_primer_homology) + len(protocol.kozak_sequence),
                ),
                ExpressionPrimerPlan(
                    name=f"{construct.construct_name}-P1-R",
                    sequence=reverse_complement(internal_overlap) + left_reverse_anneal,
                    direction="R",
                    fragment_no=1,
                    anneal_length=len(left_reverse_anneal),
                    overlap_length=len(internal_overlap),
                ),
                ExpressionPrimerPlan(
                    name=f"{construct.construct_name}-P2-F",
                    sequence=right_forward_anneal,
                    direction="F",
                    fragment_no=2,
                    anneal_length=len(right_forward_anneal),
                    overlap_length=len(right_forward_anneal),
                ),
                ExpressionPrimerPlan(
                    name=f"{construct.construct_name}-P2-R",
                    sequence=reverse_primer,
                    direction="R",
                    fragment_no=2,
                    anneal_length=len(reverse_anneal),
                    overlap_length=len(protocol.right_primer_homology),
                ),
            )
        else:
            raise ValueError("MVP expression protocol supports one or two PCR fragments")
        expected = (
            vector.sequence[: protocol.left_boundary]
            + construct.insert_sequence
            + vector.sequence[protocol.right_boundary :]
        )
        plans.append(
            ExpressionConstructVectorPlan(
                construct_id=construct.construct_id,
                construct_name=construct.construct_name,
                forward_primer=forward_primer,
                reverse_primer=reverse_primer,
                forward_anneal_length=len(forward_anneal),
                reverse_anneal_length=len(reverse_anneal),
                expected_plasmid_sequence=expected,
                expected_plasmid_checksum=normalized_circular_checksum(expected),
                left_junction=(
                    vector.sequence[protocol.left_boundary - 12 : protocol.left_boundary]
                    + construct.insert_sequence[:24]
                ),
                right_junction=(
                    construct.insert_sequence[-24:]
                    + vector.sequence[protocol.right_boundary : protocol.right_boundary + 24]
                ),
                primers=primers,
            ),
        )
    return ExpressionVectorDesignResult(
        design_version_id=design.design_version_id,
        vector_record_id=vector.vector_record_id,
        vector_checksum=vector.normalized_circular_sha256,
        protocol_version_id=protocol.protocol_version_id,
        construct_plans=tuple(plans),
    )
