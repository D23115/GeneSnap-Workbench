"""Orchestrate the complete deterministic SYN design pipeline."""

from __future__ import annotations

from dataclasses import replace
from datetime import datetime

from genesnap_workbench.domain.syn import (
    SYNDesignInput,
    SYNDesignVersion,
    SYNManualOverrideRecord,
)
from genesnap_workbench.vector_library.models import SYNVectorProtocol, VectorRecord
from genesnap_workbench.vector_library.syn import simulate_syn_plasmid

from .dna import normalize_dna, sha256_sequence
from .syn_modules import SYNModuleRules, plan_syn_modules
from .syn_oligos import (
    BiopythonThermodynamicAnalyzer,
    SYNOligoRules,
    design_assembly_oligos,
)
from .syn_qc import SYNQCRules, evaluate_syn_sequence


def create_syn_design(
    design_input: SYNDesignInput,
    vector: VectorRecord,
    protocol: SYNVectorProtocol,
    *,
    design_version_id: str,
    created_at: datetime,
    qc_rules: SYNQCRules = SYNQCRules(),
    module_rules: SYNModuleRules = SYNModuleRules(),
    oligo_rules: SYNOligoRules = SYNOligoRules(),
    manual_overrides: tuple[SYNManualOverrideRecord, ...] = (),
) -> SYNDesignVersion:
    """Create one immutable design version from input through plasmid simulation."""
    if design_input.vector_protocol_version_id != protocol.protocol_version_id:
        raise ValueError("输入绑定的 vector protocol 与当前 protocol 不一致")
    normalized = normalize_dna(design_input.raw_sequence)
    normalized_checksum = sha256_sequence(normalized)
    qc_result = evaluate_syn_sequence(
        normalized,
        qc_rules,
        design_version_id=design_version_id,
    )
    module_plan = plan_syn_modules(
        normalized,
        module_rules,
        design_version_id=design_version_id,
        qc_result=qc_result,
    )
    oligo_result = design_assembly_oligos(
        normalized,
        module_plan,
        oligo_rules,
        BiopythonThermodynamicAnalyzer(),
        design_version_id=design_version_id,
        project_id=design_input.project_id,
        target_name=design_input.target_name,
    )
    plasmid = simulate_syn_plasmid(
        vector,
        protocol,
        normalized,
        design_version_id=design_version_id,
    )
    warnings = tuple(
        dict.fromkeys(
            (*qc_result.confirmable_warnings, *oligo_result.warnings),
        ),
    )
    return SYNDesignVersion(
        design_version_id=design_version_id,
        project_id=design_input.project_id,
        version_no=1,
        created_at=created_at,
        raw_sequence_checksum=normalized_checksum,
        normalized_sequence=normalized,
        normalized_checksum=normalized_checksum,
        final_sequence=normalized,
        final_checksum=normalized_checksum,
        qc_result=qc_result,
        module_plan=oligo_result.module_plan,
        oligos=oligo_result.oligos,
        plasmid_simulation=plasmid,
        rule_versions=(
            qc_rules.rules_version,
            "syn-module-v1",
            "syn-oligo-v1",
            protocol.protocol_version_id,
        ),
        manual_overrides=manual_overrides,
        design_warnings=warnings,
        requires_confirmation=(
            bool(qc_result.confirmable_warnings)
            or oligo_result.requires_confirmation
        ),
    )


def confirm_syn_design_warnings(
    design: SYNDesignVersion,
    *,
    override_id: str,
    reason: str,
    actor: str,
    occurred_at: datetime,
) -> SYNDesignVersion:
    """Append an auditable confirmation without erasing the detected warnings."""
    if not reason.strip():
        raise ValueError("design confirmation reason must not be blank")
    if any(
        record.override_id == override_id for record in design.manual_overrides
    ):
        raise ValueError(f"Manual override already exists: {override_id}")
    confirmation = SYNManualOverrideRecord(
        override_id=override_id,
        field_path="design_confirmation",
        old_value="required",
        new_value="confirmed",
        reason=reason.strip(),
        occurred_at=occurred_at,
        actor=actor,
    )
    return replace(
        design,
        manual_overrides=design.manual_overrides + (confirmation,),
    )
