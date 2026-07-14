"""Workflow registration contracts shared by all GeneSnap workflows."""

from __future__ import annotations

from dataclasses import dataclass
from typing import Iterable


@dataclass(frozen=True, slots=True)
class WorkflowTransition:
    """Declare one controlled transition in a workflow state graph."""

    from_state: str
    action: str
    to_state: str


@dataclass(frozen=True, slots=True)
class WorkflowDefinition:
    """Describe one workflow without embedding executable protocol logic."""

    workflow_type: str
    project_category: str
    folder_suffix: str
    enabled: bool
    intake_schema_key: str
    intake_field_keys: tuple[str, ...]
    protocol_types: tuple[str, ...]
    validator_key: str
    design_engine_adapter_key: str
    state_graph: tuple[WorkflowTransition, ...]
    artifact_types: tuple[str, ...]


class WorkflowRegistry:
    """Store workflow definitions and expose only enabled entries to selectors."""

    def __init__(self, definitions: Iterable[WorkflowDefinition] = ()) -> None:
        self._definitions: dict[str, WorkflowDefinition] = {}
        for definition in definitions:
            self.register(definition)

    def register(self, definition: WorkflowDefinition) -> None:
        """Register one unique workflow definition."""
        workflow_type = definition.workflow_type
        if workflow_type in self._definitions:
            raise ValueError(f"Workflow type is already registered: {workflow_type}")
        self._definitions[workflow_type] = definition

    def get(self, workflow_type: str) -> WorkflowDefinition:
        """Return a registered workflow, including disabled definitions."""
        return self._definitions[workflow_type]

    def list_enabled(self) -> tuple[WorkflowDefinition, ...]:
        """Return enabled workflows in deterministic registration order."""
        return tuple(
            definition
            for definition in self._definitions.values()
            if definition.enabled
        )


def build_default_workflow_registry() -> WorkflowRegistry:
    """Build the first-party registry for workflows implemented in this build."""
    syn = WorkflowDefinition(
        workflow_type="de_novo_gene_synthesis",
        project_category="合成/组装类",
        folder_suffix="SYN",
        enabled=True,
        intake_schema_key="syn_v0_intake",
        intake_field_keys=(
            "project_id",
            "target_name",
            "raw_sequence",
            "input_format",
            "gene_symbol",
            "species",
            "vector_protocol_version_id",
            "received_date",
            "standard_completion_date",
        ),
        protocol_types=(
            "pUC57-EcoRV",
            "pUC57-SmaI",
            "custom_syn_vector",
        ),
        validator_key="syn_v0_validator",
        design_engine_adapter_key="syn_v0_design_engine",
        state_graph=(
            WorkflowTransition("recorded", "complete_design", "design_completed"),
            WorkflowTransition(
                "design_completed",
                "mark_materials_ordered",
                "materials_ordered",
            ),
            WorkflowTransition(
                "materials_ordered",
                "mark_materials_arrived",
                "materials_arrived",
            ),
            WorkflowTransition(
                "materials_arrived",
                "start_assembly",
                "syn_assembly_in_progress",
            ),
            WorkflowTransition(
                "syn_assembly_in_progress",
                "select_clones_for_prep",
                "plasmid_prep_in_progress",
            ),
            WorkflowTransition(
                "plasmid_prep_in_progress",
                "finish_plasmid_prep",
                "awaiting_sequencing_confirmation",
            ),
            WorkflowTransition(
                "awaiting_sequencing_confirmation",
                "additional_screening",
                "syn_assembly_in_progress",
            ),
            WorkflowTransition(
                "awaiting_sequencing_confirmation",
                "restart_assembly",
                "syn_assembly_in_progress",
            ),
            WorkflowTransition(
                "awaiting_sequencing_confirmation",
                "complete_project",
                "project_completed",
            ),
        ),
        artifact_types=(
            "design_json",
            "oligo_order_xlsx",
            "oligo_resuspension_xlsx",
            "oligo_mix_xlsx",
            "qc_report_docx",
            "assembly_plan_docx",
            "colony_pcr_xlsx",
            "expected_plasmid_genbank",
        ),
    )
    shrna = WorkflowDefinition(
        workflow_type="shrna_knockdown",
        project_category="沉默/敲低类",
        folder_suffix="KD",
        enabled=True,
        intake_schema_key="shrna_v1_intake",
        intake_field_keys=(
            "project_id",
            "gene_symbol",
            "species",
            "cds_sequence",
            "transcript_accession",
            "target_count",
            "clones_per_target",
            "vector_protocol_version_id",
            "received_date",
            "standard_completion_date",
        ),
        protocol_types=("pLKO.1-AgeI-EcoRI", "custom_shrna_vector"),
        validator_key="shrna_v1_validator",
        design_engine_adapter_key="shrna_v1_design_engine",
        state_graph=(
            WorkflowTransition("recorded", "complete_design", "design_completed"),
            WorkflowTransition("design_completed", "mark_primers_ordered", "primers_ordered"),
            WorkflowTransition("primers_ordered", "mark_primers_arrived", "primers_arrived"),
            WorkflowTransition("primers_arrived", "start_cloning", "cloning_in_progress"),
            WorkflowTransition("cloning_in_progress", "mark_sent", "sequencing_in_progress"),
            WorkflowTransition(
                "sequencing_in_progress",
                "analyze_sequencing",
                "analysis_completed",
            ),
            WorkflowTransition(
                "analysis_completed",
                "start_plasmid_prep",
                "plasmid_prep_in_progress",
            ),
            WorkflowTransition(
                "plasmid_prep_in_progress",
                "complete_project",
                "project_completed",
            ),
        ),
        artifact_types=(
            "design_json",
            "primer_order_xlsx",
            "sequencing_order_xlsx",
            "design_report_docx",
            "expected_plasmid_genbank",
            "sequencing_analysis_report",
        ),
    )
    return WorkflowRegistry((syn, shrna))
