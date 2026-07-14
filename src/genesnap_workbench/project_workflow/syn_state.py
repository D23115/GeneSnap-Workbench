"""Controlled status and assembly-substep transitions for SYN projects."""

from __future__ import annotations

from dataclasses import replace
from datetime import datetime

from genesnap_workbench.domain.syn import (
    SYNAuditEvent,
    SYNManualOverrideRecord,
    SYNProjectSnapshot,
)
from genesnap_workbench.domain.workflows import (
    WorkflowDefinition,
    build_default_workflow_registry,
)

from .syn_materials import MaterialReadiness


SYN_WORKFLOW_TYPE = "de_novo_gene_synthesis"
ASSEMBLY_SUBSTEPS = (
    "assembly_pcr",
    "amplification_pcr",
    "vector_assembly_transformation",
    "colony_pcr",
)


class SYNStateTransitionError(ValueError):
    pass


class SYNMaterialOverrideRequired(SYNStateTransitionError):
    pass


def display_status_label(status: str, workflow_type: str) -> str:
    """Return workflow-aware Chinese labels for shared internal statuses."""
    if workflow_type == SYN_WORKFLOW_TYPE:
        syn_labels = {
            "design_completed": "设计完成/待订购",
            "materials_ordered": "oligo 已订购",
            "materials_arrived": "oligo 已到货",
            "syn_assembly_in_progress": "合成组装中",
            "plasmid_prep_in_progress": "质粒抽提中",
            "awaiting_sequencing_confirmation": "待测序确认",
            "project_completed": "项目完成",
        }
        if status in syn_labels:
            return syn_labels[status]
    shared_labels = {
        "recorded": "已录入",
        "design_completed": "设计完成/待订购",
        "materials_ordered": "引物已订购",
        "materials_arrived": "引物已到货",
        "project_completed": "项目完成",
        "abnormal_or_paused": "异常/暂停",
    }
    return shared_labels.get(status, status)


class SYNStateTransitionService:
    def __init__(self, definition: WorkflowDefinition | None = None) -> None:
        self.definition = definition or build_default_workflow_registry().get(
            SYN_WORKFLOW_TYPE,
        )

    def _require_transition(
        self,
        snapshot: SYNProjectSnapshot,
        action: str,
        to_status: str,
    ) -> None:
        allowed = any(
            transition.from_state == snapshot.status
            and transition.action == action
            and transition.to_state == to_status
            for transition in self.definition.state_graph
        )
        if not allowed:
            raise SYNStateTransitionError(
                f"当前状态 {snapshot.status} 不能执行 {action} -> {to_status}",
            )

    def _event(
        self,
        *,
        event_id: str,
        event_type: str,
        occurred_at: datetime,
        actor: str,
        from_status: str,
        to_status: str,
        note: str | None = None,
    ) -> SYNAuditEvent:
        return SYNAuditEvent(
            event_id=event_id,
            event_type=event_type,
            occurred_at=occurred_at,
            actor=actor,
            note=note,
            from_status=from_status,
            to_status=to_status,
            source="user_action",
        )

    def mark_materials_arrived(
        self,
        snapshot: SYNProjectSnapshot,
        *,
        resuspension_complete: bool,
        event_id: str,
        actor: str,
        occurred_at: datetime,
    ) -> SYNProjectSnapshot:
        to_status = "materials_arrived"
        self._require_transition(
            snapshot,
            "mark_materials_arrived",
            to_status,
        )
        updated = replace(
            snapshot,
            status=to_status,
            resuspension_data_status=(
                "complete" if resuspension_complete else "missing"
            ),
        )
        return updated.append_status_event(
            self._event(
                event_id=event_id,
                event_type="mark_materials_arrived",
                occurred_at=occurred_at,
                actor=actor,
                from_status=snapshot.status,
                to_status=to_status,
            ),
        )

    def mark_materials_ordered(
        self,
        snapshot: SYNProjectSnapshot,
        *,
        event_id: str,
        actor: str,
        occurred_at: datetime,
    ) -> SYNProjectSnapshot:
        to_status = "materials_ordered"
        self._require_transition(snapshot, "mark_materials_ordered", to_status)
        updated = replace(snapshot, status=to_status)
        return updated.append_status_event(
            self._event(
                event_id=event_id,
                event_type="mark_materials_ordered",
                occurred_at=occurred_at,
                actor=actor,
                from_status=snapshot.status,
                to_status=to_status,
            ),
        )

    def start_assembly(
        self,
        snapshot: SYNProjectSnapshot,
        readiness: MaterialReadiness,
        *,
        confirm_missing: bool,
        override_reason: str | None,
        event_id: str,
        override_id: str,
        actor: str,
        occurred_at: datetime,
    ) -> SYNProjectSnapshot:
        to_status = "syn_assembly_in_progress"
        self._require_transition(snapshot, "start_assembly", to_status)
        if readiness.errors or not readiness.can_start_with_override:
            details = "；".join(readiness.errors) or "材料数据存在阻断错误"
            raise SYNStateTransitionError(details)
        if not readiness.is_ready and not confirm_missing:
            raise SYNMaterialOverrideRequired("材料数据不完整，需要确认后继续")
        if not readiness.is_ready and not (override_reason or "").strip():
            raise ValueError("override reason must not be blank")

        override_history = snapshot.manual_override_history
        if not readiness.is_ready:
            missing = (
                f"oligos={','.join(readiness.missing_oligo_ids)};"
                f"fields={','.join(readiness.missing_fields)}"
            )
            override_history += (
                SYNManualOverrideRecord(
                    override_id=override_id,
                    field_path="materials_readiness",
                    old_value=missing,
                    new_value="confirmed_start_assembly",
                    reason=override_reason.strip(),
                    occurred_at=occurred_at,
                    actor=actor,
                ),
            )

        updated = replace(
            snapshot,
            status=to_status,
            syn_assembly_substep=ASSEMBLY_SUBSTEPS[0],
            manual_override_history=override_history,
        )
        return updated.append_status_event(
            self._event(
                event_id=event_id,
                event_type="start_assembly",
                occurred_at=occurred_at,
                actor=actor,
                from_status=snapshot.status,
                to_status=to_status,
                note=override_reason if not readiness.is_ready else None,
            ),
        )

    def advance_assembly_substep(
        self,
        snapshot: SYNProjectSnapshot,
        *,
        to_substep: str,
        event_id: str,
        actor: str,
        occurred_at: datetime,
    ) -> SYNProjectSnapshot:
        if snapshot.status != "syn_assembly_in_progress":
            raise SYNStateTransitionError("项目当前不在合成组装中")
        try:
            current_index = ASSEMBLY_SUBSTEPS.index(snapshot.syn_assembly_substep)
        except ValueError:
            raise SYNStateTransitionError("当前组装子步骤无效") from None
        expected = (
            ASSEMBLY_SUBSTEPS[current_index + 1]
            if current_index + 1 < len(ASSEMBLY_SUBSTEPS)
            else None
        )
        if to_substep != expected:
            raise SYNStateTransitionError(
                f"只能从 {snapshot.syn_assembly_substep} 推进到 {expected}",
            )
        updated = replace(snapshot, syn_assembly_substep=to_substep)
        return updated.append_status_event(
            self._event(
                event_id=event_id,
                event_type="advance_assembly_substep",
                occurred_at=occurred_at,
                actor=actor,
                from_status=snapshot.status,
                to_status=snapshot.status,
                note=f"{snapshot.syn_assembly_substep} -> {to_substep}",
            ),
        )

    def correct_status(
        self,
        snapshot: SYNProjectSnapshot,
        *,
        to_status: str,
        reason: str,
        event_id: str,
        actor: str,
        occurred_at: datetime,
    ) -> SYNProjectSnapshot:
        if not reason.strip():
            raise ValueError("status correction reason must not be blank")
        known_states = {
            state
            for transition in self.definition.state_graph
            for state in (transition.from_state, transition.to_state)
        }
        if to_status not in known_states:
            raise SYNStateTransitionError(f"未知目标状态: {to_status}")
        updated = replace(
            snapshot,
            status=to_status,
            syn_assembly_substep=(
                snapshot.syn_assembly_substep
                if to_status == "syn_assembly_in_progress"
                else None
            ),
        )
        return updated.append_status_event(
            self._event(
                event_id=event_id,
                event_type="correct_status",
                occurred_at=occurred_at,
                actor=actor,
                from_status=snapshot.status,
                to_status=to_status,
                note=reason.strip(),
            ),
        )
