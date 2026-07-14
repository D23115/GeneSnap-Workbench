"""Use-case service for SYN assembly rounds, colonies, and plasmid prep."""

from __future__ import annotations

from dataclasses import dataclass, replace
from datetime import datetime

from genesnap_workbench.domain.syn import (
    PlasmidPrepRecord,
    PlasmidPrepStatus,
    SYNAuditEvent,
    SYNAssemblyAttemptRecord,
    SYNAssemblyAttemptResult,
    SYNAssemblyStepRecord,
    SYNColonyPCRRecord,
    SYNColonyPCRResult,
    SYNProjectSnapshot,
    SYNRoute,
    SYNSequencingConfirmation,
    SYNSequencingResult,
)

from .syn_state import ASSEMBLY_SUBSTEPS


class SYNWorkflowRuleError(ValueError):
    """Raised when a requested action violates the SYN workflow rules."""


class SYNRevisionConflict(SYNWorkflowRuleError):
    """Raised when the caller is acting on a stale project snapshot."""


@dataclass(frozen=True, slots=True)
class AdditionalScreeningPreview:
    preview_id: str
    project_id: str
    expected_revision: int
    attempt_id: str
    target_name: str
    clone_numbers: tuple[int, ...]
    clone_ids: tuple[str, ...]
    display_names: tuple[str, ...]


@dataclass(frozen=True, slots=True)
class SYNSequencingSummary:
    correct_count: int
    incorrect_count: int
    uncertain_count: int
    unconfirmed_count: int
    usable_clone_ids: tuple[str, ...]
    usable_clone_names: tuple[str, ...]

    @property
    def display_text(self) -> str:
        return (
            f"人工确认：正确 {self.correct_count} / "
            f"错误 {self.incorrect_count} / "
            f"不确定 {self.uncertain_count} / "
            f"未确认 {self.unconfirmed_count}"
        )


def _latest_by_clone(records):
    latest = {}
    for record in records:
        latest[record.clone_id] = record
    return latest


class SYNWorkflowService:
    """Apply controlled, append-only actions to an immutable SYN snapshot."""

    def _check_revision(
        self,
        snapshot: SYNProjectSnapshot,
        expected_revision: int,
    ) -> None:
        if snapshot.revision != expected_revision:
            raise SYNRevisionConflict(
                f"项目版本已变化：预期 {expected_revision}，实际 {snapshot.revision}",
            )

    def _require_assembly_colony_step(
        self,
        snapshot: SYNProjectSnapshot,
    ) -> None:
        if (
            snapshot.status != "syn_assembly_in_progress"
            or snapshot.syn_assembly_substep != "colony_pcr"
        ):
            raise SYNWorkflowRuleError("当前项目不在菌落 PCR 步骤")

    def _append_event(
        self,
        snapshot: SYNProjectSnapshot,
        *,
        event_id: str,
        event_type: str,
        actor: str,
        occurred_at: datetime,
        from_status: str,
        note: str | None = None,
        related_entity_id: str | None = None,
    ) -> SYNProjectSnapshot:
        return snapshot.append_status_event(
            SYNAuditEvent(
                event_id=event_id,
                event_type=event_type,
                occurred_at=occurred_at,
                actor=actor,
                note=note,
                from_status=from_status,
                to_status=snapshot.status,
                source="user_action",
                related_entity_id=related_entity_id,
            ),
        )

    def _make_colonies(
        self,
        *,
        attempt_id: str,
        round_no: int,
        target_name: str,
        clone_numbers: tuple[int, ...],
    ) -> tuple[SYNColonyPCRRecord, ...]:
        if not target_name.strip():
            raise ValueError("target_name must not be blank")
        prefix = target_name.strip()
        if round_no > 1:
            prefix = f"{prefix}-R{round_no}"
        return tuple(
            SYNColonyPCRRecord(
                clone_id=f"{attempt_id}-clone-{clone_no}",
                attempt_id=attempt_id,
                clone_no=clone_no,
                display_name=f"{prefix}-C{clone_no:02d}",
                result=SYNColonyPCRResult.PENDING,
            )
            for clone_no in clone_numbers
        )

    def start_initial_colony_screening(
        self,
        snapshot: SYNProjectSnapshot,
        *,
        target_name: str,
        route: SYNRoute,
        high_risk: bool,
        colony_count: int | None,
        expected_revision: int,
        attempt_id: str,
        event_id: str,
        actor: str,
        occurred_at: datetime,
    ) -> SYNProjectSnapshot:
        self._check_revision(snapshot, expected_revision)
        self._require_assembly_colony_step(snapshot)
        if not isinstance(route, SYNRoute):
            raise TypeError("route must be a SYNRoute")
        count = colony_count
        if count is None:
            count = 12 if high_risk or route is SYNRoute.MODULAR else 8
        if count <= 0:
            raise ValueError("colony_count must be positive")

        if snapshot.attempts:
            attempt = self._active_attempt(snapshot)
            if attempt.syn_assembly_round_no != 1:
                raise SYNWorkflowRuleError("初始菌落筛选只能用于第一轮组装")
            updated = snapshot
        else:
            attempt = SYNAssemblyAttemptRecord(
                attempt_id=attempt_id,
                project_id=snapshot.project_id,
                design_version_id=snapshot.active_design_version_id,
                syn_assembly_round_no=1,
                restart_from_substep="assembly_pcr",
                result=SYNAssemblyAttemptResult.PENDING,
                started_at=occurred_at,
            )
            updated = snapshot.append_attempt(attempt)
        colonies = self._make_colonies(
            attempt_id=attempt.attempt_id,
            round_no=1,
            target_name=target_name,
            clone_numbers=tuple(range(1, count + 1)),
        )
        updated = replace(updated, colonies=updated.colonies + colonies)
        return self._append_event(
            updated,
            event_id=event_id,
            event_type="start_initial_colony_screening",
            actor=actor,
            occurred_at=occurred_at,
            from_status=snapshot.status,
            note=f"创建 {count} 个菌落记录",
            related_entity_id=attempt.attempt_id,
        )

    def begin_initial_attempt(
        self,
        snapshot: SYNProjectSnapshot,
        *,
        expected_revision: int,
        attempt_id: str,
        event_id: str,
        actor: str,
        occurred_at: datetime,
    ) -> SYNProjectSnapshot:
        self._check_revision(snapshot, expected_revision)
        if snapshot.status != "syn_assembly_in_progress":
            raise SYNWorkflowRuleError("项目当前不在合成组装中")
        if snapshot.syn_assembly_substep != "assembly_pcr":
            raise SYNWorkflowRuleError("首轮实验必须从 Assembly PCR 开始")
        if snapshot.attempts or snapshot.syn_assembly_round_no != 0:
            raise SYNWorkflowRuleError("首轮实验记录已经存在")
        attempt = SYNAssemblyAttemptRecord(
            attempt_id=attempt_id,
            project_id=snapshot.project_id,
            design_version_id=snapshot.active_design_version_id,
            syn_assembly_round_no=1,
            restart_from_substep="assembly_pcr",
            result=SYNAssemblyAttemptResult.PENDING,
            started_at=occurred_at,
        )
        updated = snapshot.append_attempt(attempt)
        return self._append_event(
            updated,
            event_id=event_id,
            event_type="begin_initial_assembly_attempt",
            actor=actor,
            occurred_at=occurred_at,
            from_status=snapshot.status,
            note="开始第 1 轮合成组装",
            related_entity_id=attempt_id,
        )

    def create_colonies_for_active_round(
        self,
        snapshot: SYNProjectSnapshot,
        *,
        target_name: str,
        colony_count: int,
        expected_revision: int,
        event_id: str,
        actor: str,
        occurred_at: datetime,
    ) -> SYNProjectSnapshot:
        self._check_revision(snapshot, expected_revision)
        self._require_assembly_colony_step(snapshot)
        if colony_count <= 0:
            raise ValueError("colony_count must be positive")
        active_attempt = self._active_attempt(snapshot)
        active_records = (
            record
            for record in snapshot.colonies
            if record.attempt_id == active_attempt.attempt_id
        )
        max_clone_no = max((item.clone_no for item in active_records), default=0)
        clone_numbers = tuple(
            range(max_clone_no + 1, max_clone_no + colony_count + 1),
        )
        colonies = self._make_colonies(
            attempt_id=active_attempt.attempt_id,
            round_no=active_attempt.syn_assembly_round_no,
            target_name=target_name,
            clone_numbers=clone_numbers,
        )
        updated = replace(
            snapshot,
            syn_assembly_substep="colony_pcr",
            colonies=snapshot.colonies + colonies,
        )
        return self._append_event(
            updated,
            event_id=event_id,
            event_type="create_colonies_for_active_round",
            actor=actor,
            occurred_at=occurred_at,
            from_status=snapshot.status,
            note=f"新增菌落 {clone_numbers[0]}-{clone_numbers[-1]}",
            related_entity_id=active_attempt.attempt_id,
        )

    def record_colony_pcr(
        self,
        snapshot: SYNProjectSnapshot,
        *,
        clone_id: str,
        result: SYNColonyPCRResult,
        observed_note: str | None,
        expected_revision: int,
        event_id: str,
        actor: str,
        occurred_at: datetime,
    ) -> SYNProjectSnapshot:
        self._check_revision(snapshot, expected_revision)
        self._require_assembly_colony_step(snapshot)
        if not isinstance(result, SYNColonyPCRResult):
            raise TypeError("result must be a SYNColonyPCRResult")
        latest = _latest_by_clone(snapshot.colonies)
        if clone_id not in latest:
            raise SYNWorkflowRuleError(f"找不到菌落克隆：{clone_id}")
        source = latest[clone_id]
        recorded = replace(
            source,
            result=result,
            observed_note=observed_note,
            recorded_at=occurred_at,
        )
        updated = replace(snapshot, colonies=snapshot.colonies + (recorded,))
        return self._append_event(
            updated,
            event_id=event_id,
            event_type="record_colony_pcr",
            actor=actor,
            occurred_at=occurred_at,
            from_status=snapshot.status,
            note=result.value,
            related_entity_id=clone_id,
        )

    def record_assembly_step(
        self,
        snapshot: SYNProjectSnapshot,
        *,
        result: SYNAssemblyAttemptResult,
        note: str | None,
        expected_revision: int,
        record_id: str,
        event_id: str,
        actor: str,
        occurred_at: datetime,
    ) -> SYNProjectSnapshot:
        self._check_revision(snapshot, expected_revision)
        if snapshot.status != "syn_assembly_in_progress":
            raise SYNWorkflowRuleError("项目当前不在合成组装中")
        if snapshot.syn_assembly_substep not in ASSEMBLY_SUBSTEPS:
            raise SYNWorkflowRuleError("当前合成组装子步骤无效")
        if not isinstance(result, SYNAssemblyAttemptResult):
            raise TypeError("result must be a SYNAssemblyAttemptResult")
        attempt = self._active_attempt(snapshot)
        step_attempt_no = 1 + sum(
            1
            for record in snapshot.step_records
            if record.attempt_id == attempt.attempt_id
            and record.substep == snapshot.syn_assembly_substep
        )
        record = SYNAssemblyStepRecord(
            record_id=record_id,
            attempt_id=attempt.attempt_id,
            substep=snapshot.syn_assembly_substep,
            step_attempt_no=step_attempt_no,
            result=result,
            recorded_at=occurred_at,
            actor=actor,
            note=note,
        )
        updated = replace(
            snapshot,
            step_records=snapshot.step_records + (record,),
        )
        return self._append_event(
            updated,
            event_id=event_id,
            event_type="record_assembly_step",
            actor=actor,
            occurred_at=occurred_at,
            from_status=snapshot.status,
            note=f"{record.substep} 第 {step_attempt_no} 次：{result.value}",
            related_entity_id=record_id,
        )

    def select_clones_for_prep(
        self,
        snapshot: SYNProjectSnapshot,
        *,
        clone_ids: tuple[str, ...] | None,
        expected_revision: int,
        event_id: str,
        actor: str,
        occurred_at: datetime,
    ) -> SYNProjectSnapshot:
        self._check_revision(snapshot, expected_revision)
        self._require_assembly_colony_step(snapshot)
        latest = _latest_by_clone(snapshot.colonies)
        previously_selected = {
            record.clone_id
            for record in snapshot.prep_records
            if record.selected_for_prep
        }
        positive = sorted(
            (
                record
                for record in latest.values()
                if record.result is SYNColonyPCRResult.POSITIVE
                and record.clone_id not in previously_selected
            ),
            key=lambda record: (
                self._attempt_round(snapshot, record.attempt_id),
                record.clone_no,
            ),
        )
        selected_ids = clone_ids
        if selected_ids is None:
            selected_ids = tuple(record.clone_id for record in positive[:3])
        if not selected_ids:
            raise SYNWorkflowRuleError("至少需要选择 1 个菌落 PCR 阳性克隆")
        if len(selected_ids) != len(set(selected_ids)):
            raise SYNWorkflowRuleError("不能重复选择同一个克隆")
        invalid = tuple(
            clone_id
            for clone_id in selected_ids
            if clone_id not in latest
            or latest[clone_id].result is not SYNColonyPCRResult.POSITIVE
            or clone_id in previously_selected
        )
        if invalid:
            raise SYNWorkflowRuleError("只有菌落 PCR 阳性克隆可以进入小提")

        prep_records = tuple(
            PlasmidPrepRecord(
                clone_id=clone_id,
                selected_for_prep=True,
                status=PlasmidPrepStatus.PENDING,
            )
            for clone_id in selected_ids
        )
        updated = replace(
            snapshot,
            status="plasmid_prep_in_progress",
            prep_records=snapshot.prep_records + prep_records,
        )
        return self._append_event(
            updated,
            event_id=event_id,
            event_type="select_clones_for_prep",
            actor=actor,
            occurred_at=occurred_at,
            from_status=snapshot.status,
            note=f"选择 {len(selected_ids)} 个阳性克隆",
        )

    def record_plasmid_prep(
        self,
        snapshot: SYNProjectSnapshot,
        *,
        clone_id: str,
        status: PlasmidPrepStatus,
        expected_revision: int,
        event_id: str,
        actor: str,
        occurred_at: datetime,
        note: str | None = None,
    ) -> SYNProjectSnapshot:
        self._check_revision(snapshot, expected_revision)
        if snapshot.status != "plasmid_prep_in_progress":
            raise SYNWorkflowRuleError("项目当前不在质粒抽提中")
        if not isinstance(status, PlasmidPrepStatus):
            raise TypeError("status must be a PlasmidPrepStatus")
        latest = _latest_by_clone(snapshot.prep_records)
        if clone_id not in latest or not latest[clone_id].selected_for_prep:
            raise SYNWorkflowRuleError("该克隆未被选中用于质粒抽提")
        record = PlasmidPrepRecord(
            clone_id=clone_id,
            selected_for_prep=True,
            status=status,
            completed_at=(
                occurred_at if status is PlasmidPrepStatus.COMPLETED else None
            ),
            note=note,
        )
        updated = replace(
            snapshot,
            prep_records=snapshot.prep_records + (record,),
        )
        latest_after = _latest_by_clone(updated.prep_records)
        selected = (
            item for item in latest_after.values() if item.selected_for_prep
        )
        all_completed = all(
            item.status is PlasmidPrepStatus.COMPLETED for item in selected
        )
        if all_completed:
            updated = replace(
                updated,
                status="awaiting_sequencing_confirmation",
                syn_assembly_substep=None,
            )
        return self._append_event(
            updated,
            event_id=event_id,
            event_type=(
                "finish_plasmid_prep" if all_completed else "record_plasmid_prep"
            ),
            actor=actor,
            occurred_at=occurred_at,
            from_status=snapshot.status,
            note=status.value,
            related_entity_id=clone_id,
        )

    def preview_additional_screening(
        self,
        snapshot: SYNProjectSnapshot,
        *,
        target_name: str,
        expected_revision: int,
        preview_id: str,
        colony_count: int = 8,
    ) -> AdditionalScreeningPreview:
        self._check_revision(snapshot, expected_revision)
        if snapshot.status != "awaiting_sequencing_confirmation":
            raise SYNWorkflowRuleError("只有待测序确认项目可以预览追加筛选")
        if self.get_syn_sequencing_summary(snapshot).correct_count:
            raise SYNWorkflowRuleError("项目已有正确克隆，不需要追加筛选")
        if colony_count <= 0:
            raise ValueError("colony_count must be positive")
        attempt = self._active_attempt(snapshot)
        max_clone_no = max(
            (
                record.clone_no
                for record in snapshot.colonies
                if record.attempt_id == attempt.attempt_id
            ),
            default=0,
        )
        clone_numbers = tuple(
            range(max_clone_no + 1, max_clone_no + colony_count + 1),
        )
        colonies = self._make_colonies(
            attempt_id=attempt.attempt_id,
            round_no=attempt.syn_assembly_round_no,
            target_name=target_name,
            clone_numbers=clone_numbers,
        )
        return AdditionalScreeningPreview(
            preview_id=preview_id,
            project_id=snapshot.project_id,
            expected_revision=snapshot.revision,
            attempt_id=attempt.attempt_id,
            target_name=target_name,
            clone_numbers=clone_numbers,
            clone_ids=tuple(record.clone_id for record in colonies),
            display_names=tuple(record.display_name for record in colonies),
        )

    def confirm_additional_screening(
        self,
        snapshot: SYNProjectSnapshot,
        preview: AdditionalScreeningPreview,
        *,
        expected_revision: int,
        event_id: str,
        actor: str,
        occurred_at: datetime,
    ) -> SYNProjectSnapshot:
        self._check_revision(snapshot, expected_revision)
        if (
            preview.project_id != snapshot.project_id
            or preview.expected_revision != snapshot.revision
        ):
            raise SYNRevisionConflict("追加筛选预览已过期")
        if snapshot.status != "awaiting_sequencing_confirmation":
            raise SYNWorkflowRuleError("当前状态不能确认追加筛选")
        colonies = self._make_colonies(
            attempt_id=preview.attempt_id,
            round_no=snapshot.syn_assembly_round_no,
            target_name=preview.target_name,
            clone_numbers=preview.clone_numbers,
        )
        if tuple(item.clone_id for item in colonies) != preview.clone_ids:
            raise SYNWorkflowRuleError("追加筛选预览内容不一致")
        updated = replace(
            snapshot,
            status="syn_assembly_in_progress",
            syn_assembly_substep="colony_pcr",
            colonies=snapshot.colonies + colonies,
        )
        return self._append_event(
            updated,
            event_id=event_id,
            event_type="additional_screening",
            actor=actor,
            occurred_at=occurred_at,
            from_status=snapshot.status,
            note=f"追加 {len(colonies)} 个菌落",
            related_entity_id=preview.preview_id,
        )

    def restart_assembly(
        self,
        snapshot: SYNProjectSnapshot,
        *,
        restart_from_substep: str,
        expected_revision: int,
        attempt_id: str,
        event_id: str,
        actor: str,
        occurred_at: datetime,
    ) -> SYNProjectSnapshot:
        self._check_revision(snapshot, expected_revision)
        if snapshot.status != "awaiting_sequencing_confirmation":
            raise SYNWorkflowRuleError("当前状态不能重新合成组装")
        if self.get_syn_sequencing_summary(snapshot).correct_count:
            raise SYNWorkflowRuleError("项目已有正确克隆，不需要重新合成组装")
        if restart_from_substep not in ASSEMBLY_SUBSTEPS:
            raise SYNWorkflowRuleError("重启步骤必须是四个合成组装子步骤之一")
        round_no = snapshot.syn_assembly_round_no + 1
        attempt = SYNAssemblyAttemptRecord(
            attempt_id=attempt_id,
            project_id=snapshot.project_id,
            design_version_id=snapshot.active_design_version_id,
            syn_assembly_round_no=round_no,
            restart_from_substep=restart_from_substep,
            result=SYNAssemblyAttemptResult.PENDING,
            started_at=occurred_at,
        )
        updated = snapshot.append_attempt(attempt)
        updated = replace(
            updated,
            status="syn_assembly_in_progress",
            syn_assembly_substep=restart_from_substep,
        )
        return self._append_event(
            updated,
            event_id=event_id,
            event_type="restart_assembly",
            actor=actor,
            occurred_at=occurred_at,
            from_status=snapshot.status,
            note=f"从 {restart_from_substep} 开始第 {round_no} 轮",
            related_entity_id=attempt_id,
        )

    def confirm_sequencing(
        self,
        snapshot: SYNProjectSnapshot,
        *,
        clone_id: str,
        result: SYNSequencingResult,
        note: str | None,
        expected_revision: int,
        confirmation_id: str,
        event_id: str,
        actor: str,
        occurred_at: datetime,
    ) -> SYNProjectSnapshot:
        self._check_revision(snapshot, expected_revision)
        if snapshot.status != "awaiting_sequencing_confirmation":
            raise SYNWorkflowRuleError("项目当前不在待测序确认状态")
        if not isinstance(result, SYNSequencingResult):
            raise TypeError("result must be a SYNSequencingResult")
        latest_prep = _latest_by_clone(snapshot.prep_records)
        prep = latest_prep.get(clone_id)
        if prep is None or prep.status is not PlasmidPrepStatus.COMPLETED:
            raise SYNWorkflowRuleError("只有完成小提的克隆可以进行测序确认")
        if any(
            item.confirmation_id == confirmation_id
            for item in snapshot.sequencing_confirmations
        ):
            raise SYNWorkflowRuleError(f"确认记录已存在：{confirmation_id}")
        latest_confirmations = _latest_by_clone(snapshot.sequencing_confirmations)
        previous = latest_confirmations.get(clone_id)
        confirmation = SYNSequencingConfirmation(
            confirmation_id=confirmation_id,
            clone_id=clone_id,
            result=result,
            confirmed_at=occurred_at,
            actor=actor,
            note=note,
            supersedes_confirmation_id=(
                previous.confirmation_id if previous is not None else None
            ),
        )
        updated = replace(
            snapshot,
            sequencing_confirmations=(
                snapshot.sequencing_confirmations + (confirmation,)
            ),
        )
        return self._append_event(
            updated,
            event_id=event_id,
            event_type="confirm_sequencing",
            actor=actor,
            occurred_at=occurred_at,
            from_status=snapshot.status,
            note=result.value,
            related_entity_id=confirmation_id,
        )

    def get_syn_sequencing_summary(
        self,
        snapshot: SYNProjectSnapshot,
    ) -> SYNSequencingSummary:
        latest_prep = _latest_by_clone(snapshot.prep_records)
        eligible_clone_ids = tuple(
            clone_id
            for clone_id, prep in latest_prep.items()
            if prep.selected_for_prep
            and prep.status is PlasmidPrepStatus.COMPLETED
        )
        latest_confirmations = _latest_by_clone(snapshot.sequencing_confirmations)
        counts = {
            SYNSequencingResult.CORRECT: 0,
            SYNSequencingResult.INCORRECT: 0,
            SYNSequencingResult.UNCERTAIN: 0,
        }
        usable_clone_ids = []
        unconfirmed_count = 0
        for clone_id in eligible_clone_ids:
            confirmation = latest_confirmations.get(clone_id)
            if confirmation is None:
                unconfirmed_count += 1
                continue
            counts[confirmation.result] += 1
            if confirmation.result is SYNSequencingResult.CORRECT:
                usable_clone_ids.append(clone_id)

        colony_names = {
            clone_id: colony.display_name
            for clone_id, colony in _latest_by_clone(snapshot.colonies).items()
        }
        usable_names = tuple(
            colony_names.get(clone_id, clone_id) for clone_id in usable_clone_ids
        )
        return SYNSequencingSummary(
            correct_count=counts[SYNSequencingResult.CORRECT],
            incorrect_count=counts[SYNSequencingResult.INCORRECT],
            uncertain_count=counts[SYNSequencingResult.UNCERTAIN],
            unconfirmed_count=unconfirmed_count,
            usable_clone_ids=tuple(usable_clone_ids),
            usable_clone_names=usable_names,
        )

    def complete_project(
        self,
        snapshot: SYNProjectSnapshot,
        *,
        user_confirmed: bool,
        expected_revision: int,
        event_id: str,
        actor: str,
        occurred_at: datetime,
    ) -> SYNProjectSnapshot:
        self._check_revision(snapshot, expected_revision)
        if snapshot.status != "awaiting_sequencing_confirmation":
            raise SYNWorkflowRuleError("当前状态不能完成项目")
        if not user_confirmed:
            raise SYNWorkflowRuleError("需要用户显式确认完成项目")
        summary = self.get_syn_sequencing_summary(snapshot)
        if not summary.correct_count:
            raise SYNWorkflowRuleError("至少需要 1 个人工确认正确克隆")
        updated = replace(
            snapshot,
            status="project_completed",
            actual_completed_at=occurred_at,
        )
        return self._append_event(
            updated,
            event_id=event_id,
            event_type="complete_project",
            actor=actor,
            occurred_at=occurred_at,
            from_status=snapshot.status,
            note=f"可用克隆：{', '.join(summary.usable_clone_names)}",
        )

    def _active_attempt(
        self,
        snapshot: SYNProjectSnapshot,
    ) -> SYNAssemblyAttemptRecord:
        attempts = (
            attempt
            for attempt in snapshot.attempts
            if attempt.syn_assembly_round_no == snapshot.syn_assembly_round_no
        )
        attempt = next(attempts, None)
        if attempt is None:
            raise SYNWorkflowRuleError("当前组装轮次没有实验记录")
        return attempt

    def _attempt_round(self, snapshot: SYNProjectSnapshot, attempt_id: str) -> int:
        for attempt in snapshot.attempts:
            if attempt.attempt_id == attempt_id:
                return attempt.syn_assembly_round_no
        raise SYNWorkflowRuleError(f"找不到菌落所属轮次：{attempt_id}")
