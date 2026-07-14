"""Shared interruption and effective-due-date audit records."""

from __future__ import annotations

from dataclasses import dataclass
from datetime import datetime


@dataclass(frozen=True, slots=True)
class ProjectInterruptionRecord:
    interruption_id: str
    interruption_type: str
    started_at: datetime
    resumed_at: datetime
    previous_status: str
    resume_status: str
    frozen_remaining_workdays: int
    paused_workdays: int
    start_note: str
    resume_note: str

    def __post_init__(self) -> None:
        if self.interruption_type not in {"pause", "abnormal"}:
            raise ValueError("中断类型必须是 pause 或 abnormal")
        if self.resumed_at < self.started_at:
            raise ValueError("恢复时间不能早于中断时间")
        if self.paused_workdays < 0:
            raise ValueError("暂停工作日不能为负数")
