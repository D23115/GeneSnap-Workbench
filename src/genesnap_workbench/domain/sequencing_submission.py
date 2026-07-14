"""Append-only sequencing submission records shared by molecular workflows."""

from __future__ import annotations

from dataclasses import dataclass
from datetime import datetime


@dataclass(frozen=True, slots=True)
class SequencingSubmissionRecord:
    submission_id: str
    round_no: int
    submission_kind: str
    created_at: datetime
    sent_at: datetime
    sample_names: tuple[str, ...]
    internal_submission_no: str = ""
    vendor_order_no: str = ""
    template_id: str | None = None
    form_path: str | None = None
    experiment_attempt_no: int = 1
    status: str = "sent"
    note: str | None = None

    def __post_init__(self) -> None:
        if not self.submission_id.strip():
            raise ValueError("submission_id 不能为空")
        if self.round_no <= 0:
            raise ValueError("round_no 必须大于 0")
        if self.submission_kind not in {"initial", "add_on", "post_rework", "retest"}:
            raise ValueError("不支持的送测记录类型")
        if not self.sample_names or len(self.sample_names) != len(set(self.sample_names)):
            raise ValueError("送测记录必须包含不重复的样本名")
        if self.experiment_attempt_no <= 0:
            raise ValueError("experiment_attempt_no 必须大于 0")
        if self.status not in {"planned", "sent", "analyzed"}:
            raise ValueError("不支持的送测记录状态")
