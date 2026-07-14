"""Editable business-day rules seeded from the official 2026 China calendar."""

from __future__ import annotations

from dataclasses import dataclass
from datetime import date, timedelta


OFFICIAL_2026_SOURCE = (
    "https://www.gov.cn/zhengce/zhengceku/202511/content_7047091.htm"
)


def _date_range(start: date, end: date) -> set[date]:
    values = set()
    current = start
    while current <= end:
        values.add(current)
        current += timedelta(days=1)
    return values


@dataclass(frozen=True, slots=True)
class ChinaBusinessCalendar:
    holidays: frozenset[date]
    adjusted_workdays: frozenset[date]
    source: str

    @classmethod
    def for_2026(cls) -> "ChinaBusinessCalendar":
        holidays = set()
        for start, end in (
            (date(2026, 1, 1), date(2026, 1, 3)),
            (date(2026, 2, 15), date(2026, 2, 23)),
            (date(2026, 4, 4), date(2026, 4, 6)),
            (date(2026, 5, 1), date(2026, 5, 5)),
            (date(2026, 6, 19), date(2026, 6, 21)),
            (date(2026, 9, 25), date(2026, 9, 27)),
            (date(2026, 10, 1), date(2026, 10, 7)),
        ):
            holidays.update(_date_range(start, end))
        return cls(
            holidays=frozenset(holidays),
            adjusted_workdays=frozenset(
                {
                    date(2026, 1, 4),
                    date(2026, 2, 14),
                    date(2026, 2, 28),
                    date(2026, 5, 9),
                    date(2026, 9, 20),
                    date(2026, 10, 10),
                },
            ),
            source=OFFICIAL_2026_SOURCE,
        )

    def is_workday(self, value: date) -> bool:
        if value in self.adjusted_workdays:
            return True
        if value in self.holidays:
            return False
        return value.weekday() < 5

    def add_workdays(self, start: date, count: int) -> date:
        if count < 0:
            raise ValueError("count must not be negative")
        current = start
        remaining = count
        while remaining:
            current += timedelta(days=1)
            if self.is_workday(current):
                remaining -= 1
        return current

    def remaining_workdays(self, today: date, due_date: date) -> int:
        if due_date == today:
            return 0
        direction = 1 if due_date > today else -1
        current = today
        count = 0
        while current != due_date:
            current += timedelta(days=direction)
            if self.is_workday(current):
                count += direction
        return count

    def workdays_in_half_open_interval(self, start: date, end: date) -> int:
        if end < start:
            raise ValueError("end must not be earlier than start")
        current = start
        count = 0
        while current < end:
            if self.is_workday(current):
                count += 1
            current += timedelta(days=1)
        return count
