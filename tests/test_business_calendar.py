import unittest
from datetime import date

from genesnap_workbench.project_workflow.business_calendar import (
    ChinaBusinessCalendar,
)


class ChinaBusinessCalendarTests(unittest.TestCase):
    def setUp(self):
        self.calendar = ChinaBusinessCalendar.for_2026()

    def test_adjusted_sunday_is_workday_and_holiday_weekday_is_not(self):
        self.assertTrue(self.calendar.is_workday(date(2026, 1, 4)))
        self.assertFalse(self.calendar.is_workday(date(2026, 2, 16)))

    def test_add_workdays_skips_holidays_and_uses_adjusted_workday(self):
        self.assertEqual(
            self.calendar.add_workdays(date(2025, 12, 31), 2),
            date(2026, 1, 5),
        )

    def test_remaining_workdays_is_zero_when_due_today(self):
        self.assertEqual(
            self.calendar.remaining_workdays(
                date(2026, 7, 13),
                date(2026, 7, 13),
            ),
            0,
        )


if __name__ == "__main__":
    unittest.main()
