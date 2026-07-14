import unittest
from datetime import datetime, timezone
from decimal import Decimal

from genesnap_workbench.domain.syn import (
    SYNAssemblyOligo,
    SYNThermodynamicMetadata,
)
from genesnap_workbench.project_workflow.syn_materials import (
    calculate_resuspension_water_volume,
    create_mix_plan,
    create_resuspension_plan,
    validate_material_readiness,
)


NOW = datetime(2026, 7, 12, 8, 0, tzinfo=timezone.utc)


def make_oligo(oligo_id: str, pool_id: str = "pool-01") -> SYNAssemblyOligo:
    return SYNAssemblyOligo(
        oligo_id=oligo_id,
        design_version_id="design-v1",
        name=oligo_id,
        sequence="ACGT" * 15,
        strand="forward",
        start=0,
        end=60,
        pool_id=pool_id,
        module_id="module-01",
        overlap_left=None,
        overlap_right=None,
        tm_metadata=SYNThermodynamicMetadata(
            analyzer_name="test",
            analyzer_version="1",
            tm_celsius=Decimal("60"),
        ),
    )


class SYNResuspensionTests(unittest.TestCase):
    def test_25_nmol_to_100_uM_requires_250_uL_water(self):
        self.assertEqual(
            calculate_resuspension_water_volume(
                amount_nmol=Decimal("25"),
                target_stock_concentration_uM=Decimal("100"),
            ),
            Decimal("250"),
        )

    def test_plan_keeps_missing_amount_without_inventing_volume(self):
        oligos = (make_oligo("oligo-1"), make_oligo("oligo-2"))

        plan = create_resuspension_plan(
            oligos,
            actual_amounts_nmol={"oligo-1": Decimal("25")},
            target_stock_concentration_uM=Decimal("100"),
            generated_at=NOW,
        )

        self.assertEqual(plan.status.value, "partial")
        items = {item.oligo_id: item for item in plan.items}
        self.assertEqual(items["oligo-1"].water_volume_ul, Decimal("250"))
        self.assertIsNone(items["oligo-2"].actual_amount_nmol)
        self.assertIsNone(items["oligo-2"].water_volume_ul)

    def test_zero_or_negative_amount_is_rejected(self):
        for amount in (Decimal("0"), Decimal("-1")):
            with self.subTest(amount=amount):
                with self.assertRaisesRegex(ValueError, "positive"):
                    calculate_resuspension_water_volume(
                        amount_nmol=amount,
                        target_stock_concentration_uM=Decimal("100"),
                    )


class SYNOligoMixTests(unittest.TestCase):
    def make_resuspension_plan(self):
        return create_resuspension_plan(
            (make_oligo("oligo-1"), make_oligo("oligo-2")),
            actual_amounts_nmol={
                "oligo-1": Decimal("25"),
                "oligo-2": Decimal("25"),
            },
            target_stock_concentration_uM=Decimal("100"),
            actual_stock_concentrations_uM={
                "oligo-1": Decimal("100"),
                "oligo-2": Decimal("50"),
            },
            generated_at=NOW,
        )

    def test_equal_molar_mix_adjusts_only_concentration_exception(self):
        plan = create_mix_plan(
            self.make_resuspension_plan(),
            oligo_pool_ids={"oligo-1": "pool-01", "oligo-2": "pool-01"},
            standard_volume_per_oligo_ul=Decimal("2"),
            generated_at=NOW,
        )

        items = {item.oligo_id: item for item in plan.items}
        self.assertEqual(items["oligo-1"].sample_volume_ul, Decimal("2"))
        self.assertEqual(items["oligo-2"].sample_volume_ul, Decimal("4"))
        self.assertTrue(plan.is_formal_export_ready)

    def test_missing_standard_volume_keeps_formal_mix_not_ready(self):
        plan = create_mix_plan(
            self.make_resuspension_plan(),
            oligo_pool_ids={"oligo-1": "pool-01", "oligo-2": "pool-01"},
            standard_volume_per_oligo_ul=None,
            generated_at=NOW,
        )

        self.assertFalse(plan.is_formal_export_ready)
        self.assertTrue(all(item.sample_volume_ul is None for item in plan.items))

    def test_material_readiness_allows_audited_override_for_missing_data(self):
        resuspension = create_resuspension_plan(
            (make_oligo("oligo-1"), make_oligo("oligo-2")),
            actual_amounts_nmol={"oligo-1": Decimal("25")},
            generated_at=NOW,
        )
        mix = create_mix_plan(
            resuspension,
            oligo_pool_ids={"oligo-1": "pool-01", "oligo-2": "pool-01"},
            standard_volume_per_oligo_ul=None,
            generated_at=NOW,
        )

        readiness = validate_material_readiness(resuspension, mix)

        missing_item = next(item for item in mix.items if item.oligo_id == "oligo-2")
        self.assertIsNone(missing_item.actual_concentration_uM)
        self.assertIsNone(missing_item.sample_volume_ul)
        self.assertFalse(readiness.is_ready)
        self.assertTrue(readiness.can_start_with_override)
        self.assertIn("oligo-2", readiness.missing_oligo_ids)
        self.assertIn("standard_volume_per_oligo_ul", readiness.missing_fields)


if __name__ == "__main__":
    unittest.main()
