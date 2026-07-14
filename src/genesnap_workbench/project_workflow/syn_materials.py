"""Oligo resuspension, equimolar mixing, and readiness checks for SYN."""

from __future__ import annotations

from dataclasses import dataclass
from datetime import datetime, timezone
from decimal import Decimal

from genesnap_workbench.domain.syn import (
    OligoMixItem,
    OligoMixPlan,
    OligoResuspensionItem,
    OligoResuspensionPlan,
    ResuspensionStatus,
    SYNAssemblyOligo,
)


def _now() -> datetime:
    return datetime.now(timezone.utc)


def _require_positive(field_name: str, value: Decimal) -> None:
    if value <= 0:
        raise ValueError(f"{field_name} must be positive")


def calculate_resuspension_water_volume(
    *,
    amount_nmol: Decimal,
    target_stock_concentration_uM: Decimal,
) -> Decimal:
    """Calculate uL water using nmol * 1000 / target uM."""
    _require_positive("amount_nmol", amount_nmol)
    _require_positive(
        "target_stock_concentration_uM",
        target_stock_concentration_uM,
    )
    return amount_nmol * Decimal("1000") / target_stock_concentration_uM


def create_resuspension_plan(
    oligos: tuple[SYNAssemblyOligo, ...],
    actual_amounts_nmol: dict[str, Decimal],
    target_stock_concentration_uM: Decimal = Decimal("100"),
    *,
    planned_amounts_nmol: dict[str, Decimal] | None = None,
    actual_stock_concentrations_uM: dict[str, Decimal] | None = None,
    generated_at: datetime | None = None,
) -> OligoResuspensionPlan:
    """Build a plan without inventing values for oligos missing actual nmol."""
    if not oligos:
        raise ValueError("At least one oligo is required")
    design_version_ids = {oligo.design_version_id for oligo in oligos}
    if len(design_version_ids) != 1:
        raise ValueError("All oligos must belong to one design version")
    known_ids = {oligo.oligo_id for oligo in oligos}
    unknown_ids = set(actual_amounts_nmol) - known_ids
    if unknown_ids:
        raise ValueError(f"Unknown oligo IDs: {', '.join(sorted(unknown_ids))}")

    planned_amounts_nmol = planned_amounts_nmol or {}
    actual_stock_concentrations_uM = actual_stock_concentrations_uM or {}
    items: list[OligoResuspensionItem] = []
    completed_count = 0
    for oligo in oligos:
        actual_amount = actual_amounts_nmol.get(oligo.oligo_id)
        water_volume = None
        actual_concentration = None
        if actual_amount is not None:
            water_volume = calculate_resuspension_water_volume(
                amount_nmol=actual_amount,
                target_stock_concentration_uM=target_stock_concentration_uM,
            )
            actual_concentration = actual_stock_concentrations_uM.get(
                oligo.oligo_id,
                target_stock_concentration_uM,
            )
            _require_positive(
                "actual_stock_concentration_uM",
                actual_concentration,
            )
            completed_count += 1
        items.append(
            OligoResuspensionItem(
                oligo_id=oligo.oligo_id,
                planned_amount_nmol=planned_amounts_nmol.get(oligo.oligo_id),
                actual_amount_nmol=actual_amount,
                target_stock_concentration_uM=target_stock_concentration_uM,
                water_volume_ul=water_volume,
                actual_stock_concentration_uM=actual_concentration,
            ),
        )

    if completed_count == 0:
        status = ResuspensionStatus.MISSING
    elif completed_count == len(oligos):
        status = ResuspensionStatus.COMPLETE
    else:
        status = ResuspensionStatus.PARTIAL
    return OligoResuspensionPlan(
        design_version_id=next(iter(design_version_ids)),
        items=tuple(items),
        status=status,
        generated_at=generated_at or _now(),
    )


def create_mix_plan(
    resuspension_plan: OligoResuspensionPlan,
    *,
    oligo_pool_ids: dict[str, str],
    standard_volume_per_oligo_ul: Decimal | None,
    generated_at: datetime | None = None,
) -> OligoMixPlan:
    """Create one single-pool equimolar mix plan."""
    plan_oligo_ids = {item.oligo_id for item in resuspension_plan.items}
    missing_pool_ids = plan_oligo_ids - set(oligo_pool_ids)
    if missing_pool_ids:
        raise ValueError(
            f"Missing pool IDs: {', '.join(sorted(missing_pool_ids))}",
        )
    pool_ids = {oligo_pool_ids[oligo_id] for oligo_id in plan_oligo_ids}
    if len(pool_ids) != 1:
        raise ValueError("One OligoMixPlan must contain exactly one pool")
    if standard_volume_per_oligo_ul is not None:
        _require_positive(
            "standard_volume_per_oligo_ul",
            standard_volume_per_oligo_ul,
        )

    items: list[OligoMixItem] = []
    for resuspension_item in resuspension_plan.items:
        reference_concentration = resuspension_item.target_stock_concentration_uM
        actual_concentration = resuspension_item.actual_stock_concentration_uM
        sample_volume = None
        if (
            standard_volume_per_oligo_ul is not None
            and actual_concentration is not None
        ):
            sample_volume = (
                reference_concentration
                * standard_volume_per_oligo_ul
                / actual_concentration
            )
        items.append(
            OligoMixItem(
                pool_id=oligo_pool_ids[resuspension_item.oligo_id],
                oligo_id=resuspension_item.oligo_id,
                reference_concentration_uM=reference_concentration,
                actual_concentration_uM=actual_concentration,
                sample_volume_ul=sample_volume,
            ),
        )
    return OligoMixPlan(
        design_version_id=resuspension_plan.design_version_id,
        standard_volume_per_oligo_ul=standard_volume_per_oligo_ul,
        items=tuple(items),
        generated_at=generated_at or _now(),
    )


@dataclass(frozen=True, slots=True)
class MaterialReadiness:
    is_ready: bool
    can_start_with_override: bool
    missing_oligo_ids: tuple[str, ...]
    missing_fields: tuple[str, ...]
    errors: tuple[str, ...]


def validate_material_readiness(
    resuspension_plan: OligoResuspensionPlan,
    mix_plan: OligoMixPlan,
) -> MaterialReadiness:
    """Return missing-data warnings separately from blocking errors."""
    errors: list[str] = []
    if resuspension_plan.design_version_id != mix_plan.design_version_id:
        errors.append("复溶计划与混池计划的设计版本不一致")
    missing_oligo_ids = tuple(
        item.oligo_id
        for item in resuspension_plan.items
        if item.actual_amount_nmol is None or item.water_volume_ul is None
    )
    missing_fields = (
        ("standard_volume_per_oligo_ul",)
        if mix_plan.standard_volume_per_oligo_ul is None
        else ()
    )
    is_ready = (
        not errors
        and not missing_oligo_ids
        and not missing_fields
        and mix_plan.is_formal_export_ready
    )
    return MaterialReadiness(
        is_ready=is_ready,
        can_start_with_override=not errors,
        missing_oligo_ids=missing_oligo_ids,
        missing_fields=missing_fields,
        errors=tuple(errors),
    )
