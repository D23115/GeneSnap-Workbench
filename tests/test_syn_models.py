import unittest
from dataclasses import FrozenInstanceError, replace
from datetime import datetime, timezone
from decimal import Decimal

from genesnap_workbench.domain.syn import (
    OligoMixItem,
    OligoMixPlan,
    OligoResuspensionItem,
    OligoResuspensionPlan,
    PlasmidPrepRecord,
    PlasmidPrepStatus,
    ResuspensionStatus,
    SYNAssemblyAttemptRecord,
    SYNAssemblyAttemptResult,
    SYNAssemblyOligo,
    SYNAuditEvent,
    SYNColonyPCRRecord,
    SYNColonyPCRResult,
    SYNDesignInput,
    SYNDesignVersion,
    SYNModule,
    SYNModulePlan,
    SYNManualOverrideRecord,
    SYNPlasmidSimulation,
    SYNProjectSnapshot,
    SYNQCRisk,
    SYNRoute,
    SYNSequenceQCResult,
    SYNSequencingConfirmation,
    SYNSequencingResult,
    SYNThermodynamicMetadata,
)


NOW = datetime(2026, 7, 12, 8, 0, tzinfo=timezone.utc)


def make_oligo(
    oligo_id: str,
    *,
    design_version_id: str = "design-v1",
    sequence: str = "A" * 60,
    start: int = 0,
    end: int = 60,
) -> SYNAssemblyOligo:
    return SYNAssemblyOligo(
        oligo_id=oligo_id,
        design_version_id=design_version_id,
        name=oligo_id,
        sequence=sequence,
        strand="forward",
        start=start,
        end=end,
        pool_id="pool-1",
        module_id="module-1",
        overlap_left=None,
        overlap_right=None,
        tm_metadata=SYNThermodynamicMetadata(
            analyzer_name="biopython",
            analyzer_version="1.87",
            tm_celsius=Decimal("60.5"),
        ),
    )


def make_module_plan(
    *,
    design_version_id: str = "design-v1",
    oligo_ids: tuple[str, ...] = ("oligo-1", "oligo-2"),
) -> SYNModulePlan:
    return SYNModulePlan(
        design_version_id=design_version_id,
        route=SYNRoute.SINGLE_POOL,
        modules=(
            SYNModule(
                design_version_id=design_version_id,
                module_id="module-1",
                ordinal=1,
                start=0,
                end=120,
                sequence_checksum="module-checksum",
                oligo_ids=oligo_ids,
                boundary_reason="single_pool",
            ),
        ),
        requires_confirmation=False,
        routing_reason="length_at_or_below_1000_bp",
    )


def make_design_version() -> SYNDesignVersion:
    oligos = (
        make_oligo("oligo-1", start=0, end=60),
        make_oligo("oligo-2", start=60, end=120),
    )
    return SYNDesignVersion(
        design_version_id="design-v1",
        project_id="SYN-001",
        version_no=1,
        created_at=NOW,
        raw_sequence_checksum="raw-checksum",
        normalized_sequence="A" * 120,
        normalized_checksum="normalized-checksum",
        final_sequence="A" * 120,
        final_checksum="final-checksum",
        qc_result=SYNSequenceQCResult(
            design_version_id="design-v1",
            rules_version="syn-qc-v1",
            sequence_checksum="final-checksum",
            sequence_length=120,
            overall_gc_percent=Decimal("0"),
            risks=(
                SYNQCRisk(
                    rule_key="example",
                    severity="info",
                    start=0,
                    end=1,
                    observed_value="none",
                    message="example",
                    requires_confirmation=False,
                ),
            ),
            blocked_reasons=(),
            confirmable_warnings=(),
        ),
        module_plan=make_module_plan(),
        oligos=oligos,
        plasmid_simulation=SYNPlasmidSimulation(
            design_version_id="design-v1",
            vector_record_id="vector-1",
            vector_checksum="vector-checksum",
            protocol_version_id="puc57-ecorv-v1",
            linearization_sites=("EcoRV",),
            site_retention_rule="not_required",
            homology_arms=("AAAA", "TTTT"),
            junctions=("left", "right"),
            expected_plasmid_sequence="C" * 200,
            expected_plasmid_checksum="plasmid-checksum",
        ),
        rule_versions=("syn-qc-v1", "syn-oligo-v1"),
        manual_overrides=(),
    )


class SYNDesignModelTests(unittest.TestCase):
    def test_design_input_rejects_blank_required_values(self):
        valid = {
            "project_id": "SYN-001",
            "target_name": "example-syn",
            "raw_sequence": "ACGT",
            "input_format": "plain",
            "vector_protocol_version_id": "puc57-ecorv-v1",
        }

        for field_name in (
            "project_id",
            "target_name",
            "raw_sequence",
            "vector_protocol_version_id",
        ):
            with self.subTest(field_name=field_name):
                values = valid | {field_name: "   "}
                with self.assertRaisesRegex(ValueError, field_name):
                    SYNDesignInput(**values)

    def test_design_input_allows_optional_gene_and_species(self):
        design_input = SYNDesignInput(
            project_id="SYN-001",
            target_name="example-syn",
            raw_sequence="ACGT",
            input_format="plain",
            vector_protocol_version_id="puc57-ecorv-v1",
        )

        self.assertIsNone(design_input.gene_symbol)
        self.assertIsNone(design_input.species)

    def test_assembly_oligo_rejects_sequence_over_65_nt(self):
        with self.assertRaisesRegex(ValueError, "65"):
            make_oligo("too-long", sequence="A" * 66, end=66)

    def test_module_rejects_invalid_half_open_interval(self):
        with self.assertRaisesRegex(ValueError, "interval"):
            SYNModule(
                design_version_id="design-v1",
                module_id="module-1",
                ordinal=1,
                start=20,
                end=20,
                sequence_checksum="checksum",
                oligo_ids=("oligo-1",),
                boundary_reason="example",
            )

    def test_module_plan_rejects_non_contiguous_modules(self):
        module_1 = SYNModule(
            design_version_id="design-v1",
            module_id="module-1",
            ordinal=1,
            start=0,
            end=100,
            sequence_checksum="checksum-1",
            oligo_ids=("oligo-1",),
            boundary_reason="example",
        )
        module_2 = SYNModule(
            design_version_id="design-v1",
            module_id="module-2",
            ordinal=2,
            start=101,
            end=200,
            sequence_checksum="checksum-2",
            oligo_ids=("oligo-2",),
            boundary_reason="example",
        )

        with self.assertRaisesRegex(ValueError, "contiguous"):
            SYNModulePlan(
                design_version_id="design-v1",
                route=SYNRoute.MODULAR,
                modules=(module_1, module_2),
                requires_confirmation=False,
                routing_reason="over_1200_bp",
            )

    def test_design_version_rejects_nested_data_from_another_version(self):
        design = make_design_version()
        mismatched_qc = SYNSequenceQCResult(
            design_version_id="design-v2",
            rules_version="syn-qc-v1",
            sequence_checksum="final-checksum",
            sequence_length=120,
            overall_gc_percent=Decimal("0"),
            risks=(),
            blocked_reasons=(),
            confirmable_warnings=(),
        )

        with self.assertRaisesRegex(ValueError, "design_version_id"):
            SYNDesignVersion(
                design_version_id=design.design_version_id,
                project_id=design.project_id,
                version_no=design.version_no,
                created_at=design.created_at,
                raw_sequence_checksum=design.raw_sequence_checksum,
                normalized_sequence=design.normalized_sequence,
                normalized_checksum=design.normalized_checksum,
                final_sequence=design.final_sequence,
                final_checksum=design.final_checksum,
                qc_result=mismatched_qc,
                module_plan=design.module_plan,
                oligos=design.oligos,
                plasmid_simulation=design.plasmid_simulation,
                rule_versions=design.rule_versions,
                manual_overrides=design.manual_overrides,
            )

    def test_design_version_rejects_module_reference_to_missing_oligo(self):
        design = make_design_version()

        with self.assertRaisesRegex(ValueError, "missing oligo"):
            SYNDesignVersion(
                design_version_id=design.design_version_id,
                project_id=design.project_id,
                version_no=design.version_no,
                created_at=design.created_at,
                raw_sequence_checksum=design.raw_sequence_checksum,
                normalized_sequence=design.normalized_sequence,
                normalized_checksum=design.normalized_checksum,
                final_sequence=design.final_sequence,
                final_checksum=design.final_checksum,
                qc_result=design.qc_result,
                module_plan=make_module_plan(oligo_ids=("missing",)),
                oligos=design.oligos,
                plasmid_simulation=design.plasmid_simulation,
                rule_versions=design.rule_versions,
                manual_overrides=design.manual_overrides,
            )

    def test_design_version_rejects_qc_for_another_sequence_checksum(self):
        design = make_design_version()

        with self.assertRaisesRegex(ValueError, "final_checksum"):
            replace(
                design,
                qc_result=replace(
                    design.qc_result,
                    sequence_checksum="another-checksum",
                ),
            )

    def test_design_version_rejects_module_plan_that_does_not_cover_sequence(self):
        design = make_design_version()
        short_module = replace(design.module_plan.modules[0], end=119)

        with self.assertRaisesRegex(ValueError, "complete final sequence"):
            replace(
                design,
                module_plan=replace(design.module_plan, modules=(short_module,)),
            )

    def test_manual_override_keeps_old_new_values_and_reason_separate(self):
        override = SYNManualOverrideRecord(
            override_id="override-1",
            field_path="module_plan.modules[0].end",
            old_value="120",
            new_value="118",
            reason="人工调整到低风险边界",
            occurred_at=NOW,
            actor="tester",
        )

        self.assertEqual(override.old_value, "120")
        self.assertEqual(override.new_value, "118")
        self.assertEqual(override.reason, "人工调整到低风险边界")


class SYNMaterialModelTests(unittest.TestCase):
    def test_mix_plan_without_standard_volume_is_not_export_ready(self):
        plan = OligoMixPlan(
            design_version_id="design-v1",
            standard_volume_per_oligo_ul=None,
            items=(
                OligoMixItem(
                    pool_id="pool-1",
                    oligo_id="oligo-1",
                    reference_concentration_uM=Decimal("100"),
                    actual_concentration_uM=Decimal("100"),
                    sample_volume_ul=None,
                ),
            ),
            generated_at=NOW,
        )

        self.assertFalse(plan.is_formal_export_ready)

    def test_resuspension_plan_keeps_planned_and_actual_values_separate(self):
        plan = OligoResuspensionPlan(
            design_version_id="design-v1",
            items=(
                OligoResuspensionItem(
                    oligo_id="oligo-1",
                    planned_amount_nmol=Decimal("25"),
                    actual_amount_nmol=Decimal("27.4"),
                    target_stock_concentration_uM=Decimal("100"),
                    water_volume_ul=Decimal("274"),
                    actual_stock_concentration_uM=Decimal("100"),
                ),
            ),
            status=ResuspensionStatus.COMPLETE,
            generated_at=NOW,
        )

        item = plan.items[0]
        self.assertEqual(item.planned_amount_nmol, Decimal("25"))
        self.assertEqual(item.actual_amount_nmol, Decimal("27.4"))


class SYNExperimentRecordTests(unittest.TestCase):
    def make_attempt(self, attempt_id: str, round_no: int) -> SYNAssemblyAttemptRecord:
        return SYNAssemblyAttemptRecord(
            attempt_id=attempt_id,
            project_id="SYN-001",
            design_version_id="design-v1",
            syn_assembly_round_no=round_no,
            restart_from_substep="assembly_pcr",
            result=SYNAssemblyAttemptResult.PENDING,
            started_at=NOW,
        )

    def make_snapshot(self) -> SYNProjectSnapshot:
        return SYNProjectSnapshot(
            project_id="SYN-001",
            revision=1,
            status="syn_assembly_in_progress",
            resuspension_data_status="complete",
            syn_assembly_round_no=1,
            syn_assembly_substep="assembly_pcr",
            active_design_version_id="design-v1",
            attempts=(self.make_attempt("attempt-1", 1),),
            colonies=(
                SYNColonyPCRRecord(
                    clone_id="clone-1",
                    attempt_id="attempt-1",
                    clone_no=1,
                    display_name="example-C01",
                    result=SYNColonyPCRResult.PENDING,
                ),
            ),
            prep_records=(
                PlasmidPrepRecord(
                    clone_id="clone-1",
                    selected_for_prep=False,
                    status=PlasmidPrepStatus.PENDING,
                ),
            ),
            sequencing_confirmations=(
                SYNSequencingConfirmation(
                    confirmation_id="confirmation-1",
                    clone_id="clone-1",
                    result=SYNSequencingResult.UNCERTAIN,
                    confirmed_at=NOW,
                    actor="tester",
                ),
            ),
            status_history=(
                SYNAuditEvent(
                    event_id="event-1",
                    event_type="start_assembly",
                    occurred_at=NOW,
                    actor="tester",
                ),
            ),
            manual_override_history=(),
        )

    def test_snapshot_appends_attempt_without_overwriting_history(self):
        snapshot = self.make_snapshot()

        updated = snapshot.append_attempt(self.make_attempt("attempt-2", 2))

        self.assertEqual(snapshot.revision, 1)
        self.assertEqual(tuple(item.attempt_id for item in snapshot.attempts), ("attempt-1",))
        self.assertEqual(updated.revision, 2)
        self.assertEqual(
            tuple(item.attempt_id for item in updated.attempts),
            ("attempt-1", "attempt-2"),
        )

    def test_snapshot_rejects_duplicate_attempt_id(self):
        snapshot = self.make_snapshot()

        with self.assertRaisesRegex(ValueError, "attempt-1"):
            snapshot.append_attempt(self.make_attempt("attempt-1", 2))

    def test_snapshot_rejects_attempt_from_another_project_or_design(self):
        snapshot = self.make_snapshot()
        attempt = self.make_attempt("attempt-2", 2)

        with self.assertRaisesRegex(ValueError, "project_id"):
            snapshot.append_attempt(replace(attempt, project_id="SYN-OTHER"))
        with self.assertRaisesRegex(ValueError, "design_version_id"):
            snapshot.append_attempt(
                replace(attempt, design_version_id="design-v2"),
            )

    def test_snapshot_rejects_attempt_that_does_not_advance_round(self):
        snapshot = self.make_snapshot()

        with self.assertRaisesRegex(ValueError, "advance"):
            snapshot.append_attempt(self.make_attempt("attempt-2", 1))

    def test_snapshot_rejects_loaded_attempt_from_another_project(self):
        snapshot = self.make_snapshot()
        invalid_attempt = replace(snapshot.attempts[0], project_id="SYN-OTHER")

        with self.assertRaisesRegex(ValueError, "project_id"):
            replace(snapshot, attempts=(invalid_attempt,))

    def test_snapshot_appends_status_event_without_overwriting_history(self):
        snapshot = self.make_snapshot()
        event = SYNAuditEvent(
            event_id="event-2",
            event_type="advance_substep",
            occurred_at=NOW,
            actor="tester",
        )

        updated = snapshot.append_status_event(event)

        self.assertEqual(
            tuple(item.event_id for item in snapshot.status_history),
            ("event-1",),
        )
        self.assertEqual(
            tuple(item.event_id for item in updated.status_history),
            ("event-1", "event-2"),
        )
        self.assertEqual(updated.revision, 2)

    def test_snapshot_rejects_duplicate_status_event_id(self):
        snapshot = self.make_snapshot()

        with self.assertRaisesRegex(ValueError, "event-1"):
            snapshot.append_status_event(snapshot.status_history[0])

    def test_snapshot_is_immutable(self):
        snapshot = self.make_snapshot()

        with self.assertRaises(FrozenInstanceError):
            snapshot.status = "project_completed"


if __name__ == "__main__":
    unittest.main()
