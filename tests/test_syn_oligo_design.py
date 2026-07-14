import random
import unittest
from decimal import Decimal

from genesnap_workbench.domain.syn import (
    SYNQCRisk,
    SYNRoute,
    SYNSequenceQCResult,
)
from genesnap_workbench.sequence_core.dna import sha256_sequence
from genesnap_workbench.sequence_core.syn_modules import (
    SYNModuleRules,
    plan_syn_modules,
)
from genesnap_workbench.sequence_core.syn_oligos import (
    BiopythonThermodynamicAnalyzer,
    SYNOligoDesignFailure,
    SYNOligoRules,
    design_assembly_oligos,
    reconstruct_from_assembly_oligos,
)


def deterministic_dna(length: int, seed: int) -> str:
    randomizer = random.Random(seed)
    return "".join(randomizer.choice("ACGT") for _ in range(length))


class SYNModulePlannerTests(unittest.TestCase):
    def test_length_boundaries_choose_expected_route(self):
        cases = (
            (999, SYNRoute.SINGLE_POOL, False),
            (1000, SYNRoute.SINGLE_POOL, False),
            (1001, SYNRoute.SINGLE_POOL, True),
            (1200, SYNRoute.SINGLE_POOL, True),
            (1201, SYNRoute.MODULAR, False),
        )

        for length, route, requires_confirmation in cases:
            with self.subTest(length=length):
                plan = plan_syn_modules(
                    deterministic_dna(length, length),
                    SYNModuleRules(),
                    design_version_id=f"design-{length}",
                )
                self.assertEqual(plan.route, route)
                self.assertEqual(plan.requires_confirmation, requires_confirmation)

    def test_modular_plan_is_contiguous_and_respects_module_size(self):
        sequence = deterministic_dna(2600, 2600)

        plan = plan_syn_modules(
            sequence,
            SYNModuleRules(),
            design_version_id="design-v1",
        )

        self.assertEqual(plan.route, SYNRoute.MODULAR)
        self.assertEqual(plan.modules[0].start, 0)
        self.assertEqual(plan.modules[-1].end, len(sequence))
        self.assertTrue(
            all(500 <= module.end - module.start <= 900 for module in plan.modules),
        )
        self.assertTrue(
            all(
                left.end == right.start
                for left, right in zip(plan.modules, plan.modules[1:])
            ),
        )
        for left, right in zip(plan.modules, plan.modules[1:]):
            self.assertEqual(left.right_overlap, right.left_overlap)
            overlap_start, overlap_end = left.right_overlap
            overlap_sequence = sequence[overlap_start:overlap_end]
            self.assertTrue(20 <= len(overlap_sequence) <= 30)
            self.assertEqual(sequence.count(overlap_sequence), 1)

    def test_modular_plan_rejects_non_unique_module_homology(self):
        with self.assertRaisesRegex(ValueError, "唯一模块同源区"):
            plan_syn_modules(
                "A" * 1300,
                SYNModuleRules(),
                design_version_id="design-v1",
            )

    def test_modular_boundary_moves_out_of_high_risk_region(self):
        sequence = deterministic_dna(1300, 1300)
        qc_result = SYNSequenceQCResult(
            design_version_id="design-v1",
            rules_version="syn-qc-v1",
            sequence_checksum=sha256_sequence(sequence),
            sequence_length=len(sequence),
            overall_gc_percent=Decimal("50"),
            risks=(
                SYNQCRisk(
                    rule_key="local_gc",
                    severity="high_risk",
                    start=640,
                    end=661,
                    observed_value="90%",
                    message="高风险边界测试区",
                    requires_confirmation=True,
                ),
            ),
            blocked_reasons=(),
            confirmable_warnings=("高风险边界测试区",),
        )

        plan = plan_syn_modules(
            sequence,
            SYNModuleRules(),
            design_version_id="design-v1",
            qc_result=qc_result,
        )

        boundary = plan.modules[0].end
        self.assertFalse(640 <= boundary < 661)
        self.assertEqual(plan.modules[0].boundary_reason, "low_risk_adjusted_boundary")


class SYNOligoDesignerTests(unittest.TestCase):
    def setUp(self):
        self.thermodynamics = BiopythonThermodynamicAnalyzer()

    def design(
        self,
        sequence: str,
        rules: SYNOligoRules,
    ):
        module_plan = plan_syn_modules(
            sequence,
            SYNModuleRules(),
            design_version_id="design-v1",
        )
        return design_assembly_oligos(
            sequence,
            module_plan,
            rules,
            self.thermodynamics,
            design_version_id="design-v1",
            project_id="SYN-001",
            target_name="example",
        )

    def test_rules_reject_overlap_26_and_oligo_limit_66(self):
        with self.assertRaisesRegex(ValueError, "25"):
            SYNOligoRules(max_overlap_length=26)
        with self.assertRaisesRegex(ValueError, "65"):
            SYNOligoRules(max_oligo_length=66)

    def test_supports_overlap_lengths_18_20_and_25(self):
        for overlap_length in (18, 20, 25):
            with self.subTest(overlap_length=overlap_length):
                sequence_length = 4 * 59 - 3 * overlap_length
                result = self.design(
                    deterministic_dna(sequence_length, overlap_length),
                    SYNOligoRules(
                        target_overlap_length=overlap_length,
                        min_overlap_length=overlap_length,
                        max_overlap_length=overlap_length,
                    ),
                )
                self.assertEqual(
                    {overlap.end - overlap.start for overlap in result.overlaps},
                    {overlap_length},
                )

    def test_prefers_58_to_60_nt_and_reconstructs_original_sequence(self):
        sequence = deterministic_dna(176, 176)

        result = self.design(sequence, SYNOligoRules())

        self.assertTrue(all(58 <= len(oligo.sequence) <= 60 for oligo in result.oligos))
        self.assertTrue(all(len(oligo.sequence) <= 65 for oligo in result.oligos))
        self.assertEqual(
            tuple(oligo.strand for oligo in result.oligos),
            ("forward", "reverse", "forward", "reverse"),
        )
        self.assertEqual(
            reconstruct_from_assembly_oligos(result.oligos),
            sequence,
        )
        self.assertTrue(all(overlap.occurrence_count == 1 for overlap in result.overlaps))
        self.assertTrue(all(module.oligo_ids for module in result.module_plan.modules))

    def test_uses_61_to_65_nt_only_when_preferred_range_has_no_solution(self):
        sequence = deterministic_dna(185, 185)

        result = self.design(
            sequence,
            SYNOligoRules(
                target_overlap_length=20,
                min_overlap_length=20,
                max_overlap_length=20,
            ),
        )

        self.assertTrue(any(len(oligo.sequence) > 60 for oligo in result.oligos))
        self.assertTrue(all(len(oligo.sequence) <= 65 for oligo in result.oligos))

    def test_modular_oligos_reconstruct_complete_target(self):
        sequence = deterministic_dna(1301, 1301)

        result = self.design(sequence, SYNOligoRules())

        self.assertEqual(result.module_plan.route, SYNRoute.MODULAR)
        self.assertEqual(len(result.module_plan.modules), 2)
        self.assertEqual(
            reconstruct_from_assembly_oligos(result.oligos),
            sequence,
        )

    def test_repeated_overlap_is_rejected_with_explainable_failure(self):
        sequence = "A" * 176

        with self.assertRaises(SYNOligoDesignFailure) as context:
            self.design(
                sequence,
                SYNOligoRules(
                    target_overlap_length=20,
                    min_overlap_length=20,
                    max_overlap_length=20,
                ),
            )

        self.assertTrue(
            any("不唯一" in reason for reason in context.exception.reasons),
        )


if __name__ == "__main__":
    unittest.main()
