import random
import unittest
from datetime import datetime, timezone

from genesnap_workbench.domain.syn import SYNDesignInput, SYNRoute
from genesnap_workbench.sequence_core.dna import sha256_sequence
from genesnap_workbench.sequence_core.syn_design import create_syn_design
from genesnap_workbench.vector_library.starters import (
    StarterVectorConfirmationRequired,
    load_public_puc57_starter,
)


NOW = datetime(2026, 7, 12, 13, 0, tzinfo=timezone.utc)


def artificial_sequence(length: int, seed: int) -> str:
    generator = random.Random(seed)
    return "".join(generator.choice("ACGT") for _ in range(length))


class SYNDesignEngineTests(unittest.TestCase):
    def test_public_puc57_requires_explicit_sequence_confirmation(self):
        with self.assertRaises(StarterVectorConfirmationRequired):
            load_public_puc57_starter("EcoRV", user_confirmed=False)

        vector, protocol = load_public_puc57_starter(
            "EcoRV",
            user_confirmed=True,
        )
        self.assertEqual(len(vector.sequence), 2710)
        self.assertEqual(
            vector.normalized_circular_sha256,
            "4c48a8f2eb55f2948a94270467df7c20572a82e971c0d716449bf22fa4dd6aba",
        )
        self.assertEqual(vector.sequence.count("GATATC"), 1)
        self.assertEqual(protocol.experimental_validation_status, "unverified")
        self.assertIn("snapgene.com", vector.public_source_urls[0])

    def test_low_risk_input_runs_through_complete_design_pipeline(self):
        sequence = artificial_sequence(600, seed=7)
        vector, protocol = load_public_puc57_starter(
            "SmaI",
            user_confirmed=True,
        )
        design = create_syn_design(
            SYNDesignInput(
                project_id="SYN-001",
                target_name="artificial-600",
                raw_sequence=f">artificial\n{sequence}",
                input_format="fasta",
                vector_protocol_version_id=protocol.protocol_version_id,
            ),
            vector,
            protocol,
            design_version_id="SYN-001-v1",
            created_at=NOW,
        )

        self.assertEqual(design.final_sequence, sequence)
        self.assertEqual(design.final_checksum, sha256_sequence(sequence))
        self.assertEqual(design.module_plan.route, SYNRoute.SINGLE_POOL)
        self.assertTrue(design.oligos)
        self.assertEqual(
            len(design.plasmid_simulation.expected_plasmid_sequence),
            len(vector.sequence) + len(sequence),
        )
        self.assertEqual(
            design.plasmid_simulation.design_version_id,
            design.design_version_id,
        )

    def test_1300_bp_input_uses_modular_route(self):
        sequence = artificial_sequence(1300, seed=19)
        vector, protocol = load_public_puc57_starter(
            "EcoRV",
            user_confirmed=True,
        )
        design = create_syn_design(
            SYNDesignInput(
                project_id="SYN-002",
                target_name="artificial-1300",
                raw_sequence=sequence,
                input_format="plain",
                vector_protocol_version_id=protocol.protocol_version_id,
            ),
            vector,
            protocol,
            design_version_id="SYN-002-v1",
            created_at=NOW,
        )

        self.assertEqual(design.module_plan.route, SYNRoute.MODULAR)
        self.assertGreaterEqual(len(design.module_plan.modules), 2)
        self.assertTrue(
            all(module.oligo_ids for module in design.module_plan.modules),
        )


if __name__ == "__main__":
    unittest.main()
