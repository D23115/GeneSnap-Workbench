import unittest

from genesnap_workbench.sequence_core.shrna import build_shrna_oligo_pair
from genesnap_workbench.vector_library.shrna import (
    ShRNAVectorProtocolError,
    simulate_shrna_plasmid,
    validate_shrna_protocol,
)
from genesnap_workbench.vector_library.starters import load_public_plko1_puro_starter
from genesnap_workbench.vector_library.models import normalized_circular_checksum


PUBLIC_PLKO1_PURO_CHECKSUM = (
    "35b5699c39681b94a64da1a43c82e90b7e39df910596b55e4bc030b2e4c04a1c"
)


class ShRNAVectorTests(unittest.TestCase):
    def test_public_plko_starter_requires_confirmation_and_has_fixed_checksum(self):
        vector, protocol = load_public_plko1_puro_starter(user_confirmed=False)

        self.assertEqual(len(vector.sequence), 7050)
        self.assertEqual(vector.normalized_circular_sha256, PUBLIC_PLKO1_PURO_CHECKSUM)
        self.assertEqual(protocol.status, "pending_confirmation")
        self.assertEqual(protocol.experimental_validation_status, "unverified")
        self.assertEqual(protocol.left_site.name, "AgeI")
        self.assertEqual(protocol.right_site.name, "EcoRI")

    def test_confirmed_public_starter_validates_and_rebuilds_real_junctions(self):
        vector, protocol = load_public_plko1_puro_starter(user_confirmed=True)
        validation = validate_shrna_protocol(vector, protocol)
        oligos = build_shrna_oligo_pair(
            gene_symbol="TP53",
            target_no=1,
            target_id="TP53-target-1",
            target_sequence="GACTCCAGTGGTAATCTACTG",
        )

        simulation = simulate_shrna_plasmid(
            vector,
            protocol,
            target_id="TP53-target-1",
            forward_oligo=oligos.forward_sequence,
        )

        self.assertTrue(validation.is_valid)
        self.assertEqual(protocol.status, "enabled")
        self.assertIn(
            "A" + oligos.forward_sequence,
            simulation.expected_plasmid_sequence,
        )
        self.assertIn("GAATTC", simulation.expected_plasmid_sequence)
        self.assertIn(oligos.forward_sequence, simulation.expected_plasmid_sequence)
        self.assertEqual(simulation.left_site_count, 0)
        self.assertEqual(simulation.right_site_count, 1)
        self.assertEqual(
            simulation.expected_plasmid_checksum,
            normalized_circular_checksum(simulation.expected_plasmid_sequence),
        )

    def test_unconfirmed_protocol_cannot_generate_formal_simulation(self):
        vector, protocol = load_public_plko1_puro_starter(user_confirmed=False)
        oligos = build_shrna_oligo_pair(
            gene_symbol="TP53",
            target_no=1,
            target_id="TP53-target-1",
            target_sequence="GACTCCAGTGGTAATCTACTG",
        )

        with self.assertRaises(ShRNAVectorProtocolError):
            simulate_shrna_plasmid(
                vector,
                protocol,
                target_id="TP53-target-1",
                forward_oligo=oligos.forward_sequence,
            )


if __name__ == "__main__":
    unittest.main()
