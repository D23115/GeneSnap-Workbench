import unittest
from dataclasses import replace
from pathlib import Path

from Bio import SeqIO

from genesnap_workbench.vector_library.models import (
    RestrictionSite,
    SiteRetentionRule,
    SYNVectorProtocol,
    VectorRecord,
)
from genesnap_workbench.vector_library.syn import (
    SYNVectorProtocolError,
    circular_sequence_checksum,
    simulate_syn_plasmid,
    validate_syn_vector_protocol,
)


FIXTURE_PATH = (
    Path(__file__).parent
    / "fixtures"
    / "public"
    / "synthetic_puc57.gb"
)


def load_synthetic_vector() -> VectorRecord:
    record = SeqIO.read(FIXTURE_PATH, "genbank")
    return VectorRecord.from_sequence(
        vector_record_id="synthetic-puc57-v1",
        structural_display_name="人工 pUC57 测试骨架",
        sequence=str(record.seq),
        local_aliases=("synthetic-pUC57",),
        backbone_family="pUC57-like-test-only",
        public_equivalence_status="unknown",
    )


def make_protocol(
    vector: VectorRecord,
    *,
    site_name: str = "EcoRV",
    site_sequence: str = "GATATC",
    cut_offset: int = 3,
    site_retention_rule: SiteRetentionRule = SiteRetentionRule.NOT_REQUIRED,
    release_site: RestrictionSite | None = None,
) -> SYNVectorProtocol:
    return SYNVectorProtocol(
        protocol_id=f"synthetic-{site_name.lower()}",
        protocol_version_id=f"synthetic-{site_name.lower()}-v1",
        display_name=f"人工测试-{site_name}",
        status="enabled",
        experimental_validation_status="unverified",
        vector_record_id=vector.vector_record_id,
        vector_checksum=vector.normalized_circular_sha256,
        workflow_type="de_novo_gene_synthesis",
        insertion_mode="single_blunt_cut_hifi",
        linearization_site=RestrictionSite(
            name=site_name,
            sequence=site_sequence,
            cut_offset=cut_offset,
        ),
        site_retention_rule=site_retention_rule,
        release_site=release_site,
        homology_arm_length=20,
    )


class SYNVectorProtocolTests(unittest.TestCase):
    def test_fixture_is_explicitly_artificial_and_has_unique_default_sites(self):
        vector = load_synthetic_vector()

        self.assertIn("测试", vector.structural_display_name)
        self.assertEqual(vector.sequence.count("GATATC"), 1)
        self.assertEqual(vector.sequence.count("CCCGGG"), 1)

    def test_vector_checksum_mismatch_blocks_protocol(self):
        vector = load_synthetic_vector()
        protocol = replace(make_protocol(vector), vector_checksum="wrong-checksum")

        validation = validate_syn_vector_protocol(protocol, vector, insert_sequence="ACGT" * 20)

        self.assertFalse(validation.is_valid)
        self.assertIn("VECTOR_HASH_MISMATCH", validation.error_codes)

    def test_tampered_vector_record_checksum_is_detected(self):
        vector = load_synthetic_vector()
        tampered = replace(vector, sequence="A" + vector.sequence[1:])

        validation = validate_syn_vector_protocol(
            make_protocol(vector),
            tampered,
            insert_sequence="ACGT" * 20,
        )

        self.assertIn("VECTOR_RECORD_HASH_INVALID", validation.error_codes)

    def test_protocol_must_be_enabled_for_formal_design(self):
        vector = load_synthetic_vector()
        protocol = replace(make_protocol(vector), status="pending_confirmation")

        validation = validate_syn_vector_protocol(
            protocol,
            vector,
            insert_sequence="ACGT" * 20,
        )

        self.assertIn("PROTOCOL_NOT_ENABLED", validation.error_codes)

    def test_linearization_site_must_be_unique(self):
        vector = load_synthetic_vector()
        duplicated = VectorRecord.from_sequence(
            vector_record_id="duplicated-site",
            structural_display_name="重复位点测试载体",
            sequence=vector.sequence + "GATATC",
        )
        protocol = make_protocol(duplicated)

        validation = validate_syn_vector_protocol(
            protocol,
            duplicated,
            insert_sequence="ACGT" * 20,
        )

        self.assertIn("SITE_NOT_UNIQUE", validation.error_codes)

    def test_default_not_required_does_not_rebuild_linearization_site(self):
        vector = load_synthetic_vector()
        protocol = make_protocol(vector)
        insert = "ACGTCAGT" * 10

        simulation = simulate_syn_plasmid(
            vector,
            protocol,
            insert,
            design_version_id="design-v1",
        )

        self.assertEqual(simulation.design_version_id, "design-v1")
        self.assertEqual(
            len(simulation.expected_plasmid_sequence),
            len(vector.sequence) + len(insert),
        )
        self.assertNotIn("GATATC", simulation.expected_plasmid_sequence)
        self.assertEqual(simulation.site_retention_rule, "not_required")
        self.assertEqual(simulation.expected_digest_fragments_bp, ())
        insert_feature = next(
            feature for feature in simulation.features if feature.label == "SYN insert"
        )
        self.assertEqual(
            simulation.expected_plasmid_sequence[
                insert_feature.start:insert_feature.end
            ],
            insert,
        )

    def test_smai_protocol_uses_its_own_unique_cut(self):
        vector = load_synthetic_vector()
        protocol = make_protocol(
            vector,
            site_name="SmaI",
            site_sequence="CCCGGG",
            cut_offset=3,
        )

        simulation = simulate_syn_plasmid(
            vector,
            protocol,
            "ACGT" * 20,
            design_version_id="design-v1",
        )

        self.assertEqual(simulation.linearization_sites, ("SmaI",))
        self.assertNotIn("CCCGGG", simulation.expected_plasmid_sequence)

    def test_rebuild_flanking_sites_adds_two_release_sites_and_digest_plan(self):
        vector = load_synthetic_vector()
        release_site = RestrictionSite(
            name="EcoRI",
            sequence="GAATTC",
            cut_offset=1,
        )
        protocol = make_protocol(
            vector,
            site_retention_rule=SiteRetentionRule.REBUILD_FLANKING_SITES,
            release_site=release_site,
        )
        insert = "ACGTCAGT" * 10

        simulation = simulate_syn_plasmid(
            vector,
            protocol,
            insert,
            design_version_id="design-v1",
        )

        self.assertEqual(simulation.expected_plasmid_sequence.count("GAATTC"), 2)
        self.assertEqual(len(simulation.expected_digest_fragments_bp), 2)
        self.assertIn(
            len(insert) + len(release_site.sequence),
            simulation.expected_digest_fragments_bp,
        )

    def test_internal_release_site_blocks_complete_release_claim(self):
        vector = load_synthetic_vector()
        protocol = make_protocol(
            vector,
            site_retention_rule=SiteRetentionRule.REBUILD_FLANKING_SITES,
            release_site=RestrictionSite(
                name="EcoRI",
                sequence="GAATTC",
                cut_offset=1,
            ),
        )
        insert = "ACGTGAATTCTGCA"

        validation = validate_syn_vector_protocol(protocol, vector, insert)

        self.assertIn("RELEASE_SITE_INTERNAL", validation.error_codes)
        with self.assertRaises(SYNVectorProtocolError):
            simulate_syn_plasmid(
                vector,
                protocol,
                insert,
                design_version_id="design-v1",
            )

    def test_release_site_inside_backbone_blocks_two_fragment_digest_claim(self):
        original = load_synthetic_vector()
        vector = VectorRecord.from_sequence(
            vector_record_id="backbone-with-ecori",
            structural_display_name="带内部 EcoRI 的人工测试载体",
            sequence=original.sequence + "GAATTC",
        )
        protocol = make_protocol(
            vector,
            site_retention_rule=SiteRetentionRule.REBUILD_FLANKING_SITES,
            release_site=RestrictionSite(
                name="EcoRI",
                sequence="GAATTC",
                cut_offset=1,
            ),
        )

        validation = validate_syn_vector_protocol(
            protocol,
            vector,
            insert_sequence="ACGT" * 20,
        )

        self.assertIn("RELEASE_SITE_IN_BACKBONE", validation.error_codes)

    def test_final_plasmid_checksum_is_stable_across_circular_origin(self):
        vector = load_synthetic_vector()
        simulation = simulate_syn_plasmid(
            vector,
            make_protocol(vector),
            "ACGTCAGT" * 10,
            design_version_id="design-v1",
        )
        sequence = simulation.expected_plasmid_sequence
        rotated = sequence[37:] + sequence[:37]

        self.assertEqual(
            simulation.expected_plasmid_checksum,
            circular_sequence_checksum(rotated),
        )


if __name__ == "__main__":
    unittest.main()
