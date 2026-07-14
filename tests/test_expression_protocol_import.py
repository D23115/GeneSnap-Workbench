import unittest

from genesnap_workbench.sequence_core.dna import reverse_complement
from genesnap_workbench.vector_library import scan_restriction_sites
from genesnap_workbench.vector_library.expression_import import (
    resolve_manual_homology,
    resolve_restriction_insertion,
    scan_expression_restriction_sites,
)
from genesnap_workbench.vector_library.models import ExpressionVectorProtocol


# LV-037 已确认插入区两侧的脱敏局部结构：NheI 位于 Kozak 前，BamHI 位于右侧侧翼。
LV037_LOCAL_INSERT_REGION = (
    "AAAAAAAAAAAAAAAAAAAAAAAAAAAAAA"
    "GAAGCTTGGGCTGCAGGTCGACGCTAGCGCCACC"
    "TTTTTTTTTTTTTTTTTTTT"
    "GAATTCCGAGGATCCATGGAC"
    "AAAAAAAAAAAAAAAAAAAAAAAAAAAAAA"
)

SCORING_VECTOR = (
    "AAAAAAAAAAAAAAAAAAAAAAAAAAAAAA"
    "ATGCATGCATGCATGCATG"
    "GCTAGC"
    "TTTTTTTTTTTTTTTTTTTT"
    "GGATCCATGCATGCATGCATG"
    "AAAAAAAAAAAAAAAAAAAAAAAAAAAAAA"
)


def _occurrence(sequence: str, enzyme_name: str):
    return next(
        item
        for item in scan_restriction_sites(sequence)
        if item.enzyme_name == enzyme_name
    )


class ExpressionProtocolImportTests(unittest.TestCase):
    def test_scans_lv037_nhei_and_bamhi_at_unique_stable_positions(self):
        occurrences = scan_restriction_sites(LV037_LOCAL_INSERT_REGION)
        nhei = _occurrence(LV037_LOCAL_INSERT_REGION, "NheI")
        bamhi = _occurrence(LV037_LOCAL_INSERT_REGION, "BamHI")

        self.assertEqual(
            sum(item.enzyme_name == "NheI" for item in occurrences),
            1,
        )
        self.assertEqual(
            sum(item.enzyme_name == "BamHI" for item in occurrences),
            1,
        )
        self.assertEqual(nhei.start, LV037_LOCAL_INSERT_REGION.index("GCTAGC"))
        self.assertEqual(nhei.cut_position, nhei.start + 1)
        self.assertEqual(bamhi.start, LV037_LOCAL_INSERT_REGION.index("GGATCC"))
        self.assertEqual(bamhi.cut_position, bamhi.start + 1)
        self.assertEqual(
            tuple((item.start, item.enzyme_name) for item in occurrences),
            tuple(sorted((item.start, item.enzyme_name) for item in occurrences)),
        )
        self.assertEqual(
            scan_expression_restriction_sites(LV037_LOCAL_INSERT_REGION),
            occurrences,
        )

    def test_scans_circular_site_across_vector_origin_once_with_normalized_cut(self):
        circular_sequence = "CTAGC" + "A" * 30 + "G"

        nhei_occurrences = tuple(
            item
            for item in scan_restriction_sites(circular_sequence)
            if item.enzyme_name == "NheI"
        )

        self.assertEqual(len(nhei_occurrences), 1)
        self.assertEqual(nhei_occurrences[0].start, len(circular_sequence) - 1)
        self.assertEqual(nhei_occurrences[0].cut_position, 0)

    def test_restriction_resolution_rejects_circular_origin_spanning_site_in_chinese(self):
        circular_sequence = "CTAGC" + "A" * 30 + "G"
        nhei = _occurrence(circular_sequence, "NheI")

        with self.assertRaisesRegex(ValueError, "跨越载体首尾"):
            resolve_restriction_insertion(
                circular_sequence,
                left_occurrence=nhei,
                right_occurrence=nhei,
            )

    def test_restriction_resolution_selects_deterministic_arms_in_primer_directions(self):
        nhei = _occurrence(LV037_LOCAL_INSERT_REGION, "NheI")
        bamhi = _occurrence(LV037_LOCAL_INSERT_REGION, "BamHI")

        resolution = resolve_restriction_insertion(
            LV037_LOCAL_INSERT_REGION,
            left_occurrence=nhei,
            right_occurrence=bamhi,
        )

        self.assertEqual(resolution.source, "restriction_sites")
        self.assertEqual(resolution.left_boundary, nhei.cut_position)
        self.assertEqual(resolution.right_boundary, bamhi.cut_position)
        self.assertEqual(
            resolution.left_primer_homology,
            LV037_LOCAL_INSERT_REGION[
                nhei.cut_position - len(resolution.left_primer_homology) : nhei.cut_position
            ],
        )
        self.assertEqual(
            reverse_complement(resolution.right_primer_homology),
            LV037_LOCAL_INSERT_REGION[
                bamhi.cut_position : bamhi.cut_position
                + len(resolution.right_primer_homology)
            ],
        )
        self.assertTrue(18 <= len(resolution.left_primer_homology) <= 25)
        self.assertTrue(18 <= len(resolution.right_primer_homology) <= 25)
        self.assertEqual(resolution.left_site.enzyme_name, "NheI")
        self.assertEqual(resolution.right_site.enzyme_name, "BamHI")

    def test_resolvers_accept_planned_positional_arguments(self):
        nhei = _occurrence(LV037_LOCAL_INSERT_REGION, "NheI")
        bamhi = _occurrence(LV037_LOCAL_INSERT_REGION, "BamHI")
        forward_homology = LV037_LOCAL_INSERT_REGION[nhei.cut_position - 20 : nhei.cut_position]
        right_top_strand = LV037_LOCAL_INSERT_REGION[
            bamhi.cut_position : bamhi.cut_position + 20
        ]

        restriction_resolution = resolve_restriction_insertion(
            LV037_LOCAL_INSERT_REGION,
            nhei,
            bamhi,
        )
        manual_resolution = resolve_manual_homology(
            LV037_LOCAL_INSERT_REGION,
            forward_homology,
            reverse_complement(right_top_strand),
        )

        self.assertEqual(
            (restriction_resolution.left_boundary, restriction_resolution.right_boundary),
            (manual_resolution.left_boundary, manual_resolution.right_boundary),
        )

    def test_homology_selection_prioritizes_gc_tm_and_twenty_base_length(self):
        nhei = _occurrence(SCORING_VECTOR, "NheI")
        bamhi = _occurrence(SCORING_VECTOR, "BamHI")

        resolution = resolve_restriction_insertion(
            SCORING_VECTOR,
            left_occurrence=nhei,
            right_occurrence=bamhi,
        )

        self.assertEqual(resolution.left_primer_homology, "ATGCATGCATGCATGCATGG")
        self.assertEqual(resolution.left_homology_tm, 60)
        self.assertEqual(
            reverse_complement(resolution.right_primer_homology),
            "GATCCATGCATGCATGCATG",
        )

    def test_single_restriction_site_resolves_equal_boundaries_and_protocol_accepts_them(self):
        nhei = _occurrence(LV037_LOCAL_INSERT_REGION, "NheI")

        resolution = resolve_restriction_insertion(
            LV037_LOCAL_INSERT_REGION,
            left_occurrence=nhei,
            right_occurrence=nhei,
        )

        self.assertEqual(resolution.left_boundary, resolution.right_boundary)
        protocol = ExpressionVectorProtocol(
            protocol_id="single-site-test",
            protocol_version_id="single-site-test-v1",
            display_name="single-site test",
            status="enabled",
            experimental_validation_status="unverified",
            vector_record_id="lv037-local",
            vector_checksum="checksum",
            workflow_type="expression",
            insertion_mode="single_restriction_site",
            left_boundary=resolution.left_boundary,
            right_boundary=resolution.right_boundary,
            left_primer_homology=resolution.left_primer_homology,
            right_primer_homology=resolution.right_primer_homology,
            kozak_sequence="GCCACC",
            stop_codon_rule="preserve",
            c_terminal_fusion_name=None,
        )
        self.assertEqual(protocol.left_boundary, protocol.right_boundary)

    def test_manual_homology_resolves_boundaries_without_coordinates(self):
        nhei = _occurrence(LV037_LOCAL_INSERT_REGION, "NheI")
        bamhi = _occurrence(LV037_LOCAL_INSERT_REGION, "BamHI")
        forward_homology = LV037_LOCAL_INSERT_REGION[nhei.cut_position - 20 : nhei.cut_position]
        right_top_strand = LV037_LOCAL_INSERT_REGION[
            bamhi.cut_position : bamhi.cut_position + 20
        ]

        resolution = resolve_manual_homology(
            LV037_LOCAL_INSERT_REGION,
            forward_primer_homology=forward_homology,
            reverse_primer_homology=reverse_complement(right_top_strand),
        )

        self.assertEqual(resolution.source, "manual_homology")
        self.assertEqual(resolution.left_boundary, nhei.cut_position)
        self.assertEqual(resolution.right_boundary, bamhi.cut_position)
        self.assertEqual(resolution.left_primer_homology, forward_homology)
        self.assertEqual(
            resolution.right_primer_homology,
            reverse_complement(right_top_strand),
        )
        self.assertIsNone(resolution.left_site)
        self.assertIsNone(resolution.right_site)

    def test_manual_homology_rejects_multiple_left_matches_even_with_one_ordered_pair(self):
        left_homology = "GATTACAGATTACAGATTAC"
        right_top_strand = "CCTGACCTGACCTGACCTGA"
        vector_sequence = (
            left_homology
            + "TTTT"
            + right_top_strand
            + "CCCC"
            + left_homology
        )

        with self.assertRaisesRegex(ValueError, "F 引物.*多个匹配"):
            resolve_manual_homology(
                vector_sequence,
                forward_primer_homology=left_homology,
                reverse_primer_homology=reverse_complement(right_top_strand),
            )

    def test_manual_homology_rejects_multiple_right_matches_even_with_one_ordered_pair(self):
        left_homology = "GATTACAGATTACAGATTAC"
        right_top_strand = "CCTGACCTGACCTGACCTGA"
        vector_sequence = (
            right_top_strand
            + "TTTT"
            + left_homology
            + "CCCC"
            + right_top_strand
        )

        with self.assertRaisesRegex(ValueError, "R 引物.*多个匹配"):
            resolve_manual_homology(
                vector_sequence,
                forward_primer_homology=left_homology,
                reverse_primer_homology=reverse_complement(right_top_strand),
            )


if __name__ == "__main__":
    unittest.main()
