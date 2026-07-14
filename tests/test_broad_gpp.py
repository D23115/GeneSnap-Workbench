import unittest

from genesnap_workbench.integrations.broad_gpp import (
    parse_hairpin_candidates,
    parse_oligo_detail,
)


RESULT_HTML = """
<html><body>
<table class="grid">
  <tr><th></th><th>Start Pos.</th><th>Intrinsic Score</th>
      <th>Target Sequence</th><th>Oligo Design</th><th>Existing Clone</th></tr>
  <tr class="z1"><td>1</td><td>31</td><td>15.000</td>
      <td>GTTTCGGGTATCCGGTTAAAT</td>
      <td><a href="/gpp/public/oligo/design?seq=GTTTCGGGTATCCGGTTAAAT">oligo</a></td><td></td></tr>
  <tr class="z2"><td>2</td><td>210</td><td>13.200</td>
      <td>GAATCGGGTATCCCGTTAAAT</td>
      <td><a href="/gpp/public/oligo/design?seq=GAATCGGGTATCCCGTTAAAT">oligo</a></td><td></td></tr>
</table>
</body></html>
"""


OLIGO_HTML = """
<html><body>
<h2>Forward Oligo:</h2>
<p>5'-CCGGGTTTCGGGTATCCGGTTAAATCTCGAGATTTAACCGGATACCCGAAACTTTTTG-3'</p>
<h2>Reverse Oligo:</h2>
<p>5'-AATTCAAAAAGTTTCGGGTATCCGGTTAAATCTCGAGATTTAACCGGATACCCGAAAC-3'</p>
</body></html>
"""


class BroadGPPParserTests(unittest.TestCase):
    def test_parses_ranked_candidates_and_oligo_links(self):
        candidates = parse_hairpin_candidates(RESULT_HTML)

        self.assertEqual(len(candidates), 2)
        self.assertEqual(candidates[0].source_rank, 1)
        self.assertEqual(candidates[0].start_position, 31)
        self.assertEqual(str(candidates[0].intrinsic_score), "15.000")
        self.assertEqual(candidates[0].target_sequence, "GTTTCGGGTATCCGGTTAAAT")
        self.assertIn("oligo/design", candidates[0].oligo_url)

    def test_parses_broad_full_oligo_sequences(self):
        oligos = parse_oligo_detail(OLIGO_HTML)

        self.assertEqual(
            oligos.forward_sequence,
            "CCGGGTTTCGGGTATCCGGTTAAATCTCGAGATTTAACCGGATACCCGAAACTTTTTG",
        )
        self.assertEqual(
            oligos.reverse_sequence,
            "AATTCAAAAAGTTTCGGGTATCCGGTTAAATCTCGAGATTTAACCGGATACCCGAAAC",
        )


if __name__ == "__main__":
    unittest.main()
