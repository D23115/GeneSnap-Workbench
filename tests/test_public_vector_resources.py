from __future__ import annotations

import hashlib
import unittest
from pathlib import Path

from Bio import SeqIO


ROOT = Path(__file__).resolve().parents[1]
VECTOR_ROOT = ROOT / "src" / "genesnap_workbench" / "resources" / "vectors"
EXPECTED_SEQUENCE_HASHES = {
    "plko1_puro_snapgene_public.gb": (
        "13459c0789caffafc35b6f79f285da19c240a4c563c87d2da14951455d1e013e"
    ),
    "puc57_snapgene_public.gb": (
        "3ee8c04a82fa801b4ce9b17b434208301c1b59091f7425299930f2ea96d7734f"
    ),
}


class PublicVectorResourceTests(unittest.TestCase):
    def test_bundled_references_keep_sequence_but_not_third_party_annotations(self):
        for filename, expected_hash in EXPECTED_SEQUENCE_HASHES.items():
            with self.subTest(filename=filename):
                record = SeqIO.read(VECTOR_ROOT / filename, "genbank")
                sequence_hash = hashlib.sha256(
                    str(record.seq).upper().encode("ascii")
                ).hexdigest()

                self.assertEqual(sequence_hash, expected_hash)
                self.assertEqual([feature.type for feature in record.features], ["source"])
                self.assertIn("Source:", record.annotations.get("comment", ""))


if __name__ == "__main__":
    unittest.main()
