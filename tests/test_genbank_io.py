from __future__ import annotations

import builtins
import os
import tempfile
import unittest
from pathlib import Path
from unittest.mock import patch

from Bio.Seq import Seq
from Bio.SeqFeature import FeatureLocation, SeqFeature
from Bio.SeqRecord import SeqRecord

from genesnap_workbench.template_engine.genbank_io import write_genbank_utf8


class GenBankIOTests(unittest.TestCase):
    def test_unicode_record_is_written_when_system_text_encoding_is_cp1252(self):
        record = SeqRecord(Seq("ATGGCC"), id="TP53-delta", name="TP53-delta")
        record.annotations["molecule_type"] = "DNA"
        record.features.append(
            SeqFeature(
                FeatureLocation(0, len(record.seq)),
                type="misc_feature",
                qualifiers={"label": ["中文标签（Δ）"]},
            )
        )

        real_open = builtins.open

        def cp1252_default_open(file, mode="r", *args, **kwargs):
            if isinstance(file, (str, os.PathLike)) and "b" not in mode:
                kwargs.setdefault("encoding", "cp1252")
            return real_open(file, mode, *args, **kwargs)

        with tempfile.TemporaryDirectory() as temp_dir:
            output = Path(temp_dir) / "unicode.gb"
            with patch("Bio.File.open", side_effect=cp1252_default_open, create=True):
                write_genbank_utf8(record, output)

            payload = output.read_bytes()

        self.assertIn("中文标签（Δ）".encode("utf-8"), payload)


if __name__ == "__main__":
    unittest.main()
