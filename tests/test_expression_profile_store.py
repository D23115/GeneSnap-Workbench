import json
import tempfile
import unittest
from pathlib import Path

from genesnap_workbench.vector_library.expression_profiles import (
    ExpressionProfileIntegrityError,
    LocalExpressionProtocolStore,
)
from tests.test_expression_vector_protocol import vector_and_protocol


class ExpressionProfileStoreTests(unittest.TestCase):
    def setUp(self):
        self.temp_dir = tempfile.TemporaryDirectory()
        self.addCleanup(self.temp_dir.cleanup)
        self.store = LocalExpressionProtocolStore(Path(self.temp_dir.name))

    def test_saved_profile_reopens_with_exact_vector_and_protocol(self):
        vector, protocol = vector_and_protocol()

        saved = self.store.save_profile(vector, protocol)
        reopened = LocalExpressionProtocolStore(Path(self.temp_dir.name))
        loaded_vector, loaded_protocol = reopened.load_profile(saved.profile_id)

        self.assertEqual(loaded_vector, vector)
        self.assertEqual(loaded_protocol, protocol)
        self.assertEqual(reopened.list_profiles(), (saved,))
        self.assertEqual(saved.vector_checksum, vector.normalized_circular_sha256)
        self.assertEqual(saved.distribution_scope, "local_only")

    def test_tampered_local_sequence_is_blocked_on_load(self):
        vector, protocol = vector_and_protocol()
        saved = self.store.save_profile(vector, protocol)
        path = self.store.profile_path(saved.profile_id)
        payload = json.loads(path.read_text(encoding="utf-8"))
        payload["vector"]["sequence"] = payload["vector"]["sequence"][:-1] + "G"
        path.write_text(json.dumps(payload, ensure_ascii=False), encoding="utf-8")

        with self.assertRaises(ExpressionProfileIntegrityError):
            self.store.load_profile(saved.profile_id)


if __name__ == "__main__":
    unittest.main()
