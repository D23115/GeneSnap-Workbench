import unittest

import genesnap_workbench


class SmokeTests(unittest.TestCase):
    def test_package_imports(self):
        self.assertEqual(genesnap_workbench.__version__, "0.3.4")

    def test_vector_library_can_be_imported_before_design_engine(self):
        from genesnap_workbench.vector_library.models import VectorRecord
        from genesnap_workbench.sequence_core.syn_design import create_syn_design

        self.assertIsNotNone(VectorRecord)
        self.assertTrue(callable(create_syn_design))


if __name__ == "__main__":
    unittest.main()
