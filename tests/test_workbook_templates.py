import tempfile
import unittest
from pathlib import Path

from openpyxl import Workbook, load_workbook

from genesnap_workbench.template_engine.workbook_templates import (
    ContactProfile,
    LocalContactProfileStore,
    LocalWorkbookTemplateStore,
    inspect_workbook_template,
)


class WorkbookTemplateTests(unittest.TestCase):
    def setUp(self):
        self.temp_dir = tempfile.TemporaryDirectory()
        self.addCleanup(self.temp_dir.cleanup)
        self.root = Path(self.temp_dir.name)
        self.source = self.root / "vendor-primer.xlsx"
        workbook = Workbook()
        sheet = workbook.active
        sheet.title = "订购单"
        sheet["A1"] = "客户姓名"
        sheet["C1"] = "客户单位"
        sheet["A4"] = "引物名称"
        sheet["B4"] = "引物序列（5'-3'）"
        sheet["C4"] = "纯化方式"
        sheet["D4"] = "备注"
        sheet["A5"] = "示例行"
        workbook.save(self.source)

    def test_inspection_guesses_table_and_contact_mapping(self):
        inspected = inspect_workbook_template(self.source, kind="primer_order")

        self.assertEqual(inspected.sheet_name, "订购单")
        self.assertEqual(inspected.header_row, 4)
        self.assertEqual(inspected.data_start_row, 5)
        self.assertEqual(inspected.table_columns["primer_name"], 1)
        self.assertEqual(inspected.table_columns["sequence"], 2)
        self.assertEqual(inspected.table_columns["purification"], 3)
        self.assertEqual(inspected.contact_cells["customer_name"], "B1")
        self.assertEqual(inspected.contact_cells["organization"], "D1")

    def test_inspection_recognizes_tsingke_style_primer_headers_and_contacts(self):
        source = self.root / "tsingke-primer.xlsx"
        workbook = Workbook()
        sheet = workbook.active
        sheet.title = "引物合成订单"
        sheet["A5"] = " *客户姓名："
        sheet["A6"] = " *负责人姓名："
        sheet["A7"] = " *客户单位："
        sheet["A9"] = " *客户地址："
        sheet["A10"] = " *联系电话："
        sheet["A11"] = " *客户Email："
        headers = (
            "ID",
            "引物名称",
            "引物序列(5'to3')",
            "碱基数",
            "纯化方法",
            "OD/管",
            "管数",
            "5'修饰",
            "3'修饰",
            "Gene ID",
            "NM号",
            "扩增产物长度",
            "中间修饰",
            "备注",
        )
        for column, value in enumerate(headers, start=1):
            sheet.cell(12, column).value = value
        workbook.save(source)

        inspected = inspect_workbook_template(source, kind="primer_order")

        self.assertEqual(inspected.sheet_name, "引物合成订单")
        self.assertEqual(inspected.header_row, 12)
        self.assertEqual(
            inspected.table_columns,
            {
                "primer_name": 2,
                "sequence": 3,
                "length": 4,
                "purification": 5,
                "scale": 6,
                "note": 14,
            },
        )
        self.assertEqual(
            inspected.contact_cells,
            {
                "customer_name": "B5",
                "responsible_name": "B6",
                "organization": "B7",
                "address": "B9",
                "phone": "B10",
                "email": "B11",
            },
        )

    def test_saved_mapping_reopens_and_renders_without_reconfirmation(self):
        store = LocalWorkbookTemplateStore(self.root / "templates")
        inspected = inspect_workbook_template(self.source, kind="primer_order")
        saved = store.save_import(
            self.source,
            display_name="擎科引物订购表",
            inspected=inspected,
        )
        reopened = LocalWorkbookTemplateStore(self.root / "templates")

        self.assertEqual(reopened.list_profiles("primer_order"), (saved,))
        loaded = reopened.load_profile(saved.template_id)
        self.assertEqual(loaded.table_columns, inspected.table_columns)

        output = self.root / "filled.xlsx"
        reopened.render(
            saved.template_id,
            records=(
                {
                    "primer_name": "TP53-1-F",
                    "sequence": "ACGTACGT",
                    "purification": "PAGE",
                },
                {
                    "primer_name": "TP53-1-R",
                    "sequence": "TGCATGCA",
                },
            ),
            contact=ContactProfile(
                customer_name="示例用户",
                organization="示例单位",
            ),
            output_path=output,
        )

        workbook = load_workbook(output, data_only=False)
        sheet = workbook["订购单"]
        self.assertEqual(sheet["B1"].value, "示例用户")
        self.assertEqual(sheet["D1"].value, "示例单位")
        self.assertEqual(sheet["A5"].value, "TP53-1-F")
        self.assertEqual(sheet["B6"].value, "TGCATGCA")
        self.assertIsNone(sheet["C6"].value)

    def test_unmapped_optional_fields_do_not_block_rendering(self):
        store = LocalWorkbookTemplateStore(self.root / "templates")
        inspected = inspect_workbook_template(self.source, kind="primer_order")
        inspected.table_columns.pop("purification")
        saved = store.save_import(
            self.source,
            display_name="精简模板",
            inspected=inspected,
        )

        store.render(
            saved.template_id,
            records=({"primer_name": "P1", "sequence": "ACGT", "scale": "5 OD"},),
            contact=ContactProfile(customer_name="测试"),
            output_path=self.root / "optional.xlsx",
        )

        self.assertTrue((self.root / "optional.xlsx").exists())

    def test_single_default_contact_profile_can_be_changed(self):
        store = LocalContactProfileStore(self.root / "contact.json")
        self.assertEqual(store.load(), ContactProfile())

        store.save(ContactProfile(customer_name="示例用户", phone="00000000000"))

        reopened = LocalContactProfileStore(self.root / "contact.json")
        self.assertEqual(reopened.load().customer_name, "示例用户")
        self.assertEqual(reopened.load().phone, "00000000000")


if __name__ == "__main__":
    unittest.main()
