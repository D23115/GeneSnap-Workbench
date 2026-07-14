import unittest
from dataclasses import FrozenInstanceError

from genesnap_workbench.domain.workflows import (
    WorkflowDefinition,
    WorkflowRegistry,
    WorkflowTransition,
    build_default_workflow_registry,
)


def make_workflow(
    workflow_type: str,
    *,
    enabled: bool = True,
) -> WorkflowDefinition:
    return WorkflowDefinition(
        workflow_type=workflow_type,
        project_category="测试类",
        folder_suffix="TEST",
        enabled=enabled,
        intake_schema_key=f"{workflow_type}_intake_v1",
        intake_field_keys=("project_id", "target_name"),
        protocol_types=("test_protocol",),
        validator_key=f"{workflow_type}_validator",
        design_engine_adapter_key=f"{workflow_type}_engine",
        state_graph=(
            WorkflowTransition(
                from_state="recorded",
                action="complete",
                to_state="completed",
            ),
        ),
        artifact_types=("design_json",),
    )


class WorkflowRegistryTests(unittest.TestCase):
    def test_default_registry_exposes_enabled_syn_workflow(self):
        registry = build_default_workflow_registry()

        enabled = registry.list_enabled()

        self.assertEqual(
            tuple(definition.workflow_type for definition in enabled),
            ("de_novo_gene_synthesis", "shrna_knockdown"),
        )
        syn = registry.get("de_novo_gene_synthesis")
        self.assertEqual(syn.project_category, "合成/组装类")
        self.assertEqual(syn.folder_suffix, "SYN")
        self.assertEqual(syn.intake_schema_key, "syn_v0_intake")
        self.assertIn("target_name", syn.intake_field_keys)
        self.assertIn("raw_sequence", syn.intake_field_keys)
        self.assertEqual(
            syn.protocol_types,
            ("pUC57-EcoRV", "pUC57-SmaI", "custom_syn_vector"),
        )
        self.assertEqual(syn.validator_key, "syn_v0_validator")
        self.assertEqual(syn.design_engine_adapter_key, "syn_v0_design_engine")
        self.assertIn(
            WorkflowTransition(
                from_state="recorded",
                action="complete_design",
                to_state="design_completed",
            ),
            syn.state_graph,
        )
        shrna = registry.get("shrna_knockdown")
        self.assertEqual(shrna.project_category, "沉默/敲低类")
        self.assertEqual(shrna.folder_suffix, "KD")
        self.assertEqual(shrna.intake_schema_key, "shrna_v1_intake")
        self.assertEqual(shrna.design_engine_adapter_key, "shrna_v1_design_engine")
        self.assertIn("pLKO.1-AgeI-EcoRI", shrna.protocol_types)

    def test_duplicate_workflow_type_is_rejected(self):
        registry = WorkflowRegistry()
        registry.register(make_workflow("example"))

        with self.assertRaisesRegex(ValueError, "example"):
            registry.register(make_workflow("example"))

    def test_workflow_definition_is_immutable(self):
        definition = make_workflow("immutable_workflow")

        with self.assertRaises(FrozenInstanceError):
            definition.enabled = False

    def test_disabled_cko_workflow_is_not_listed_as_enabled(self):
        registry = WorkflowRegistry()
        registry.register(make_workflow("simple_cko_flox", enabled=False))

        self.assertEqual(registry.list_enabled(), ())
        self.assertEqual(
            registry.get("simple_cko_flox").workflow_type,
            "simple_cko_flox",
        )

    def test_new_workflow_can_be_registered_without_registry_changes(self):
        registry = WorkflowRegistry()
        definition = make_workflow("future_workflow")

        registry.register(definition)

        self.assertEqual(registry.get("future_workflow"), definition)
        self.assertEqual(registry.list_enabled(), (definition,))


if __name__ == "__main__":
    unittest.main()
