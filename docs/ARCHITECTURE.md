# Architecture Notes

This clean implementation keeps the main subsystems separate:

- `domain`: stable data models.
- `sequence_core`: sequence IO, normalization, coordinate mapping, primer/insert primitives.
- `vector_library`: user-imported vectors, insertion sites, and construction protocols.
- `sequencing`: AB1/SEQ matching and PASS/FAIL/WARNING judgments.
- `template_engine`: xlsx/docx export using user-provided templates.
- `project_workflow`: project intake, deadlines, status transitions, and folder structure.
- `integrations`: NCBI, Ensembl, Broad, and optional private/local integrations.
- `app`: PySide6 UI.

The legacy `project-review-package-v2` source is a reference implementation for two protocol prototypes:

- pLKO/shRNA.
- LV-037/OE.

Private reference materials are outside the package and must not become default public fixtures.
