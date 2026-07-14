"""Biological sequence operations and design primitives."""

from .dna import (
    InvalidDNASequenceError,
    normalize_dna,
    reverse_complement,
    sha256_sequence,
)
from .syn_qc import SYNQCRules, evaluate_syn_sequence
from .syn_modules import SYNModuleRules, plan_syn_modules
from .syn_oligos import (
    BiopythonThermodynamicAnalyzer,
    SYNOligoDesignFailure,
    SYNOligoRules,
    design_assembly_oligos,
    reconstruct_from_assembly_oligos,
)
from .expression import create_expression_design
from .reporter import create_reporter_design


_SHRNA_EXPORTS = {
    "advance_blast_selection",
    "build_shrna_oligo_pair",
    "create_shrna_design",
    "select_initial_candidates",
}


def __getattr__(name: str):
    if name in _SHRNA_EXPORTS:
        from . import shrna

        return getattr(shrna, name)
    raise AttributeError(name)

__all__ = [
    "InvalidDNASequenceError",
    "BiopythonThermodynamicAnalyzer",
    "SYNModuleRules",
    "SYNOligoDesignFailure",
    "SYNOligoRules",
    "SYNQCRules",
    "design_assembly_oligos",
    "build_shrna_oligo_pair",
    "advance_blast_selection",
    "create_shrna_design",
    "create_expression_design",
    "create_reporter_design",
    "evaluate_syn_sequence",
    "normalize_dna",
    "plan_syn_modules",
    "reconstruct_from_assembly_oligos",
    "reverse_complement",
    "select_initial_candidates",
    "sha256_sequence",
]
