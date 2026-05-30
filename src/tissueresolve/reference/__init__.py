"""
Reference construction and gene selection for TissueResolve.

Public API
----------
``ReferenceBuilder``
    Build a :class:`~tissueresolve.results.ReferenceSignature` from a count
    matrix (DataFrame, AnnData, h5ad, or CSV/TSV).

``GeneFilterSet``
    Instance-based gene filter lists (intronic-dominant, dissociation-stress,
    etc.).  No module-level mutable singletons.

``GeneSelector``, ``MarkerSelectionResult``
    Composite marker gene selection combining CHIMERA-style scoring and
    SpatCAR-style pairwise discriminability.

``compute_separability``, ``merge_nonseparable_types``, ``SeparabilityWarning``
    Pairwise cell-type separability diagnostics (Bhattacharyya coefficient,
    Jeffreys divergence, Pearson r).  The merge function includes the P0-3
    dispersion bug fix from the SpatCAR legacy code.
"""
from tissueresolve.reference.build import ReferenceBuilder
from tissueresolve.reference.gene_filters import GeneFilterSet
from tissueresolve.reference.markers import GeneSelector, MarkerSelectionResult
from tissueresolve.reference.separability import (
    SeparabilityWarning,
    compute_separability,
    merge_nonseparable_types,
    separability_heatmap_data,
)

__all__ = [
    "ReferenceBuilder",
    "GeneFilterSet",
    "GeneSelector",
    "MarkerSelectionResult",
    "SeparabilityWarning",
    "compute_separability",
    "merge_nonseparable_types",
    "separability_heatmap_data",
]
