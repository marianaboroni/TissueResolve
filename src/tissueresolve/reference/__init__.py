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
from tissueresolve.reference.resolution import (
    CellTypeFamily,
    ResolutionConfig,
    ResolutionReport,
    apply_unresolved_mode,
    assign_resolution_families,
    build_nonseparable_graph,
    build_resolution_report,
    classify_resolvability,
    family_label_for_members,
    infer_cell_type_families,
    recommend_cell_type_merges,
    suggest_merges,
    summarize_resolution_report,
    write_recommended_merges,
)
from tissueresolve.reference.pairwise_markers import (
    augment_marker_panel_for_confusable_pairs,
    select_pairwise_discriminative_genes,
)
from tissueresolve.reference.hierarchy import (
    aggregate_predictions_by_family,
    aggregate_reference_by_group,
    build_hierarchical_reference,
    compare_fine_vs_merged_predictions,
    decompose_within_family,
    infer_broad_groups_from_labels,
    merge_reference_cell_types,
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
    # Stage 6 — resolution-aware layer
    "ResolutionConfig",
    "ResolutionReport",
    "CellTypeFamily",
    "build_resolution_report",
    "classify_resolvability",
    "infer_cell_type_families",
    "suggest_merges",
    "apply_unresolved_mode",
    "build_nonseparable_graph",
    "recommend_cell_type_merges",
    "assign_resolution_families",
    "family_label_for_members",
    "write_recommended_merges",
    "summarize_resolution_report",
    "select_pairwise_discriminative_genes",
    "augment_marker_panel_for_confusable_pairs",
    "infer_broad_groups_from_labels",
    "aggregate_reference_by_group",
    "merge_reference_cell_types",
    "build_hierarchical_reference",
    "decompose_within_family",
    "aggregate_predictions_by_family",
    "compare_fine_vs_merged_predictions",
]
