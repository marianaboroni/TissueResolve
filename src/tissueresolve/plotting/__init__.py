"""
Publication-quality plotting for TissueResolve.

Every plot function saves the underlying data alongside the figure.
Figures without saved data violate the DESIGN_SPEC.

Submodules (implemented in Stage 5):

``style``
    Shared matplotlib rcParams, colour palettes, figure factory.

``bulk_plots``
    Proportion heatmaps, R² scatter, CI bar charts.

``spatial_plots``
    Spatial proportion maps, convergence traces, mismatch factor plots,
    niche maps.

``qc_plots``
    Separability heatmap, marker recall bar charts, Moran's I plots,
    protocol risk summary.

``benchmark_plots``
    Method comparison figures, ablation tables.

``captions``
    Auto-generated figure captions with full parameter provenance.
"""
