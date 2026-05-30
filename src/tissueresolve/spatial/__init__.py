"""
Spatial deconvolution workflow for 10x Visium (SpatCAR algorithm).

This subpackage contains all spatial-specific algorithmic components.
It must not import from ``bulk`` or mix spatial and bulk model assumptions.

Submodules (implemented in Stage 4):

``model``
    SpatCARModel — negative-binomial MAP with CAR spatial proximal mixing.
    Block coordinate descent: NB multiplicative update → spatial blend →
    mismatch factor update.

``graph``
    SpatialGraph, build_hex_graph, build_expression_weighted_graph.
    Hexagonal neighbourhood graph from Visium array coordinates.

``pipeline``
    SpatialPipeline — orchestrates: load Visium → build reference →
    select marker genes → separability check → build graph → fit model →
    QC → neighbourhood stats → report.

``qc``
    Spot-level QC: entropy, NB log-likelihood, spatial residual, dominant
    type.  Moran's I per cell type.  Model-level convergence QC.

``neighbourhood``
    NeighbourhoodStats — permutation-based co-occurrence testing with
    BH-FDR correction.  detect_spatial_niches.  distance_decay_cooccurrence.

``benchmark``
    Synthetic Visium simulation and four-method comparison (SpatCAR, DWLS,
    Spatial-NNLS, NNLS).
"""
