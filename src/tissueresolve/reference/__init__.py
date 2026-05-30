"""
Reference construction and gene selection for TissueResolve.

Submodules (implemented in Stage 1):

``build``
    Unified ReferenceBuilder: donor-aware aggregation (from CHIMERA) +
    batched streaming (from SpatCAR) + NB overdispersion estimation.
    Produces ReferenceSignature with phi, R_cpm, phi_g, donor_cv.

``markers``
    GeneSelector: composite log-additive panel scoring (from CHIMERA) +
    pairwise discriminability augmentation (from SpatCAR v1.1).

``gene_filters``
    Gene list loading for protocol-aware filtering (intronic-dominant,
    dissociation-stress, length-biased, hypervariable, blacklist,
    protein-coding).  Instance-based — no module-level singletons.

``separability``
    compute_separability, merge_nonseparable_types, SeparabilityReport.
    Shared between bulk and spatial workflows.
"""
