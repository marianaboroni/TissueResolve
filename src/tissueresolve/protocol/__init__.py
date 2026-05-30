"""
Protocol metadata and mismatch handling for TissueResolve.

Two distinct protocol concepts are handled here — they must not be confused:

``metadata`` + ``risk``  (bulk)
    Protocol *risk* — systematic gene-level biases introduced by mismatched
    sequencing protocols (polyA bulk + snRNA reference, etc.).  Handled by
    excluding or down-weighting biased genes from the panel before solving.
    Uses curated gene lists in ``data/gene_lists/``.

``mismatch``  (spatial)
    Protocol *mismatch* — per-gene multiplicative scale factors ``d_g``
    estimated jointly with proportions during SpatCAR optimisation.  Accounts
    for the systematic expression difference between Visium and the sc/snRNA
    reference.

Submodules (implemented in Stage 2):

``metadata``
    BulkProtocol, RefModality, RefCapture, RefCounting enums.
    ProtocolMetadata dataclass.

``risk``
    ProtocolRiskAssessor, ProtocolRiskReport.
    Computes per-gene risk scores from gene lists.

``mismatch``
    ProtocolMismatch, compute_discordance, update_mismatch_factors.
    Spatial-only per-gene scale factor estimation.
"""
