"""
Protocol layer for TissueResolve.

Two conceptually distinct workflows live here.  They must not be confused.

Bulk protocol risk  (``protocol.metadata`` + ``protocol.risk``)
    Assesses gene-level biases introduced by mismatched sequencing protocols
    between bulk RNA-seq and single-cell references.  Used **before** bulk
    deconvolution to exclude or down-weight biased genes from the marker panel.

    Key classes: :class:`ProtocolMetadata`, :class:`ProtocolRiskAssessor`,
    :class:`ProtocolRiskReport`.

Spatial reference-query mismatch  (``protocol.mismatch``)
    Estimates per-gene multiplicative scale factors ``d_g`` that correct for
    the empirical expression difference between a Visium dataset and its
    pseudo-bulk reference.  Estimated **during** spatial deconvolution.

    Key classes: :class:`SpatialMismatch` (alias: :data:`ProtocolMismatch`),
    :func:`compute_spatial_discordance`, :func:`update_mismatch_factors`.
"""
# Metadata
from tissueresolve.protocol.metadata import (
    BulkProtocol,
    ProtocolMetadata,
    RefCapture,
    RefCounting,
    RefModality,
    SpatialPlatform,
)

# Bulk risk
from tissueresolve.protocol.risk import (
    ProtocolRiskAssessor,
    ProtocolRiskReport,
    RISK_HARD_THRESHOLD,
)

# Spatial mismatch
from tissueresolve.protocol.mismatch import (
    ProtocolMismatch,
    SpatialMismatch,
    compute_discordance,
    compute_spatial_discordance,
    update_mismatch_factors,
)

__all__ = [
    # metadata enums
    "BulkProtocol",
    "RefModality",
    "RefCapture",
    "RefCounting",
    "SpatialPlatform",
    # metadata container
    "ProtocolMetadata",
    # bulk risk
    "ProtocolRiskAssessor",
    "ProtocolRiskReport",
    "RISK_HARD_THRESHOLD",
    # spatial mismatch
    "SpatialMismatch",
    "ProtocolMismatch",
    "compute_spatial_discordance",
    "compute_discordance",
    "update_mismatch_factors",
]
