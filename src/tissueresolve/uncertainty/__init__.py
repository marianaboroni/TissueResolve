"""
Uncertainty estimation for TissueResolve.

Public API (Stage 3)
--------------------
``bootstrap.BulkBootstrapCI``
    Gene-panel bootstrap for bulk deconvolution.  Returns empirical
    percentile CIs for each sample × cell-type.  Non-suppressible warning
    emitted when fewer than ``MIN_PANEL_GENES`` genes are available.

Future (Stage 4)
----------------
``bootstrap.SpatialBootstrapCI``
    Parametric bootstrap for spatial deconvolution (from SpatCAR).

``stability``
    Cross-run stability checks.
"""
from tissueresolve.uncertainty.bootstrap import BulkBootstrapCI, SpatialBootstrapCI

__all__ = ["BulkBootstrapCI", "SpatialBootstrapCI"]
