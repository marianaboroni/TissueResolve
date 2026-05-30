"""
Uncertainty estimation for TissueResolve.

Submodules (implemented in Stages 3–4):

``bootstrap``
    BootstrapCI — gene-panel bootstrap for bulk deconvolution (from CHIMERA).
    bootstrap_proportions — parametric bootstrap for spatial deconvolution
    (from SpatCAR).  Both return (ci_lo, ci_hi, metadata) with documented
    nominal-vs-empirical coverage.

``stability``
    Cross-run stability checks (future).
"""
