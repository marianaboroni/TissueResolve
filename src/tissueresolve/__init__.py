"""
TissueResolve — unified cell-type and cell-state deconvolution.

Workflows
---------
``tissueresolve bulk``    — bulk RNA-seq deconvolution (CHIMERA algorithm).
``tissueresolve spatial`` — 10x Visium spatial deconvolution (SpatCAR algorithm).

Scientific output convention
-----------------------------
Bulk:
    ``BulkDeconvResult.proportions`` contains **mRNA proportions**,
    not cell fractions.  The estimate type is recorded in every output file.

Spatial:
    ``SpatialDeconvResult.proportions`` contains **spot-level RNA-derived
    cellular composition estimates**, not direct single-cell counts unless
    explicitly calibrated via mRNA content correction.
"""
from importlib.metadata import PackageNotFoundError, version as _version

try:
    __version__: str = _version("tissueresolve")
except PackageNotFoundError:
    __version__ = "0.1.0-dev"

from tissueresolve.api import (
    build_reference,
    deconv_bulk,
    deconv_spatial,
    generate_report,
    plot_results,
)

__all__ = [
    "__version__",
    "build_reference",
    "deconv_bulk",
    "deconv_spatial",
    "generate_report",
    "plot_results",
]
