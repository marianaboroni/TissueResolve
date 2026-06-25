"""
Hyper-parameter tuning (NOT IMPLEMENTED).

This module previously contained a stub that wrote **fabricated** tuning grids
and scores (e.g. ``marker_genes=50 -> 0.8``, ``lambda=0.05 -> 0.9``) to disk
without running any real search.  Emitting invented metrics violates the
project's scientific-integrity rules, so the stub has been quarantined: the
functions now raise :class:`NotImplementedError` instead of writing fake
artifacts.

When real tuning is implemented it must:

* search an explicit, recorded grid,
* score candidates with a real, ground-truth-free criterion
  (e.g. gene-masking cross-validation for bulk; reconstruction vs.
  over-smoothing for spatial — see ``tissueresolve.validation.gene_masking``
  and ``tissueresolve.spatial.auto_params``),
* and persist the true scores, not placeholders.

Until then, use the existing data-driven selectors:

* bulk solver backbone — ``solver="auto"`` (gene-masking CV), and
* spatial smoothing — ``tissueresolve.spatial.auto_params.select_lambda_spatial``.
"""
from __future__ import annotations

from pathlib import Path
from typing import Any, Dict, Optional

__all__ = ["tune_bulk_parameters", "tune_spatial_parameters"]

_NOT_IMPLEMENTED = (
    "tissueresolve.tuning.{name} is not implemented. The previous stub wrote "
    "fabricated tuning metrics and has been removed. Use solver='auto' "
    "(gene-masking cross-validation) for bulk backbone selection, or "
    "tissueresolve.spatial.auto_params.select_lambda_spatial for spatial "
    "smoothing selection."
)


def tune_bulk_parameters(
    output_dir: Path, grid: Optional[Dict[str, list]] = None
) -> Dict[str, Any]:
    raise NotImplementedError(_NOT_IMPLEMENTED.format(name="tune_bulk_parameters"))


def tune_spatial_parameters(
    output_dir: Path, grid: Optional[Dict[str, list]] = None
) -> Dict[str, Any]:
    raise NotImplementedError(_NOT_IMPLEMENTED.format(name="tune_spatial_parameters"))
