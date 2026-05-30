"""
Auto-generated figure captions for TissueResolve plots.

Captions encode the **estimate type** and **relevant parameters** so that a
figure is never separated from what it actually represents.  Two scientific
warnings are baked in by design:

* Bulk proportions are **RNA-derived mRNA proportions**, not absolute cell
  fractions.
* Spatial proportions are **spot-level RNA-derived composition estimates**,
  not direct single-cell counts.

These mirror the non-negotiable rules in ``CLAUDE.md`` and the estimate-type
constants on :class:`~tissueresolve.results.BulkDeconvResult` and
:class:`~tissueresolve.results.SpatialDeconvResult`.
"""
from __future__ import annotations

__all__ = [
    "BULK_PROP_WARNING",
    "SPATIAL_COMP_WARNING",
    "bulk_composition_caption",
    "bulk_qc_caption",
    "bulk_ci_caption",
    "spatial_abundance_caption",
    "dominant_type_caption",
    "spatial_qc_caption",
    "separability_heatmap_caption",
    "benchmark_estimated_vs_true_caption",
]

BULK_PROP_WARNING = (
    "Values are RNA-derived mRNA proportions, NOT absolute cell fractions; "
    "cell types differ in mRNA content per cell."
)
SPATIAL_COMP_WARNING = (
    "Values are spot-level RNA-derived cellular composition estimates, NOT "
    "direct single-cell counts."
)


def _join(parts: list[str]) -> str:
    return " ".join(p.strip() for p in parts if p and p.strip())


def bulk_composition_caption(
    n_samples: int, n_cell_types: int, *, gene_panel_size: int | None = None
) -> str:
    """Caption for the bulk cell-type composition plot."""
    panel = f" Deconvolution used {gene_panel_size} panel genes." if gene_panel_size else ""
    return _join([
        f"Bulk cell-type composition across {n_samples} sample(s) and "
        f"{n_cell_types} cell type(s).{panel}",
        BULK_PROP_WARNING,
    ])


def bulk_qc_caption(*, r2_warn: float | None = None) -> str:
    """Caption for the bulk QC (R² / profile correlation) plot."""
    thr = f" R² warning threshold = {r2_warn}." if r2_warn is not None else ""
    return _join([
        "Per-sample bulk QC: reconstruction R² and bulk-vs-reference profile "
        f"correlation.{thr} Thresholds are heuristic, not universal.",
        BULK_PROP_WARNING,
    ])


def bulk_ci_caption(sample: str, *, ci_level: float = 0.95, n_bootstrap: int | None = None) -> str:
    """Caption for the bulk bootstrap CI interval plot."""
    nb = f" from {n_bootstrap} bootstrap resamples" if n_bootstrap else ""
    return _join([
        f"Bootstrap {int(ci_level * 100)}% confidence intervals{nb} for sample "
        f"{sample!r}.",
        BULK_PROP_WARNING,
    ])


def spatial_abundance_caption(
    cell_type: str | list[str], *, lambda_spatial: float, n_spots: int
) -> str:
    """Caption for a spatial abundance map (single or multi-panel)."""
    if isinstance(cell_type, list):
        what = f"cell types {', '.join(cell_type)}"
    else:
        what = f"cell type {cell_type!r}"
    return _join([
        f"Spatial abundance map of {what} across {n_spots} spots "
        f"(λ_spatial = {lambda_spatial}; smoothing parameter recorded).",
        SPATIAL_COMP_WARNING,
    ])


def dominant_type_caption(*, lambda_spatial: float, n_spots: int) -> str:
    """Caption for the dominant cell-type map."""
    return _join([
        f"Dominant (argmax) cell type per spot across {n_spots} spots "
        f"(λ_spatial = {lambda_spatial}).",
        SPATIAL_COMP_WARNING,
        "The dominant label does not imply the spot is pure.",
    ])


def spatial_qc_caption(metric: str, *, lambda_spatial: float | None = None) -> str:
    """Caption for a spatial QC map."""
    lam = f" (λ_spatial = {lambda_spatial})" if lambda_spatial is not None else ""
    return _join([
        f"Per-spot spatial QC metric {metric!r}{lam}, plotted at the recorded "
        "array coordinates.",
        SPATIAL_COMP_WARNING,
    ])


def separability_heatmap_caption(
    *, n_critical: int = 0, n_high: int = 0, warn_threshold: float = 0.90
) -> str:
    """Caption for the pairwise separability heatmap."""
    flag = ""
    if n_critical or n_high:
        flag = (
            f" WARNING: {n_critical} CRITICAL and {n_high} HIGH poorly-separable "
            f"pair(s) (Bhattacharyya coefficient > {warn_threshold}); estimates "
            "for those types are unreliable."
        )
    return _join([
        "Pairwise cell-type separability (Bhattacharyya coefficient; 1 = "
        f"indistinguishable, 0 = orthogonal).{flag}",
    ])


def benchmark_estimated_vs_true_caption(method: str, scenario: str) -> str:
    """Caption for an estimated-vs-true benchmark scatter."""
    return _join([
        f"Estimated vs ground-truth proportions for method {method!r} on the "
        f"{scenario!r} synthetic scenario.  The dashed line is the identity "
        "(perfect recovery).  Synthetic benchmark only — not real-data "
        "validation.",
    ])
