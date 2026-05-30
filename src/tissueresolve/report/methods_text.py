"""
Auto-generated manuscript-style methods text for TissueResolve runs.

Each function returns a plain-text paragraph describing one analysis step,
with the actual parameter values filled in.  ``compose_bulk_methods`` and
``compose_spatial_methods`` assemble a full methods section from a pipeline
result.

Scientific-honesty rules (enforced by wording, asserted by tests):

* Bulk outputs are described as **RNA-derived mRNA proportions**, explicitly
  *not* absolute cell fractions unless mRNA-content correction was applied.
* Spatial outputs are described as **spot-level RNA-derived composition
  estimates**, explicitly *not* direct single-cell counts.
* The text does not overclaim: heuristic thresholds are named as heuristic,
  and synthetic benchmarks are named as synthetic.
"""
from __future__ import annotations

from typing import Any, Optional

__all__ = [
    "reference_construction_methods",
    "bulk_deconvolution_methods",
    "bootstrap_uncertainty_methods",
    "protocol_risk_methods",
    "spatial_graph_methods",
    "spatial_deconvolution_methods",
    "spatial_qc_methods",
    "separability_methods",
    "estimate_type_statement",
    "compose_bulk_methods",
    "compose_spatial_methods",
]


def estimate_type_statement(modality: str) -> str:
    """One-sentence statement of what the estimates represent."""
    if modality == "bulk":
        return (
            "Reported bulk values are RNA-derived mRNA proportions and are not "
            "equivalent to absolute cell fractions; conversion to cell fractions "
            "requires explicit mRNA-content correction with appropriate per-type "
            "mRNA content."
        )
    if modality == "spatial":
        return (
            "Reported spatial values are spot-level RNA-derived cellular "
            "composition estimates and do not represent direct single-cell "
            "counts."
        )
    raise ValueError(f"modality must be 'bulk' or 'spatial', got {modality!r}.")


def reference_construction_methods(
    *, n_genes: int, n_cell_types: int, genome: str = "hg38",
    donor_aware: bool = False, overdispersion: bool = False,
) -> str:
    donor = (
        "Cell-type profiles were computed as the mean of per-donor averages "
        "(donor-aware aggregation) and a cross-donor coefficient of variation "
        "was recorded per gene. "
        if donor_aware else
        "Cell-type profiles were computed as the mean expression across cells. "
    )
    od = (
        "Per-gene negative-binomial overdispersion was estimated for the "
        "spatial count model. "
        if overdispersion else ""
    )
    return (
        f"A reference signature was constructed from single-cell/single-nucleus "
        f"data ({genome}) over {n_cell_types} cell types and {n_genes} genes. "
        f"{donor}{od}"
        "The reference is stored both as an L1-normalised probability matrix "
        "(for weighted NNLS) and as a CPM matrix (for the negative-binomial "
        "spatial model)."
    )


def bulk_deconvolution_methods(*, n_genes_panel: int, condition_number: Optional[float] = None) -> str:
    kappa = (
        f" The panel condition number κ(Φ) was {condition_number:.1f}."
        if isinstance(condition_number, (int, float)) else ""
    )
    return (
        "Bulk cell-type composition was estimated per sample by weighted "
        "non-negative least squares (wNNLS) against a protocol-aware marker "
        f"panel of {n_genes_panel} genes, with proportions constrained to the "
        "simplex.{kappa} "
        "Outputs are RNA-derived mRNA proportions, not absolute cell fractions."
    ).replace("{kappa}", kappa)


def bootstrap_uncertainty_methods(
    *, n_bootstrap: int, ci_level: float = 0.95, coverage_note: Optional[str] = None
) -> str:
    note = f" {coverage_note}" if coverage_note else ""
    return (
        f"Uncertainty was quantified by gene-panel bootstrap ({n_bootstrap} "
        f"resamples), reporting {int(ci_level * 100)}% percentile confidence "
        f"intervals.{note}"
    )


def protocol_risk_methods(*, risk_level: str, n_excluded: int, mismatch_types: Optional[list] = None) -> str:
    mism = (
        f" Detected protocol mismatch(es): {', '.join(mismatch_types)}."
        if mismatch_types else ""
    )
    return (
        "Protocol compatibility between the bulk and reference assays was "
        f"assessed (overall risk: {risk_level}).{mism} "
        f"{n_excluded} gene(s) judged biased by the relevant protocol "
        "combination were excluded from the panel; affected-but-usable genes "
        "were down-weighted. No genes were removed silently — all excluded and "
        "down-weighted genes are listed in the output."
    )


def spatial_graph_methods(*, n_spots: int, mean_degree: Optional[float] = None) -> str:
    deg = f" (mean degree {mean_degree:.1f})" if isinstance(mean_degree, (int, float)) else ""
    return (
        f"A hexagonal neighbourhood graph was built from the Visium array "
        f"coordinates over {n_spots} spots{deg}, using a row-normalised "
        "adjacency matrix for the spatial prior."
    )


def spatial_deconvolution_methods(*, lambda_spatial: float, n_marker_genes: int,
                                  converged: bool, n_iter: int) -> str:
    conv = "converged" if converged else "did NOT converge (estimates may be suboptimal)"
    return (
        "Spatial cell-type composition was estimated with a negative-binomial "
        "count model and a conditional-autoregressive (CAR) spatial prior "
        f"(λ_spatial = {lambda_spatial}; the smoothing strength is recorded and "
        "never applied silently) over a panel of "
        f"{n_marker_genes} marker genes. Block coordinate descent {conv} after "
        f"{n_iter} iteration(s). "
        "Outputs are spot-level RNA-derived composition estimates, not single-cell counts."
    )


def spatial_qc_methods(*, has_morans_i: bool = False) -> str:
    mi = " Moran's I spatial autocorrelation was computed per cell type." if has_morans_i else ""
    return (
        "Per-spot quality control included proportion entropy, the "
        "negative-binomial log-likelihood, the dominant cell type, and the "
        f"spatial residual.{mi}"
    )


def separability_methods(*, n_critical: int = 0, n_high: int = 0, warn_threshold: float = 0.90) -> str:
    flag = ""
    if n_critical or n_high:
        flag = (
            f" {n_critical} CRITICAL and {n_high} HIGH poorly-separable pair(s) "
            "were flagged; estimates for those cell types are reported with a "
            "reliability warning and were not presented as confident."
        )
    return (
        "Pairwise cell-type separability was assessed with the Bhattacharyya "
        f"coefficient (poorly-separable threshold {warn_threshold}), the Jeffreys "
        f"divergence, and profile correlation.{flag}"
    )


def compose_bulk_methods(result: Any) -> str:
    """Assemble a full bulk methods section from a ``BulkPipelineResult``."""
    deconv = result.deconv
    qc = result.qc
    meta = getattr(result, "run_metadata", {}) or {}
    paras: list[str] = []
    paras.append(reference_construction_methods(
        n_genes=meta.get("n_reference_genes", len(deconv.gene_panel)),
        n_cell_types=deconv.n_cell_types,
        genome=meta.get("genome", "hg38"),
        donor_aware=bool(meta.get("donor_aware", False)),
    ))
    risk = getattr(result, "protocol_risk", None)
    if risk is not None:
        paras.append(protocol_risk_methods(
            risk_level=getattr(risk, "risk_level", "unknown"),
            n_excluded=getattr(risk, "n_genes_excluded", 0),
            mismatch_types=getattr(risk, "mismatch_types", None),
        ))
    paras.append(bulk_deconvolution_methods(
        n_genes_panel=len(deconv.gene_panel),
        condition_number=getattr(qc, "condition_number", None),
    ))
    if deconv.lower_ci is not None:
        paras.append(bootstrap_uncertainty_methods(
            n_bootstrap=int(meta.get("n_bootstrap", 0)),
            ci_level=float(meta.get("ci_level", 0.95)),
        ))
    paras.append(estimate_type_statement("bulk"))
    return "\n\n".join(paras)


def compose_spatial_methods(result: Any, *, separability: Any = None) -> str:
    """Assemble a full spatial methods section from a ``SpatialPipelineResult``."""
    deconv = result.deconv
    meta = getattr(result, "run_metadata", {}) or {}
    paras: list[str] = []
    paras.append(reference_construction_methods(
        n_genes=meta.get("n_reference_genes", len(deconv.marker_genes)),
        n_cell_types=deconv.n_cell_types,
        genome=meta.get("genome", "hg38"),
        overdispersion=True,
    ))
    paras.append(spatial_graph_methods(n_spots=deconv.n_spots))
    paras.append(spatial_deconvolution_methods(
        lambda_spatial=deconv.lambda_spatial,
        n_marker_genes=len(deconv.marker_genes),
        converged=deconv.converged, n_iter=deconv.n_iter,
    ))
    paras.append(spatial_qc_methods(
        has_morans_i=getattr(result, "morans_i", None) is not None
    ))
    if separability is not None:
        paras.append(separability_methods(
            n_critical=getattr(separability, "n_critical", 0),
            n_high=getattr(separability, "n_high", 0),
        ))
    paras.append(estimate_type_statement("spatial"))
    return "\n\n".join(paras)
