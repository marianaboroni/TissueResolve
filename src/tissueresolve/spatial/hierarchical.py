"""
Hierarchical (broad → fine) spatial / Visium deconvolution.

Mirrors :mod:`tissueresolve.bulk.hierarchical` for the spatial NB-CAR model
**without modifying the core solver**:

A. **Family level** — run the spatial model against a family-aggregated
   reference to estimate broad family composition per spot.
B. **Within family** — run the spatial model against the full fine reference,
   then renormalise to conditional ``P(subtype | family)`` per spot.
C. **Combine** — ``fine_spot(subtype) = family_spot(broad) × P(subtype | broad)``.
D. **Unresolved mass** — families whose subtypes are not separable keep their
   per-spot mass at the broad level (``unresolved_<family>``); the report and
   maps must show this rather than inventing per-subtype precision.

The CAR spatial-smoothing parameter ``lambda_spatial`` is recorded for both the
family-level and fine-level solves (never hidden).
"""
from __future__ import annotations

from dataclasses import dataclass, field
from pathlib import Path
from typing import Any, Optional

import numpy as np
import pandas as pd

from tissueresolve.reference.hierarchy import (
    HierarchicalEstimates,
    aggregate_reference_by_family,
    assemble_hierarchical_estimates,
)
from tissueresolve.results import QCReport, ReferenceSignature, SpatialDeconvResult

__all__ = [
    "HierarchicalSpatialResult",
    "run_hierarchical_spatial",
    "save_hierarchical_spatial_outputs",
]


@dataclass
class HierarchicalSpatialResult:
    """Result of a hierarchical spatial deconvolution run.

    ``deconv`` is a :class:`SpatialDeconvResult` whose ``proportions`` are the
    combined final per-spot estimates (resolved subtypes + ``unresolved_*``;
    rows sum to 1).  ``estimates`` holds the full hierarchical breakdown.
    """

    deconv: SpatialDeconvResult
    qc: QCReport
    estimates: HierarchicalEstimates
    family_result: Any
    fine_result: Any
    spot_qc: Optional[pd.DataFrame] = None
    morans_i: Optional[pd.Series] = None
    run_metadata: dict[str, Any] = field(default_factory=dict)


def run_hierarchical_spatial(
    Y,
    fine_ref: ReferenceSignature,
    array_row,
    array_col,
    lib_sizes,
    gene_names: list[str],
    mapping: dict[str, str],
    spot_ids: Optional[list[str]] = None,
    *,
    config=None,
    allow_unresolved: bool = True,
    unresolved_threshold: float = 0.10,
    min_discriminating_genes: int = 10,
    within_family_spillover_threshold: float = 0.30,
    allow_partial_resolution: bool = True,
    subtype_confidence_threshold: float = 0.10,
    hierarchical_gating: str = "soft",
    gating_version: str = "soft_gating-1.0",
    **run_kwargs,
) -> HierarchicalSpatialResult:
    """Run broad-to-fine hierarchical spatial deconvolution.

    Parameters mirror :meth:`SpatialPipeline.run`; *mapping* is the
    ``{fine_cell_type: broad_family}`` hierarchy.  Gating parameters control
    the within-family unresolved-mass decision.  ``**run_kwargs`` are forwarded
    to both spatial solves (``marker_genes``, ``run_neighbourhood`` …).
    """
    from tissueresolve.config import TissueResolveConfig
    from tissueresolve.spatial.pipeline import SpatialPipeline

    cfg = config or TissueResolveConfig()
    family_ref = aggregate_reference_by_family(fine_ref, mapping)

    pipeline = SpatialPipeline(cfg)
    # A. broad / family-level spatial deconvolution
    family_result = pipeline.run(
        Y, family_ref, array_row, array_col, lib_sizes, gene_names,
        spot_ids, **run_kwargs)
    # B. flat fine-level spatial deconvolution
    fine_result = SpatialPipeline(cfg).run(
        Y, fine_ref, array_row, array_col, lib_sizes, gene_names,
        spot_ids, **run_kwargs)

    estimates = assemble_hierarchical_estimates(
        family_result.deconv.proportions,
        fine_result.deconv.proportions,
        fine_ref,
        mapping,
        allow_unresolved=allow_unresolved,
        unresolved_threshold=unresolved_threshold,
        min_discriminating_genes=min_discriminating_genes,
        within_family_spillover_threshold=within_family_spillover_threshold,
        allow_partial_resolution=allow_partial_resolution,
        subtype_confidence_threshold=subtype_confidence_threshold,
        gating=hierarchical_gating,
        gating_version=gating_version,
        extra_metadata={
            "modality": "spatial",
            "lambda_spatial_family": float(family_result.deconv.lambda_spatial),
            "lambda_spatial_fine": float(fine_result.deconv.lambda_spatial),
        },
    )

    _trusted = estimates.metadata.get("trusted_resolution", {})
    _lam = float(family_result.deconv.lambda_spatial)
    run_metadata = {
        **family_result.run_metadata,
        "resolution_mode": "hierarchical",
        "merge_stage": "hierarchical",
        "hierarchy_mapping": dict(mapping),
        "broad_families": list(family_ref.cell_types),
        "fine_cell_types": list(fine_ref.cell_types),
        "lambda_spatial_family": _lam,
        "lambda_spatial_fine": float(fine_result.deconv.lambda_spatial),
        # modality provenance (spatial-specific; rule: explicit modality differences)
        "modality": "spatial",
        "prediction_unit": "spot",
        "gene_weighting_mode": "marker_unweighted_nb_car",
        "spatial_smoothing_used": bool(_lam > 0),
        "coordinates_used": True,
        "h_and_e_available": False,
        "fine_predictions_trusted": [f for f, s in _trusted.items()
                                     if s in ("full_fine", "selected_fine")],
        "fine_predictions_diagnostic_only": [f for f, s in _trusted.items()
                                             if s == "broad_only"],
        **estimates.metadata,
    }

    deconv = SpatialDeconvResult(
        proportions=estimates.combined_fine,
        cell_types=list(estimates.combined_fine.columns),
        marker_genes=family_result.deconv.marker_genes,
        n_iter=family_result.deconv.n_iter,
        converged=family_result.deconv.converged,
        convergence_trace=family_result.deconv.convergence_trace,
        lambda_spatial=family_result.deconv.lambda_spatial,
        mismatch_factors=family_result.deconv.mismatch_factors,
        lower_ci=None,
        upper_ci=None,
        bootstrap_coverage_note=family_result.deconv.bootstrap_coverage_note,
        n_smooth=family_result.deconv.n_smooth,
        run_metadata=run_metadata,
    )

    return HierarchicalSpatialResult(
        deconv=deconv,
        qc=family_result.qc,
        estimates=estimates,
        family_result=family_result,
        fine_result=fine_result,
        spot_qc=getattr(family_result, "spot_qc", None),
        morans_i=getattr(family_result, "morans_i", None),
        run_metadata=run_metadata,
    )


def save_hierarchical_spatial_outputs(
    result: HierarchicalSpatialResult, out_dir: Path | str,
) -> dict[str, Path]:
    """Write all hierarchical spatial TSVs under *out_dir*.

    Files: ``spatial_family_proportions.tsv`` /
    ``spatial_broad_proportions.tsv``,
    ``spatial_conditional_fine_proportions.tsv``,
    ``spatial_hierarchical_fine_proportions.tsv``,
    ``spatial_hierarchical_combined_proportions.tsv``,
    ``spatial_unresolved_family_mass.tsv``, ``spatial_hierarchical_qc.tsv``,
    ``cell_type_hierarchy.tsv``.
    """
    out = Path(out_dir)
    out.mkdir(parents=True, exist_ok=True)
    est = result.estimates
    written: dict[str, Path] = {}

    def _w(df: pd.DataFrame, name: str, comments: list[str] | None = None) -> None:
        path = out / name
        with path.open("w", encoding="utf-8") as fh:
            for c in comments or []:
                fh.write(f"# {c}\n")
            df.to_csv(fh, sep="\t")
        written[name] = path

    note = ["estimate_type: spot_rna_composition",
            f"lambda_spatial: {result.deconv.lambda_spatial}",
            "Hierarchical broad-to-fine estimates; see spatial_hierarchical_qc.tsv."]
    _w(est.family_proportions, "spatial_family_proportions.tsv", note)
    _w(est.family_proportions, "spatial_broad_proportions.tsv", note)
    _w(est.conditional_proportions, "spatial_conditional_fine_proportions.tsv",
       ["Conditional P(subtype | family) per spot; sums to 1 within family."])
    _w(est.fine_proportions, "spatial_hierarchical_fine_proportions.tsv",
       note + ["Resolved subtypes only (unresolved families are 0 here)."])
    _w(est.combined_fine, "spatial_hierarchical_combined_proportions.tsv",
       note + ["Resolved subtypes + unresolved_<family> columns; rows sum to 1."])
    _w(est.unresolved_mass, "spatial_unresolved_family_mass.tsv",
       ["Per-spot family mass not split into subtypes (not separable)."])

    qc_path = out / "spatial_hierarchical_qc.tsv"
    with qc_path.open("w", encoding="utf-8") as fh:
        fh.write("# Per-family within-family resolvability and decision.\n")
        est.qc.to_csv(fh, sep="\t", index=False)
    written["spatial_hierarchical_qc.tsv"] = qc_path

    from tissueresolve.reference.hierarchy import save_hierarchy_mapping
    written["cell_type_hierarchy.tsv"] = save_hierarchy_mapping(
        est.mapping, out / "cell_type_hierarchy.tsv")
    return written
