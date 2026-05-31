"""
Hierarchical (broad → fine) bulk RNA-seq deconvolution.

This module orchestrates the two-stage broad-to-fine strategy **without
modifying the core wNNLS solver**:

A. **Family level** — deconvolve the bulk against a family-aggregated
   reference to estimate broad cell-type-family mRNA proportions.
B. **Within family** — deconvolve the bulk against the full fine reference to
   estimate subtype proportions, which are renormalised to conditional
   ``P(subtype | family)`` proportions within each family.
C. **Combine** — ``fine(subtype) = family(broad) × P(subtype | family)``.
D. **Unresolved mass** — families whose subtypes are not separable (low
   within-family separability, too few discriminating genes, or high
   within-family spillover) keep their mass at the broad level as
   ``unresolved_<family>`` rather than being split into unsupported subtypes.

Both stages reuse :class:`~tissueresolve.bulk.pipeline.BulkPipeline`; the
hierarchy arithmetic and gating live in
:mod:`tissueresolve.reference.hierarchy`.
"""
from __future__ import annotations

from dataclasses import dataclass, field
from pathlib import Path
from typing import Any, Optional

import pandas as pd

from tissueresolve.reference.hierarchy import (
    HierarchicalEstimates,
    aggregate_reference_by_family,
    assemble_hierarchical_estimates,
)
from tissueresolve.results import BulkDeconvResult, QCReport, ReferenceSignature

__all__ = [
    "HierarchicalBulkResult",
    "run_hierarchical_bulk",
    "save_hierarchical_bulk_outputs",
]


@dataclass
class HierarchicalBulkResult:
    """Result of a hierarchical bulk deconvolution run.

    Attributes
    ----------
    deconv:
        :class:`BulkDeconvResult` whose ``proportions`` are the **combined
        final** estimates (resolved subtypes + ``unresolved_<family>`` columns;
        rows sum to 1).  Provided so the report/plotting/CLI code that expects a
        standard pipeline result keeps working.
    qc:
        Family-level :class:`QCReport` (from the broad-level solve).
    estimates:
        Full :class:`HierarchicalEstimates` (family / conditional / fine /
        unresolved / per-family resolvability).
    family_result, fine_result:
        The underlying broad-level and flat-fine-level pipeline results.
    run_metadata:
        Merged run metadata, including the hierarchy mapping and gating
        thresholds.
    """

    deconv: BulkDeconvResult
    qc: QCReport
    estimates: HierarchicalEstimates
    family_result: Any
    fine_result: Any
    gene_selection: Optional[Any] = None
    protocol_risk: Optional[Any] = None
    run_metadata: dict[str, Any] = field(default_factory=dict)


def run_hierarchical_bulk(
    bulk: pd.DataFrame,
    fine_ref: ReferenceSignature,
    mapping: dict[str, str],
    *,
    config=None,
    allow_unresolved: bool = True,
    unresolved_threshold: float = 0.10,
    min_discriminating_genes: int = 10,
    within_family_spillover_threshold: float = 0.30,
    allow_partial_resolution: bool = True,
    subtype_confidence_threshold: float = 0.10,
    **run_kwargs,
) -> HierarchicalBulkResult:
    """Run broad-to-fine hierarchical bulk deconvolution.

    Parameters
    ----------
    bulk:
        Bulk counts (genes × samples or samples × genes — the pipeline
        auto-orients).
    fine_ref:
        Fine-level :class:`ReferenceSignature` (subtypes as ``cell_types``).
    mapping:
        ``{fine_cell_type: broad_family}`` covering every reference cell type.
    config:
        Optional :class:`~tissueresolve.config.TissueResolveConfig`.
    allow_unresolved, unresolved_threshold, min_discriminating_genes,
    within_family_spillover_threshold:
        Within-family resolvability gating parameters.
    **run_kwargs:
        Forwarded verbatim to both :meth:`BulkPipeline.run` calls
        (``gene_panel``, ``n_bootstrap``, ``protocol_meta`` …).

    Returns
    -------
    HierarchicalBulkResult
    """
    from tissueresolve.config import TissueResolveConfig
    from tissueresolve.bulk.pipeline import BulkPipeline

    cfg = config or TissueResolveConfig()

    family_ref = aggregate_reference_by_family(fine_ref, mapping)

    # A. broad / family-level deconvolution
    family_result = BulkPipeline(cfg).run(bulk, family_ref, **run_kwargs)
    # B. flat fine-level deconvolution (used for conditional within-family split)
    fine_result = BulkPipeline(cfg).run(bulk, fine_ref, **run_kwargs)

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
        extra_metadata={"modality": "bulk"},
    )

    run_metadata = {
        **family_result.run_metadata,
        "resolution_mode": "hierarchical",
        "merge_stage": "hierarchical",
        "hierarchy_mapping": dict(mapping),
        "broad_families": list(family_ref.cell_types),
        "fine_cell_types": list(fine_ref.cell_types),
        **estimates.metadata,
    }

    deconv = BulkDeconvResult(
        proportions=estimates.combined_fine,
        coverage_r2=family_result.deconv.coverage_r2,
        gene_panel=family_result.deconv.gene_panel,
        gene_weights=family_result.deconv.gene_weights,
        lower_ci=None,
        upper_ci=None,
        cell_fractions=None,
        run_metadata=run_metadata,
    )

    return HierarchicalBulkResult(
        deconv=deconv,
        qc=family_result.qc,
        estimates=estimates,
        family_result=family_result,
        fine_result=fine_result,
        gene_selection=getattr(family_result, "gene_selection", None),
        protocol_risk=getattr(family_result, "protocol_risk", None),
        run_metadata=run_metadata,
    )


def save_hierarchical_bulk_outputs(
    result: HierarchicalBulkResult, out_dir: Path | str,
) -> dict[str, Path]:
    """Write all hierarchical bulk TSVs under *out_dir*.

    Files (per the spec):

    * ``bulk_family_proportions.tsv`` / ``bulk_broad_proportions.tsv``
    * ``bulk_conditional_fine_proportions.tsv``
    * ``bulk_hierarchical_fine_proportions.tsv``  (resolved subtypes only)
    * ``bulk_hierarchical_combined_proportions.tsv``  (resolved + unresolved)
    * ``bulk_unresolved_family_mass.tsv``
    * ``bulk_hierarchical_qc.tsv``
    * ``cell_type_hierarchy.tsv``
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

    note = ["estimate_type: mRNA_proportion",
            "Hierarchical broad-to-fine estimates; see bulk_hierarchical_qc.tsv."]
    _w(est.family_proportions, "bulk_family_proportions.tsv", note)
    # alias name used in the revised spec
    _w(est.family_proportions, "bulk_broad_proportions.tsv", note)
    _w(est.conditional_proportions, "bulk_conditional_fine_proportions.tsv",
       ["Conditional P(subtype | family); sums to 1 within each family."])
    _w(est.fine_proportions, "bulk_hierarchical_fine_proportions.tsv",
       note + ["Resolved subtypes only (unresolved families are 0 here)."])
    _w(est.combined_fine, "bulk_hierarchical_combined_proportions.tsv",
       note + ["Resolved subtypes + unresolved_<family> columns; rows sum to 1."])
    _w(est.unresolved_mass, "bulk_unresolved_family_mass.tsv",
       ["Family mass not split into subtypes (subtypes not separable)."])

    qc_path = out / "bulk_hierarchical_qc.tsv"
    with qc_path.open("w", encoding="utf-8") as fh:
        fh.write("# Per-family within-family resolvability and decision.\n")
        est.qc.to_csv(fh, sep="\t", index=False)
    written["bulk_hierarchical_qc.tsv"] = qc_path

    from tissueresolve.reference.hierarchy import save_hierarchy_mapping
    written["cell_type_hierarchy.tsv"] = save_hierarchy_mapping(
        est.mapping, out / "cell_type_hierarchy.tsv")
    return written
