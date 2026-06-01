"""
State-aware three-level hierarchical bulk deconvolution (broad → cell type →
state), inspired by BayesPrism's state↔type aggregation — but built entirely
from TissueResolve's existing, tested two-level machinery (no Bayesian Gibbs
sampling, no new core solver).

Strategy
--------
1. Aggregate the state-level reference up to cell-type and broad references.
2. Deconvolve at broad and cell-type granularity (existing ``BulkPipeline``).
3. **Stage 1** — ``assemble_hierarchical_estimates`` with the cell-type→broad
   mapping and the *cell-type-within-family* gene panels → cell-type proportions
   with ``unresolved_<broad_family>`` where a family is not type-separable.
4. **Stage 2** — ``assemble_hierarchical_estimates`` again, parent = the resolved
   cell-type proportions, child = state deconvolution, with the
   *state-within-cell-type* panels → state proportions with
   ``unresolved_<cell_type>`` where a cell type is not state-separable.
5. Re-attach the broad-level unresolved mass so the final state vector sums to 1.

Unresolved mass is therefore assigned at the **correct level** and never hidden;
states are estimated only where the panels make them separable (the gating lives
in ``assemble_hierarchical_estimates``).  No prediction-value semantics of the
existing two-level path are changed.
"""
from __future__ import annotations

from dataclasses import dataclass, field
from pathlib import Path
from typing import Any, Optional

import pandas as pd

from tissueresolve.reference.hierarchy import (
    aggregate_reference_by_family,
    assemble_hierarchical_estimates,
)
from tissueresolve.reference.three_level_hierarchy import ThreeLevelHierarchy
from tissueresolve.results import ReferenceSignature

__all__ = ["StateAwareBulkResult", "run_state_aware_hierarchical_bulk",
           "write_state_aware_outputs"]


@dataclass
class StateAwareBulkResult:
    broad_proportions: pd.DataFrame
    cell_type_proportions: pd.DataFrame          # cell types + unresolved_<broad>
    state_proportions: Optional[pd.DataFrame]    # states + unresolved_<cell_type> + unresolved_<broad>
    unresolved_by_level: pd.DataFrame            # tidy: level, label, mean_mass
    allocation_metadata: pd.DataFrame
    celltype_estimates: Any = None               # stage-1 HierarchicalEstimates
    state_estimates: Any = None                  # stage-2 HierarchicalEstimates
    metadata: dict = field(default_factory=dict)


def _unresolved_cols(df: pd.DataFrame) -> list[str]:
    return [c for c in df.columns if str(c).startswith("unresolved_")]


def run_state_aware_hierarchical_bulk(
    bulk: pd.DataFrame,
    state_ref: ReferenceSignature,
    hierarchy: ThreeLevelHierarchy,
    *,
    config=None,
    celltype_panels: Optional[dict] = None,
    state_panels: Optional[dict] = None,
    allow_unresolved: bool = True,
    unresolved_threshold: float = 0.10,
    min_discriminating_genes: int = 10,
    within_family_spillover_threshold: float = 0.30,
    allow_partial_resolution: bool = True,
    subtype_confidence_threshold: float = 0.10,
    **run_kwargs,
) -> StateAwareBulkResult:
    """Run broad → cell type → state hierarchical bulk deconvolution.

    *state_ref* is the finest (state-level) :class:`ReferenceSignature`
    (``cell_types`` = states; if ``hierarchy.has_states`` is False they are cell
    types and the state stage is skipped).  *celltype_panels* /
    *state_panels* are ``{group: [genes]}`` dicts from
    :func:`granular_signatures.build_multigranularity_panels`.
    """
    from tissueresolve.config import TissueResolveConfig
    from tissueresolve.bulk.pipeline import BulkPipeline

    cfg = config or TissueResolveConfig()
    gating = dict(allow_unresolved=allow_unresolved,
                  unresolved_threshold=unresolved_threshold,
                  min_discriminating_genes=min_discriminating_genes,
                  within_family_spillover_threshold=within_family_spillover_threshold,
                  allow_partial_resolution=allow_partial_resolution,
                  subtype_confidence_threshold=subtype_confidence_threshold)

    # --- references at each granularity (aggregate upward) ---
    if hierarchy.has_states:
        celltype_ref = aggregate_reference_by_family(
            state_ref, hierarchy.state_to_celltype)
    else:
        celltype_ref = state_ref
    broad_ref = aggregate_reference_by_family(
        celltype_ref, hierarchy.celltype_to_broad)

    # --- flat deconvolutions ---
    broad_props = BulkPipeline(cfg).run(bulk, broad_ref, **run_kwargs).deconv.proportions
    celltype_props = BulkPipeline(cfg).run(bulk, celltype_ref, **run_kwargs).deconv.proportions

    # --- Stage 1: cell types within broad families (unresolved_<broad>) ---
    ct_est = assemble_hierarchical_estimates(
        broad_props, celltype_props, celltype_ref,
        dict(hierarchy.celltype_to_broad),
        family_gene_panels=celltype_panels,
        extra_metadata={"level": "cell_type", "within_family_panels":
                        bool(celltype_panels)}, **gating)

    cell_type_proportions = ct_est.combined_fine
    broad_unresolved = ct_est.unresolved_mass            # unresolved_<broad>

    state_proportions = None
    st_est = None
    if hierarchy.has_states:
        state_props = BulkPipeline(cfg).run(bulk, state_ref, **run_kwargs).deconv.proportions
        # parent for the state split = resolved cell-type columns (exclude
        # unresolved_<broad>, which stays at the broad level)
        ct_cols = [c for c in cell_type_proportions.columns
                   if not str(c).startswith("unresolved_")]
        parent_for_state = cell_type_proportions[ct_cols]
        st_est = assemble_hierarchical_estimates(
            parent_for_state, state_props, state_ref,
            dict(hierarchy.state_to_celltype),
            family_gene_panels=state_panels,
            extra_metadata={"level": "state", "within_cell_type_panels":
                            bool(state_panels)}, **gating)
        # final state vector = states + unresolved_<cell_type> + re-attached
        # unresolved_<broad> (so the row sums back to 1)
        state_proportions = pd.concat([st_est.combined_fine, broad_unresolved], axis=1)

    # --- tidy unresolved-by-level summary ---
    rows = []
    for c in _unresolved_cols(broad_unresolved):
        rows.append({"level": "broad_family", "label": c,
                     "mean_mass": float(broad_unresolved[c].mean())})
    if st_est is not None:
        for c in _unresolved_cols(st_est.unresolved_mass):
            rows.append({"level": "cell_type", "label": c,
                         "mean_mass": float(st_est.unresolved_mass[c].mean())})
    unresolved_by_level = pd.DataFrame(rows, columns=["level", "label", "mean_mass"])

    # --- allocation metadata (which level resolved, panel sizes) ---
    meta_rows = []
    ct_res = ct_est.resolvability
    for fam in ct_res.index:
        n_panel = int(ct_res.loc[fam, "n_panel_genes"]) if "n_panel_genes" in ct_res.columns else -1
        meta_rows.append({
            "level": "cell_type_within_broad", "group": str(fam),
            "resolved": bool(ct_res.loc[fam, "resolvable"]),
            "n_panel_genes": n_panel,
            "panel_source": "cell_type_within_family" if celltype_panels else "global",
            "reason": str(ct_res.loc[fam, "reason"]),
        })
    if st_est is not None:
        st_res = st_est.resolvability
        for ct in st_res.index:
            n_panel = int(st_res.loc[ct, "n_panel_genes"]) if "n_panel_genes" in st_res.columns else -1
            meta_rows.append({
                "level": "state_within_cell_type", "group": str(ct),
                "resolved": bool(st_res.loc[ct, "resolvable"]),
                "n_panel_genes": n_panel,
                "panel_source": "state_within_cell_type" if state_panels else "global",
                "reason": str(st_res.loc[ct, "reason"]),
            })
    allocation_metadata = pd.DataFrame(meta_rows)

    return StateAwareBulkResult(
        broad_proportions=ct_est.family_proportions,
        cell_type_proportions=cell_type_proportions,
        state_proportions=state_proportions,
        unresolved_by_level=unresolved_by_level,
        allocation_metadata=allocation_metadata,
        celltype_estimates=ct_est,
        state_estimates=st_est,
        metadata={"has_states": hierarchy.has_states,
                  "used_celltype_panels": bool(celltype_panels),
                  "used_state_panels": bool(state_panels)},
    )


def write_state_aware_outputs(result: StateAwareBulkResult, out_dir) -> dict[str, Path]:
    """Write the Part-4 deconvolution tables; returns ``{name: path}``."""
    out = Path(out_dir)
    out.mkdir(parents=True, exist_ok=True)
    tables = {
        "broad_proportions": result.broad_proportions,
        "cell_type_proportions": result.cell_type_proportions,
        "unresolved_mass_by_level": result.unresolved_by_level,
        "hierarchical_allocation_metadata": result.allocation_metadata,
    }
    if result.state_proportions is not None:
        tables["state_proportions"] = result.state_proportions
    paths = {}
    for name, df in tables.items():
        p = out / f"{name}.tsv"
        df.to_csv(p, sep="\t")
        paths[name] = p
    return paths
