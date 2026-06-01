"""
Multi-granularity gene panels: broad family / cell type / cell state.

Three panels are built, each with genes appropriate to *its* level:

* **broad** — genes that distinguish broad families (across all cells);
* **cell type within family** — genes that distinguish cell types inside each
  broad family;
* **state within cell type** — genes that distinguish states inside each cell
  type.

This reuses the committed within-family machinery
(:func:`tissueresolve.reference.within_family_markers.build_family_specific_gene_panels`,
which selects HVG ∪ pairwise-discriminative genes for the *children* of a
*parent* grouping) by calling it with the right (parent, child) column pair at
each level — so the scoring/weighting/exclusion logic is shared, not duplicated.
"""
from __future__ import annotations

from pathlib import Path
from typing import Mapping, Optional

import pandas as pd

from tissueresolve.reference.within_family_markers import (
    FamilyPanels,
    build_family_specific_gene_panels,
)

__all__ = [
    "build_broad_gene_panel",
    "build_celltype_within_family_gene_panels",
    "build_state_within_celltype_gene_panels",
    "compute_granularity_specific_gene_weights",
    "compute_pairwise_state_marker_support",
    "build_multigranularity_panels",
    "write_granular_signature_outputs",
]

_ROOT = "__tr_root__"


def build_broad_gene_panel(adata_ref, broad_col: str, **kw) -> list:
    """Genes that distinguish broad families (a single panel across all cells).

    Implemented by treating the whole dataset as one parent group whose children
    are the broad families.  Returns the gene list (also retrievable as the
    ``"ALL"`` entry of the underlying :class:`FamilyPanels`).
    """
    added = _ROOT not in adata_ref.obs.columns
    if added:
        adata_ref.obs[_ROOT] = "ALL"
    try:
        panels = build_family_specific_gene_panels(adata_ref, _ROOT, broad_col, **kw)
    finally:
        if added:
            del adata_ref.obs[_ROOT]
    return panels.family_panels.get("ALL", [])


def build_celltype_within_family_gene_panels(
    adata_ref, broad_col: str, cell_type_col: str, **kw
) -> FamilyPanels:
    """Per-broad-family panels of genes distinguishing the cell types within it."""
    return build_family_specific_gene_panels(adata_ref, broad_col, cell_type_col, **kw)


def build_state_within_celltype_gene_panels(
    adata_ref, cell_type_col: str, state_col: str, **kw
) -> FamilyPanels:
    """Per-cell-type panels of genes distinguishing the states within it."""
    return build_family_specific_gene_panels(adata_ref, cell_type_col, state_col, **kw)


def compute_pairwise_state_marker_support(
    state_panels: FamilyPanels,
) -> pd.DataFrame:
    """Pairwise state marker support (cell_type, state_a, state_b, n_markers)."""
    df = state_panels.marker_support_by_pair.copy()
    if df.empty:
        return pd.DataFrame(columns=["cell_type", "state_a", "state_b",
                                     "n_markers", "mean_score"])
    return df.rename(columns={"family": "cell_type", "type_a": "state_a",
                              "type_b": "state_b"})


def compute_granularity_specific_gene_weights(
    broad_panel: list,
    celltype_panels: FamilyPanels,
    state_panels: Optional[FamilyPanels] = None,
) -> pd.DataFrame:
    """Stack per-level gene weights into one tidy table.

    Columns: ``granularity`` (broad / cell_type / state), ``parent`` (the parent
    group: "ALL" for broad, broad family for cell_type, cell type for state),
    ``gene``, ``final_weight``.
    """
    rows = []
    for g in broad_panel:
        rows.append({"granularity": "broad", "parent": "ALL", "gene": g,
                     "final_weight": float("nan")})
    if not celltype_panels.gene_weights.empty:
        for _, r in celltype_panels.gene_weights.iterrows():
            rows.append({"granularity": "cell_type", "parent": r["family"],
                         "gene": r["gene"],
                         "final_weight": float(r["final_within_family_weight"])})
    if state_panels is not None and not state_panels.gene_weights.empty:
        for _, r in state_panels.gene_weights.iterrows():
            rows.append({"granularity": "state", "parent": r["family"],
                         "gene": r["gene"],
                         "final_weight": float(r["final_within_family_weight"])})
    return pd.DataFrame(rows, columns=["granularity", "parent", "gene",
                                       "final_weight"])


def build_multigranularity_panels(
    adata_ref, broad_col: str, cell_type_col: str,
    state_col: Optional[str] = None, **kw
) -> dict:
    """Build all available levels at once.

    Returns ``{"broad_panel", "celltype_panels", "state_panels",
    "gene_weights", "pairwise_state_markers"}``.  State level is skipped when
    ``state_col`` is None (returns ``state_panels=None``).
    """
    broad_panel = build_broad_gene_panel(adata_ref, broad_col, **kw)
    celltype_panels = build_celltype_within_family_gene_panels(
        adata_ref, broad_col, cell_type_col, **kw)
    state_panels = None
    pairwise_state = pd.DataFrame()
    if state_col and state_col in adata_ref.obs.columns:
        state_panels = build_state_within_celltype_gene_panels(
            adata_ref, cell_type_col, state_col, **kw)
        pairwise_state = compute_pairwise_state_marker_support(state_panels)
    weights = compute_granularity_specific_gene_weights(
        broad_panel, celltype_panels, state_panels)
    return {
        "broad_panel": broad_panel,
        "celltype_panels": celltype_panels,
        "state_panels": state_panels,
        "gene_weights": weights,
        "pairwise_state_markers": pairwise_state,
    }


def _panels_to_frame(panels: FamilyPanels, parent_name: str) -> pd.DataFrame:
    rows = []
    for parent, genes in (panels.family_panels or {}).items():
        for g in genes:
            rows.append({parent_name: parent, "gene": g})
    return pd.DataFrame(rows, columns=[parent_name, "gene"])


def write_granular_signature_outputs(result: dict, out_dir) -> dict[str, Path]:
    """Write the Part-3 signature tables; returns ``{name: path}``."""
    out = Path(out_dir)
    out.mkdir(parents=True, exist_ok=True)
    paths: dict[str, Path] = {}

    broad_df = pd.DataFrame({"gene": result["broad_panel"]})
    tables = {
        "broad_gene_panel": broad_df,
        "celltype_within_family_gene_panels":
            _panels_to_frame(result["celltype_panels"], "broad_family"),
        "granularity_gene_weights": result["gene_weights"],
    }
    if result.get("state_panels") is not None:
        tables["state_within_celltype_gene_panels"] = \
            _panels_to_frame(result["state_panels"], "cell_type")
        tables["pairwise_state_marker_support"] = result["pairwise_state_markers"]
    else:
        tables["state_within_celltype_gene_panels"] = pd.DataFrame(
            columns=["cell_type", "gene"])
        tables["pairwise_state_marker_support"] = pd.DataFrame(
            columns=["cell_type", "state_a", "state_b", "n_markers", "mean_score"])
    for name, df in tables.items():
        p = out / f"{name}.tsv"
        df.to_csv(p, sep="\t", index=False)
        paths[name] = p
    return paths
