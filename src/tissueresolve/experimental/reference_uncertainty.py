"""Reference uncertainty & donor-aware signature reliability (P2c).

Makes the *reliability of the reference itself* explicit. A mean profile estimated
from 8 cells of 1 donor is treated by the solver exactly like one from 8000 cells of
20 donors — but it is far less trustworthy. This module quantifies that, per state
and per family, from the reference AnnData.

It is **diagnostic / reporting** first: the outputs are intended to feed gating and
reports (e.g. down-weight unreliable rare states, raise unresolved mass, explain why
a state is diagnostic-only). It does **not** change solver defaults.

Per-state quantities:
  n_cells, n_donors, donor_variability, mean_gene_cv, centroid_se, marker_stability,
  outlier_donor_score, reliability (∈[0,1], higher = more trustworthy).

Falls back safely when no donor column is available (cell-level variability only).
"""
from __future__ import annotations

from dataclasses import dataclass, field
from pathlib import Path
from typing import Mapping, Optional, Sequence

import numpy as np
import pandas as pd

__all__ = ["ReferenceUncertainty", "estimate_reference_uncertainty", "FEATURE_STATUS"]

FEATURE_STATUS = "experimental"
_EPS = 1e-9


@dataclass
class ReferenceUncertainty:
    per_state: pd.DataFrame
    family: pd.DataFrame = field(default_factory=pd.DataFrame)
    metadata: dict = field(default_factory=dict)

    def reliability(self) -> pd.Series:
        return self.per_state["reliability"]

    def write(self, out_dir) -> dict:
        out = Path(out_dir)
        out.mkdir(parents=True, exist_ok=True)
        paths = {}
        p1 = out / "reference_uncertainty.tsv"
        self.per_state.to_csv(p1, sep="\t")
        paths["reference_uncertainty"] = p1
        p2 = out / "state_reliability.tsv"
        self.per_state[["reliability"]].to_csv(p2, sep="\t")
        paths["state_reliability"] = p2
        if not self.family.empty:
            p3 = out / "family_reliability.tsv"
            self.family.to_csv(p3, sep="\t")
            paths["family_reliability"] = p3
        return paths


def _dense_lognorm(adata, genes: Optional[Sequence[str]]) -> tuple[np.ndarray, list]:
    """Return (cells × genes) log1p-CPM matrix and the gene list used."""
    import scipy.sparse as sp
    if genes is not None:
        keep = [g for g in genes if g in set(map(str, adata.var_names))]
        ad = adata[:, keep]
        gene_list = keep
    else:
        ad = adata
        gene_list = list(map(str, adata.var_names))
    X = ad.X
    X = X.toarray() if sp.issparse(X) else np.asarray(X)
    X = np.asarray(X, dtype=np.float64)
    lib = np.maximum(X.sum(axis=1, keepdims=True), 1.0)
    return np.log1p(X / lib * 1e4), gene_list


def estimate_reference_uncertainty(
    reference_adata,
    cell_type_col: str,
    donor_col: Optional[str] = None,
    genes: Optional[Sequence[str]] = None,
    mapping: Optional[Mapping[str, str]] = None,
    *,
    top_n_markers: int = 50,
) -> ReferenceUncertainty:
    """Per-state and per-family reference reliability from the reference AnnData.

    Parameters
    ----------
    reference_adata:
        AnnData with raw counts in ``X``, ``obs[cell_type_col]`` labels, optional
        ``obs[donor_col]``.
    mapping:
        Optional ``{state: family}`` for the family-level aggregation.
    """
    obs = reference_adata.obs
    if cell_type_col not in obs.columns:
        raise ValueError(f"cell_type_col {cell_type_col!r} not in obs.")
    has_donor = donor_col is not None and donor_col in obs.columns

    Xln, gene_list = _dense_lognorm(reference_adata, genes)
    labels = obs[cell_type_col].astype(str).to_numpy()
    donors = obs[donor_col].astype(str).to_numpy() if has_donor else None
    states = sorted(set(labels))

    rows = []
    for st in states:
        mask = labels == st
        n_cells = int(mask.sum())
        Xs = Xln[mask]                                  # (n_cells, G)
        centroid = Xs.mean(axis=0)
        # cell-level dispersion
        cell_sd = Xs.std(axis=0)
        mean_gene_cv = float(np.mean(cell_sd / (np.abs(centroid) + _EPS)))
        # centroid standard error (∝ 1/sqrt(n_cells)) — few cells ⇒ large SE
        centroid_se = float(np.median(cell_sd / np.sqrt(max(n_cells, 1))))
        # top markers by centroid expression
        top_idx = np.argsort(centroid)[::-1][:min(top_n_markers, centroid.size)]

        n_donors = 0
        donor_variability = float("nan")
        outlier_donor_score = float("nan")
        marker_stability = float("nan")
        if has_donor:
            d_here = donors[mask]
            uniq = sorted(set(d_here))
            n_donors = len(uniq)
            if n_donors >= 2:
                dmeans = np.vstack([Xs[d_here == d].mean(axis=0) for d in uniq])  # (D, G)
                # donor-to-donor variability: mean over genes of SD across donor means
                donor_sd = dmeans.std(axis=0)
                donor_variability = float(np.mean(donor_sd))
                # marker stability: 1 - mean CV of top markers across donors
                top_cv = donor_sd[top_idx] / (np.abs(centroid[top_idx]) + _EPS)
                marker_stability = float(np.clip(1.0 - np.mean(top_cv), 0.0, 1.0))
                # outlier donor: max distance of a donor mean to centroid / median distance
                dist = np.linalg.norm(dmeans - centroid[None, :], axis=1)
                med = float(np.median(dist)) + _EPS
                outlier_donor_score = float(np.max(dist) / med)
        if not has_donor or n_donors < 2:
            # cell-level fallback for stability
            top_cv = cell_sd[top_idx] / (np.abs(centroid[top_idx]) + _EPS)
            marker_stability = float(np.clip(1.0 - np.mean(top_cv), 0.0, 1.0))

        # reliability ∈ [0,1]: more cells, more donors, more stable ⇒ higher.
        cell_support = float(np.tanh(n_cells / 50.0))
        donor_support = float(np.tanh(n_donors / 3.0)) if has_donor else 0.5
        stability = marker_stability if np.isfinite(marker_stability) else 0.5
        var_penalty = (1.0 / (1.0 + donor_variability)) if np.isfinite(donor_variability) else 1.0
        reliability = float(np.clip(
            np.mean([cell_support, donor_support, stability, var_penalty]), 0.0, 1.0))

        rows.append({
            "state": st, "n_cells": n_cells, "n_donors": n_donors,
            "donor_variability": donor_variability, "mean_gene_cv": mean_gene_cv,
            "centroid_se": centroid_se, "marker_stability": marker_stability,
            "outlier_donor_score": outlier_donor_score,
            "reliability": reliability,
            "family": str(mapping.get(st, st)) if mapping else st,
        })

    per_state = pd.DataFrame(rows).set_index("state")

    family = pd.DataFrame()
    if mapping is not None:
        agg = per_state.groupby("family").agg(
            n_states=("reliability", "size"),
            mean_reliability=("reliability", "mean"),
            min_reliability=("reliability", "min"),
            total_cells=("n_cells", "sum"),
            mean_donor_variability=("donor_variability", "mean"))
        family = agg

    return ReferenceUncertainty(
        per_state=per_state, family=family,
        metadata={"feature_status": FEATURE_STATUS, "has_donor": bool(has_donor),
                  "n_states": len(states), "n_genes": len(gene_list),
                  "donor_fallback": not has_donor, "used_for": "diagnostics_and_gating"})
