"""
Spillover (cross-type leakage) simulation and estimation.

Even with a perfect algorithm, signal from one cell type leaks into the
estimate of a similar cell type.  This module *measures* that leakage so it
can be reported honestly:

* simulate pure / pairwise / realistic mixtures from a reference,
* deconvolve them with the real bulk solver and build a **spillover matrix**
  (row = true type, column = predicted type; rows sum to 1),
* summarise per-type spillover risk and the main leaking partner.

For the spatial workflow, where running the full NB-CAR model per simulated
spot is expensive, :func:`expression_spillover_proxy` gives a cheap
expression-only baseline (profile cosine similarity) usable without spot
coordinates.

This module only *uses* the deconvolution API; it never modifies the core
algorithms.
"""
from __future__ import annotations

from typing import Optional

import numpy as np
import pandas as pd

from tissueresolve.results import ReferenceSignature

__all__ = [
    "simulate_pure_profiles",
    "simulate_pairwise_mixtures",
    "simulate_realistic_mixtures",
    "estimate_spillover_matrix",
    "expression_spillover_proxy",
    "compute_spillover_risk",
    "summarize_spillover_partners",
]


# ---------------------------------------------------------------------------
# Simulation
# ---------------------------------------------------------------------------


def simulate_pure_profiles(ref: ReferenceSignature, *, lib_size: float = 1e5) -> pd.DataFrame:
    """Pure pseudobulk expression per cell type as a **genes × types** frame.

    Each column is the reference CPM profile of one type scaled to *lib_size*
    expected counts (non-negative; suitable as input to the bulk solver).
    """
    R = ref.as_R_cpm().astype(np.float64)              # (K, G)
    profiles = (R / 1e6) * float(lib_size)             # (K, G) expected counts
    df = pd.DataFrame(profiles.T, index=[str(g) for g in ref.gene_names],
                      columns=list(ref.cell_types))
    df.index.name = "gene"
    return df


def simulate_pairwise_mixtures(
    ref: ReferenceSignature, *, fractions=(0.5,), lib_size: float = 1e5,
) -> tuple[pd.DataFrame, pd.DataFrame]:
    """Two-type mixtures for every pair, at the given mixing *fractions*.

    Returns ``(counts_genes_x_samples, truth_samples_x_types)``.
    """
    pure = simulate_pure_profiles(ref, lib_size=lib_size)   # genes × types
    cts = list(ref.cell_types)
    cols: dict[str, np.ndarray] = {}
    truth_rows: list[dict[str, float]] = []
    sample_ids: list[str] = []
    for i in range(len(cts)):
        for j in range(i + 1, len(cts)):
            for f in fractions:
                a, b = cts[i], cts[j]
                sid = f"{a}__{b}__{f:.2f}"
                sample_ids.append(sid)
                cols[sid] = f * pure[a].to_numpy() + (1 - f) * pure[b].to_numpy()
                row = {ct: 0.0 for ct in cts}
                row[a], row[b] = f, 1 - f
                truth_rows.append(row)
    counts = pd.DataFrame(cols, index=pure.index)
    truth = pd.DataFrame(truth_rows, index=sample_ids)[cts]
    return counts, truth


def simulate_realistic_mixtures(
    ref: ReferenceSignature, *, n: int = 50, alpha: float = 0.3,
    lib_size: float = 1e5, seed: int = 0,
) -> tuple[pd.DataFrame, pd.DataFrame]:
    """Dirichlet-mixed pseudobulk samples.  Returns ``(counts, truth)``."""
    rng = np.random.default_rng(seed)
    pure = simulate_pure_profiles(ref, lib_size=lib_size)   # genes × types
    cts = list(ref.cell_types)
    P = pure.to_numpy()                                     # (G, K)
    cols: dict[str, np.ndarray] = {}
    truth = np.zeros((n, len(cts)))
    for s in range(n):
        w = rng.dirichlet(np.full(len(cts), alpha))
        truth[s] = w
        cols[f"mix_{s:03d}"] = P @ w
    counts = pd.DataFrame(cols, index=pure.index)
    truth_df = pd.DataFrame(truth, index=list(cols.keys()), columns=cts)
    return counts, truth_df


# ---------------------------------------------------------------------------
# Spillover estimation
# ---------------------------------------------------------------------------


def estimate_spillover_matrix(
    ref: ReferenceSignature, *, lib_size: float = 1e5, config=None,
) -> pd.DataFrame:
    """Deconvolve pure profiles and return the spillover matrix.

    ``spillover[i, j]`` = estimated proportion of type *j* when the input is
    pure type *i*.  Rows are indexed by the true type, columns by the
    predicted type, and each row sums to 1 (it is a proportion vector).  A
    perfectly identifiable reference yields the identity matrix; off-diagonal
    mass is leakage.
    """
    from tissueresolve.api import deconv_bulk

    pure = simulate_pure_profiles(ref, lib_size=lib_size)          # genes × types
    result = deconv_bulk(pure, ref, config=config, n_bootstrap=0)
    est = result.deconv.proportions                                # samples × types
    cts = list(ref.cell_types)
    # Align to (true types × predicted types).
    est = est.reindex(index=cts, columns=cts).fillna(0.0)
    est.index.name = "true_type"
    est.columns.name = "predicted_type"
    return est


def expression_spillover_proxy(ref: ReferenceSignature) -> pd.DataFrame:
    """Expression-only spillover proxy: row-normalised profile cosine similarity.

    A cheap, deconvolution-free baseline for the spatial workflow.  Each row is
    normalised to sum to 1, so it reads like a leakage distribution: the
    diagonal is self-retention and off-diagonal mass approximates confusability.
    """
    R = ref.as_R_log().astype(np.float64)              # (K, G)
    norm = np.linalg.norm(R, axis=1, keepdims=True)
    norm[norm == 0] = 1.0
    cos = (R @ R.T) / (norm @ norm.T)                  # (K, K) in [-1, 1]
    cos = np.clip(cos, 0.0, None)
    row_sums = cos.sum(axis=1, keepdims=True)
    row_sums[row_sums == 0] = 1.0
    M = cos / row_sums
    cts = list(ref.cell_types)
    df = pd.DataFrame(M, index=cts, columns=cts)
    df.index.name = "true_type"
    df.columns.name = "predicted_type"
    return df


def compute_spillover_risk(spillover_matrix: pd.DataFrame) -> pd.DataFrame:
    """Per-type self-retention, spillover risk, and main leaking partner.

    ``spillover_risk = 1 − self_retention`` (off-diagonal mass).  The main
    partner is the type receiving the most leaked mass.
    """
    cts = list(spillover_matrix.index)
    M = spillover_matrix.to_numpy(dtype=float)
    rows = []
    for i, ct in enumerate(cts):
        self_ret = float(M[i, i]) if ct in spillover_matrix.columns else float(M[i, i])
        off = M[i].copy()
        off[i] = -np.inf
        j = int(np.argmax(off))
        partner = list(spillover_matrix.columns)[j]
        partner_frac = float(M[i, j]) if np.isfinite(off[j]) else 0.0
        rows.append({
            "cell_type": ct,
            "self_retention": round(self_ret, 4),
            "spillover_risk": round(1.0 - self_ret, 4),
            "main_partner": partner,
            "main_partner_fraction": round(partner_frac, 4),
        })
    return pd.DataFrame(rows).set_index("cell_type")


def summarize_spillover_partners(spillover_matrix: pd.DataFrame,
                                 *, min_fraction: float = 0.05) -> pd.DataFrame:
    """Long-form table of every off-diagonal leakage above *min_fraction*."""
    cts_r = list(spillover_matrix.index)
    cts_c = list(spillover_matrix.columns)
    rows = []
    for i, a in enumerate(cts_r):
        for j, b in enumerate(cts_c):
            if a == b:
                continue
            frac = float(spillover_matrix.iloc[i, j])
            if frac >= min_fraction:
                rows.append({"true_type": a, "leaks_into": b,
                             "fraction": round(frac, 4)})
    out = pd.DataFrame(rows, columns=["true_type", "leaks_into", "fraction"])
    if not out.empty:
        out = out.sort_values("fraction", ascending=False).reset_index(drop=True)
    return out
