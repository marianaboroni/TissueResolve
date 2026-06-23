"""Spatial-fidelity metrics for spatial deconvolution benchmarks (PART 12/13).

Beyond per-spot accuracy, spatial deconvolution must preserve *spatial
structure*: autocorrelation (Moran's I / Geary's C), tissue domains, and local
(neighbourhood) accuracy — and must not over-smooth.  These operate on a
proportions matrix (spots × cell types) plus spot coordinates, and are robust to
the zeros pervasive in deconvolution output.  All deterministic.
"""
from __future__ import annotations

import numpy as np
import pandas as pd

__all__ = ["knn_weights", "morans_i", "gearys_c", "morans_i_per_column",
           "morans_i_preservation", "oversmoothing_score", "local_rmse",
           "domain_recovery_ari", "spatial_fidelity_metrics",
           "boundary_metrics", "rare_niche_sensitivity"]


def knn_weights(coords: np.ndarray, k: int = 6) -> tuple:
    """Row-standardised k-NN spatial weights as (indices, weights) per spot.

    Returns ``(nbr_idx, w)`` where ``nbr_idx[i]`` are the k nearest neighbours of
    spot *i* (excluding itself) and ``w`` is ``1/k`` (row-standardised).
    """
    from scipy.spatial import cKDTree
    coords = np.asarray(coords, float)
    n = coords.shape[0]
    k = min(k, n - 1)
    tree = cKDTree(coords)
    _, idx = tree.query(coords, k=k + 1)      # includes self at col 0
    nbr = idx[:, 1:]
    return nbr, 1.0 / k


def morans_i(x: np.ndarray, nbr: np.ndarray, w: float) -> float:
    """Moran's I spatial autocorrelation for one variable with k-NN weights."""
    x = np.asarray(x, float)
    n = x.size
    xb = x.mean()
    d = x - xb
    denom = np.sum(d * d)
    if denom <= 0:
        return float("nan")
    # row-standardised weights: W = n (each row sums to 1), so n/W = 1
    num = np.sum(d[:, None] * (w * d[nbr]))   # sum_i d_i * sum_j w_ij d_j
    return float((n / (n * 1.0)) * num / denom)  # n/W * num/denom, W=n


def gearys_c(x: np.ndarray, nbr: np.ndarray, w: float) -> float:
    """Geary's C (≈0 strong positive autocorr, 1 none, >1 negative)."""
    x = np.asarray(x, float)
    n = x.size
    d = x - x.mean()
    denom = 2.0 * np.sum(d * d)               # 2*W*var with W=n -> 2*sum d^2 (W=n cancels)
    if denom <= 0:
        return float("nan")
    diff2 = (x[:, None] - x[nbr]) ** 2        # (x_i - x_j)^2 over k-NN
    # row-standardised weights ⇒ total weight W = n, so denom = 2 n sum d^2
    return float((n - 1) * np.sum(w * diff2) / (2.0 * n * np.sum(d * d)))


def morans_i_per_column(df: pd.DataFrame, coords, k: int = 6) -> pd.Series:
    nbr, w = knn_weights(coords, k)
    return pd.Series({c: morans_i(df[c].to_numpy(float), nbr, w) for c in df.columns},
                     name="morans_i")


def morans_i_preservation(true_df: pd.DataFrame, pred_df: pd.DataFrame,
                          coords, k: int = 6) -> dict:
    """How well predicted spatial autocorrelation matches the truth, per cell type."""
    cols = [c for c in true_df.columns if c in pred_df.columns]
    it = morans_i_per_column(true_df[cols], coords, k)
    ip = morans_i_per_column(pred_df[cols], coords, k)
    diff = (ip - it).abs()
    valid = it.notna() & ip.notna()
    corr = float(np.corrcoef(it[valid], ip[valid])[0, 1]) if valid.sum() > 2 else float("nan")
    return {"morans_i_mae": float(diff[valid].mean()),
            "morans_i_corr_true_pred": corr,
            "mean_true_morans_i": float(it[valid].mean()),
            "mean_pred_morans_i": float(ip[valid].mean())}


def oversmoothing_score(true_df, pred_df, coords, k: int = 6) -> float:
    """Ratio of mean predicted to mean true Moran's I (>1 over-smoothed, <1 under).

    Uses the **aggregate** mean autocorrelation rather than a per-cell-type ratio:
    per-type ratios blow up when a type's true Moran's I is near zero, so the
    aggregate ratio is the stable, interpretable form.
    """
    cols = [c for c in true_df.columns if c in pred_df.columns]
    it = morans_i_per_column(true_df[cols], coords, k)
    ip = morans_i_per_column(pred_df[cols], coords, k)
    valid = it.notna() & ip.notna()
    mt = float(it[valid].mean())
    if abs(mt) < 1e-9:
        return float("nan")
    return float(ip[valid].mean() / mt)


def local_rmse(true_df: pd.DataFrame, pred_df: pd.DataFrame, coords, k: int = 6) -> float:
    """Mean over spots of the RMSE within each spot's k-NN neighbourhood.

    Captures *local* accuracy / boundary fidelity rather than global RMSE.
    """
    cols = [c for c in true_df.columns if c in pred_df.columns]
    t = true_df[cols].to_numpy(float)
    p = pred_df[cols].to_numpy(float)
    nbr, _ = knn_weights(coords, k)
    errs = []
    for i in range(t.shape[0]):
        block = np.r_[i, nbr[i]]
        errs.append(np.sqrt(np.mean((t[block] - p[block]) ** 2)))
    return float(np.mean(errs))


def domain_recovery_ari(pred_df: pd.DataFrame, domain_labels, *,
                        n_domains: int = None, seed: int = 0) -> float:
    """Adjusted Rand Index between KMeans clusters of predicted composition and
    the true spatial domain labels (how well composition recovers tissue domains).
    """
    from sklearn.cluster import KMeans
    from sklearn.metrics import adjusted_rand_score
    labels = np.asarray(domain_labels)
    k = n_domains or len(np.unique(labels))
    km = KMeans(n_clusters=k, n_init=10, random_state=seed).fit(pred_df.to_numpy(float))
    return float(adjusted_rand_score(labels, km.labels_))


def boundary_metrics(true_df: pd.DataFrame, pred_df: pd.DataFrame, coords,
                     domain_labels, *, mapping=None, k: int = 6) -> dict:
    """Boundary preservation + edge blurring (needs domain labels).

    A spot is a TRUE boundary if any k-NN neighbour has a different domain label.
    A spot is a PRED boundary if its predicted dominant (broad) cell type differs
    from any neighbour's.  `boundary_f1` = F1 of predicted vs true boundary spots.
    `edge_blurring` = boundary local-RMSE minus interior local-RMSE (↑ = more
    blurring at edges; lower/negative = sharper).
    """
    cols = [c for c in true_df.columns if c in pred_df.columns
            and not str(c).startswith("unresolved_")]
    t = true_df[cols]; p = pred_df.reindex(index=t.index, columns=cols).fillna(0.0)
    # aggregate to broad for the dominant-domain comparison if a mapping is given
    if mapping is not None:
        fam = {c: str(mapping.get(c, c)) for c in cols}
        fams = sorted(set(fam.values()))
        tb = pd.DataFrame({f: t[[c for c in cols if fam[c] == f]].sum(1) for f in fams})
        pb = pd.DataFrame({f: p[[c for c in cols if fam[c] == f]].sum(1) for f in fams})
    else:
        tb, pb = t, p
    nbr, _ = knn_weights(coords, k)
    dom = np.asarray(domain_labels)
    pred_dom = pb.to_numpy().argmax(1)
    n = len(dom)
    true_b = np.array([any(dom[j] != dom[i] for j in nbr[i]) for i in range(n)])
    pred_b = np.array([any(pred_dom[j] != pred_dom[i] for j in nbr[i]) for i in range(n)])
    from sklearn.metrics import f1_score
    f1 = float(f1_score(true_b, pred_b)) if true_b.any() and pred_b.any() else float("nan")
    tv, pv = t.to_numpy(float), p.to_numpy(float)
    lr = []
    for i in range(n):
        blk = np.r_[i, nbr[i]]
        lr.append(np.sqrt(np.mean((tv[blk] - pv[blk]) ** 2)))
    lr = np.array(lr)
    edge_blur = float(lr[true_b].mean() - lr[~true_b].mean()) if true_b.any() and (~true_b).any() else float("nan")
    return {"boundary_f1": f1, "edge_blurring": edge_blur,
            "boundary_local_rmse": float(lr[true_b].mean()) if true_b.any() else float("nan"),
            "interior_local_rmse": float(lr[~true_b].mean()) if (~true_b).any() else float("nan")}


def rare_niche_sensitivity(true_df: pd.DataFrame, pred_df: pd.DataFrame, domain_labels,
                           rare_type, *, niche_label="niche", thresh=0.01) -> float:
    """Recall of the rare cell type in the niche spots (fraction of niche spots where
    the rare type is predicted above *thresh*)."""
    dom = np.asarray(domain_labels)
    niche = dom == niche_label
    if not niche.any() or rare_type not in pred_df.columns:
        return float("nan")
    return float((pred_df[rare_type].to_numpy()[niche] > thresh).mean())


def spatial_fidelity_metrics(true_df: pd.DataFrame, pred_df: pd.DataFrame,
                             coords, domain_labels=None, k: int = 6) -> dict:
    """Bundle: Moran's I preservation, oversmoothing, local RMSE, domain ARI."""
    out = {}
    out.update(morans_i_preservation(true_df, pred_df, coords, k))
    out["oversmoothing_score"] = oversmoothing_score(true_df, pred_df, coords, k)
    out["local_rmse"] = local_rmse(true_df, pred_df, coords, k)
    if domain_labels is not None:
        out["domain_recovery_ari"] = domain_recovery_ari(pred_df, domain_labels)
    return out
