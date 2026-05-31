"""
Benchmark metrics — bulk (with ground truth), spatial (no ground truth), and
synthetic spatial (with ground truth).
"""
from __future__ import annotations

from typing import Optional

import numpy as np
import pandas as pd

__all__ = [
    "align", "accuracy_metrics", "per_celltype_rmse", "dominant_accuracy",
    "calibration", "concordance_metrics", "pairwise_method_correlation",
    "proportion_entropy", "near_zero_fraction",
    "family_of_column", "aggregate_to_families", "family_level_metrics",
    "unresolved_columns", "unresolved_aware_metrics", "fine_metrics_resolvable_only",
    "hierarchical_fair_metrics",
]


def align(true_df: pd.DataFrame, est_df: pd.DataFrame):
    """Restrict to common rows and columns."""
    rows = [r for r in true_df.index if r in set(est_df.index)]
    cols = [c for c in true_df.columns if c in set(est_df.columns)]
    if not rows or not cols:
        raise ValueError("no overlap between true and estimated proportions")
    return true_df.loc[rows, cols], est_df.loc[rows, cols]


def _safe_corr(a, b, kind="pearson") -> float:
    from scipy.stats import pearsonr, spearmanr
    a, b = np.asarray(a, float), np.asarray(b, float)
    if a.std() < 1e-12 or b.std() < 1e-12:
        return float("nan")
    return float((pearsonr if kind == "pearson" else spearmanr)(a, b)[0])


def accuracy_metrics(true_df: pd.DataFrame, est_df: pd.DataFrame) -> dict:
    """Overall accuracy on the flattened obs × cell-type matrix."""
    t, e = align(true_df, est_df)
    tv, ev = t.to_numpy(float).ravel(), e.to_numpy(float).ravel()
    diff = ev - tv
    return {
        "pearson": _safe_corr(tv, ev, "pearson"),
        "spearman": _safe_corr(tv, ev, "spearman"),
        "rmse": float(np.sqrt(np.mean(diff ** 2))),
        "mae": float(np.mean(np.abs(diff))),
        "signed_bias": float(np.mean(diff)),
        "n_obs": int(t.shape[0]),
        "n_cell_types": int(t.shape[1]),
    }


def per_celltype_rmse(true_df: pd.DataFrame, est_df: pd.DataFrame) -> pd.DataFrame:
    t, e = align(true_df, est_df)
    rows = []
    for c in t.columns:
        d = e[c].to_numpy(float) - t[c].to_numpy(float)
        rows.append({"cell_type": c,
                     "rmse": float(np.sqrt(np.mean(d ** 2))),
                     "mae": float(np.mean(np.abs(d))),
                     "signed_bias": float(np.mean(d)),
                     "pearson": _safe_corr(t[c], e[c])})
    return pd.DataFrame(rows).set_index("cell_type")


def dominant_accuracy(true_df: pd.DataFrame, est_df: pd.DataFrame, top_k: int = 1) -> float:
    """Fraction of obs whose top-k true types include the top-1 estimated type."""
    t, e = align(true_df, est_df)
    correct = 0
    for i in range(t.shape[0]):
        true_top = set(t.iloc[i].nlargest(top_k).index)
        est_top = e.iloc[i].idxmax()
        correct += int(est_top in true_top)
    return correct / max(t.shape[0], 1)


def calibration(true_df: pd.DataFrame, est_df: pd.DataFrame) -> dict:
    """Slope/intercept of est ~ true (1.0/0.0 is perfect)."""
    t, e = align(true_df, est_df)
    tv, ev = t.to_numpy(float).ravel(), e.to_numpy(float).ravel()
    if tv.std() < 1e-12:
        return {"slope": float("nan"), "intercept": float("nan")}
    slope, intercept = np.polyfit(tv, ev, 1)
    return {"slope": float(slope), "intercept": float(intercept)}


# --- no-ground-truth (real spatial) -----------------------------------------


def concordance_metrics(preds: dict[str, pd.DataFrame]) -> dict:
    """Cross-method concordance summary (no ground truth)."""
    return {"n_methods": len(preds),
            "mean_pairwise_pearson": float(np.nanmean(list(
                pairwise_method_correlation(preds).values()))) if len(preds) > 1 else float("nan")}


def pairwise_method_correlation(preds: dict[str, pd.DataFrame]) -> dict:
    """Pearson r between flattened predictions of each method pair."""
    out = {}
    names = list(preds)
    for i in range(len(names)):
        for j in range(i + 1, len(names)):
            a, b = names[i], names[j]
            try:
                ta, tb = align(preds[a], preds[b])
                out[f"{a}__vs__{b}"] = _safe_corr(
                    ta.to_numpy(float).ravel(), tb.to_numpy(float).ravel())
            except ValueError:
                out[f"{a}__vs__{b}"] = float("nan")
    return out


def proportion_entropy(df: pd.DataFrame) -> pd.Series:
    p = df.to_numpy(float)
    p = np.clip(p, 1e-12, None)
    p = p / p.sum(axis=1, keepdims=True)
    ent = -(p * np.log(p)).sum(axis=1)
    return pd.Series(ent, index=df.index, name="entropy")


def near_zero_fraction(df: pd.DataFrame, eps: float = 1e-3) -> float:
    v = df.to_numpy(float)
    return float((v < eps).mean())


# --- fair hierarchical / family-level / unresolved-aware metrics ------------
#
# Hierarchical methods abstain on non-separable families (mass parked in
# ``unresolved_<family>`` columns).  Judging them only at the fine level
# unfairly scores that honest abstention as a wrong fine prediction.  These
# helpers add (B) family-level metrics for *all* methods, (C) unresolved-aware
# metrics for hierarchical methods, and (D) fine metrics restricted to
# resolvable families.

_UNRESOLVED_PREFIX = "unresolved_"


def family_of_column(col: str, mapping: dict) -> str:
    """Map a prediction column to its broad family.

    ``unresolved_<family>`` columns map to ``<family>``; fine subtypes map via
    *mapping* (unmapped → themselves)."""
    col = str(col)
    if col.startswith(_UNRESOLVED_PREFIX):
        return col[len(_UNRESOLVED_PREFIX):]
    return mapping.get(col, col)


def aggregate_to_families(df: pd.DataFrame, mapping: dict) -> pd.DataFrame:
    """Sum prediction columns (fine + unresolved_*) into broad families.

    Total mass per row is preserved.  Works for flat methods (subtypes → family)
    and hierarchical methods (subtypes + unresolved_<family> → family)."""
    fam_of = {c: family_of_column(c, mapping) for c in df.columns}
    families = []
    for c in df.columns:
        f = fam_of[c]
        if f not in families:
            families.append(f)
    out = pd.DataFrame(index=df.index)
    for f in families:
        members = [c for c in df.columns if fam_of[c] == f]
        out[f] = df[members].sum(axis=1)
    return out


def family_level_metrics(true_fine: pd.DataFrame, est: pd.DataFrame,
                         mapping: dict) -> dict:
    """Accuracy after aggregating both truth and estimate to broad families."""
    true_fam = aggregate_to_families(true_fine, mapping)
    est_fam = aggregate_to_families(est, mapping)
    return accuracy_metrics(true_fam, est_fam)


def unresolved_columns(df: pd.DataFrame) -> list[str]:
    return [c for c in df.columns if str(c).startswith(_UNRESOLVED_PREFIX)]


def unresolved_aware_metrics(est: pd.DataFrame, mapping: dict,
                             high_risk_families: set) -> dict:
    """Quantify abstention behavior for a hierarchical prediction.

    *high_risk_families* are families known (from the reference) to be
    non-separable — abstaining on those is *appropriate*, not an error.

    Returns unresolved-mass fraction and the precision/recall of the method's
    unresolved families against the high-risk set.  Honest abstention on
    high-risk families is rewarded, not penalised.
    """
    ucols = unresolved_columns(est)
    unresolved_fams = {c[len(_UNRESOLVED_PREFIX):] for c in ucols}
    total = est.to_numpy(float).sum()
    umass = est[ucols].to_numpy(float).sum() if ucols else 0.0
    hr = set(high_risk_families)
    tp = len(unresolved_fams & hr)
    precision = tp / len(unresolved_fams) if unresolved_fams else float("nan")
    recall = tp / len(hr) if hr else float("nan")
    return {
        "unresolved_mass_fraction": float(umass / total) if total else 0.0,
        "n_unresolved_families": len(unresolved_fams),
        "n_high_risk_families": len(hr),
        "unresolved_precision": precision,   # of abstentions, fraction truly high-risk
        "unresolved_recall": recall,         # of high-risk families, fraction abstained
        "abstains": bool(ucols),
    }


def fine_metrics_resolvable_only(true_fine: pd.DataFrame, est: pd.DataFrame,
                                 mapping: dict, unresolved_families: set) -> dict:
    """Fine-level accuracy restricted to families the method did NOT abstain on.

    This is the fair fine-level view for a hierarchical method: it is judged on
    the subtypes it actually claims, not on the ones it honestly abstained on.
    """
    fam_of = {c: family_of_column(c, mapping) for c in est.columns}
    keep = [c for c in est.columns
            if not str(c).startswith(_UNRESOLVED_PREFIX)
            and fam_of[c] not in set(unresolved_families)]
    keep = [c for c in keep if c in true_fine.columns]
    if not keep:
        return {"pearson": float("nan"), "rmse": float("nan"),
                "n_cell_types": 0, "note": "no resolvable fine subtypes"}
    return accuracy_metrics(true_fine[keep], est[keep])


def hierarchical_fair_metrics(true_fine: pd.DataFrame, est: pd.DataFrame,
                              mapping: dict, high_risk_families: set) -> dict:
    """Bundle fine (resolvable-only), family-level, and unresolved-aware metrics."""
    unresolved_fams = {c[len(_UNRESOLVED_PREFIX):] for c in unresolved_columns(est)}
    fam = family_level_metrics(true_fine, est, mapping)
    fine_res = fine_metrics_resolvable_only(true_fine, est, mapping, unresolved_fams)
    unaware = unresolved_aware_metrics(est, mapping, high_risk_families)
    return {
        "family_pearson": fam["pearson"], "family_rmse": fam["rmse"],
        "fine_resolvable_pearson": fine_res["pearson"],
        "fine_resolvable_rmse": fine_res["rmse"],
        **unaware,
    }
