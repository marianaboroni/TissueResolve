"""
Transparent, scenario-separated spatial-benchmark ranking.

Real Visium without spot-level ground truth and synthetic data with known truth
are **never** ranked together: real-data ranking uses marker recovery,
method concordance, spatial structure, stability and runtime/completeness (NO
accuracy); synthetic ranking uses accuracy (RMSE / JSD / Pearson / Spearman /
CCC / dominant accuracy).  Skipped / exported-only / failed tools are excluded
from ranking; imported tools are ranked but flagged.

Every score is a weighted sum of **normalised, direction-oriented** dimensions
with the weights renormalised over the dimensions actually present, so the
calculation is fully transparent (see :func:`explain_spatial_ranking`).
"""
from __future__ import annotations

from typing import Mapping, Optional, Sequence

import numpy as np
import pandas as pd

__all__ = [
    "DEFAULT_REAL_WEIGHTS",
    "DEFAULT_SYNTHETIC_WEIGHTS",
    "jensen_shannon_divergence",
    "marker_recovery_score",
    "normalize_metric_direction",
    "compute_real_spatial_score",
    "compute_synthetic_spatial_accuracy_score",
    "compute_spatial_overall_scorecard",
    "rank_spatial_methods",
    "explain_spatial_ranking",
]

# Real (no ground truth): concordance/structure/usability only — never accuracy.
DEFAULT_REAL_WEIGHTS = {
    "marker_recovery": 0.30, "method_concordance": 0.20,
    "spatial_structure": 0.20, "stability": 0.10,
    "runtime_completeness": 0.10, "interpretability": 0.10,
}
# Synthetic (known truth): accuracy.
DEFAULT_SYNTHETIC_WEIGHTS = {
    "inverse_error": 0.35,          # from RMSE / JSD / MAE (lower error → higher)
    "correlation": 0.30,            # Pearson / Spearman / CCC
    "dominant_accuracy": 0.15,
    "unresolved_aware": 0.10,
    "runtime_completeness": 0.10,
}

_SCORED = ("executed", "success", "imported", "executed_imported")
_NOT_RANKED = ("skipped", "exported_only", "exported_not_run", "failed")
_EPS = 1e-12


# --------------------------------------------------------------------------- #
# metric primitives
# --------------------------------------------------------------------------- #
def jensen_shannon_divergence(p, q, *, base: float = 2.0) -> float:
    """Jensen–Shannon divergence between two non-negative vectors.

    Vectors are L1-normalised to distributions; zeros are handled (terms with
    zero probability contribute zero, no ``log(0)``).  Returns a value in
    ``[0, 1]`` for ``base=2``.  Symmetric; 0 for identical distributions.
    """
    p = np.asarray(p, dtype=np.float64).ravel()
    q = np.asarray(q, dtype=np.float64).ravel()
    if p.shape != q.shape or p.size == 0:
        raise ValueError("p and q must be non-empty, same-length vectors.")
    sp, sq = p.sum(), q.sum()
    if sp <= 0 or sq <= 0:
        return float("nan")
    p = p / sp
    q = q / sq
    m = 0.5 * (p + q)

    def _kl(a, b):
        mask = a > 0
        return float(np.sum(a[mask] * (np.log(a[mask] / b[mask]) / np.log(base))))

    return float(0.5 * _kl(p, m) + 0.5 * _kl(q, m))


def marker_recovery_score(predicted_abundance: pd.DataFrame,
                          marker_expression: pd.DataFrame,
                          marker_sets: Mapping[str, Sequence[str]]) -> pd.Series:
    """Per-cell-type marker-recovery score (a proxy, NOT accuracy validation).

    For each cell type, correlate its predicted per-spot abundance with the mean
    per-spot expression of its marker genes (Pearson over spots).  Real-data
    spot-level ground truth is unavailable, so this measures whether a method's
    spatial abundance tracks known marker expression — a proxy only.

    Parameters
    ----------
    predicted_abundance : spots × cell_types
    marker_expression   : spots × genes  (same spot index)
    marker_sets         : {cell_type: [marker genes]}
    """
    spots = predicted_abundance.index
    expr = marker_expression.reindex(spots)
    out = {}
    for ct in predicted_abundance.columns:
        genes = [g for g in marker_sets.get(str(ct), []) if g in expr.columns]
        if not genes:
            out[ct] = np.nan
            continue
        marker_mean = expr[genes].mean(axis=1)
        pred = predicted_abundance[ct]
        if pred.std() < _EPS or marker_mean.std() < _EPS:
            out[ct] = np.nan
            continue
        out[ct] = float(np.corrcoef(pred.to_numpy(), marker_mean.to_numpy())[0, 1])
    return pd.Series(out, name="marker_recovery")


def normalize_metric_direction(values, *, higher_is_better: bool = True) -> pd.Series:
    """Min–max normalise a metric to ``[0, 1]`` oriented so 1 = best.

    When ``higher_is_better`` is False (e.g. RMSE, JSD, runtime), the metric is
    inverted after normalisation.  A constant column maps to 1.0 (no
    discrimination); NaNs are preserved.
    """
    s = pd.Series(values, dtype="float64")
    valid = s.dropna()
    if valid.empty:
        return s
    lo, hi = valid.min(), valid.max()
    if hi - lo < _EPS:
        norm = pd.Series(1.0, index=s.index)
        norm[s.isna()] = np.nan
        return norm
    norm = (s - lo) / (hi - lo)
    if not higher_is_better:
        norm = 1.0 - norm
    return norm


# --------------------------------------------------------------------------- #
# scoring
# --------------------------------------------------------------------------- #
def _weighted_over_present(dim_scores: pd.DataFrame, weights: Mapping[str, float]) -> pd.Series:
    """Weighted sum over the dimension columns actually present (weights renorm)."""
    present = [d for d in weights if d in dim_scores.columns]
    if not present:
        return pd.Series(np.nan, index=dim_scores.index)
    w = {d: weights[d] for d in present}
    tot = sum(w.values()) or 1.0
    w = {d: v / tot for d, v in w.items()}
    out = pd.Series(0.0, index=dim_scores.index)
    wsum = pd.Series(0.0, index=dim_scores.index)
    for d in present:
        col = dim_scores[d].astype(float)
        mask = col.notna()
        out[mask] += w[d] * col[mask]
        wsum[mask] += w[d]
    # renormalise per-row over the dimensions that were non-NaN for that row
    final = out / wsum.replace(0.0, np.nan)
    return final.round(6)


def _scored_mask(df: pd.DataFrame) -> pd.Series:
    if "status" not in df.columns:
        return pd.Series(True, index=df.index)
    return df["status"].astype(str).isin(_SCORED)


def compute_real_spatial_score(metrics: pd.DataFrame, *,
                               weights: Optional[Mapping[str, float]] = None) -> pd.DataFrame:
    """Real-data (no ground truth) spatial score — concordance/structure only.

    *metrics* is indexed by method with any of the dimension columns
    ``marker_recovery, method_concordance, spatial_structure, stability,
    runtime_completeness, interpretability`` (already higher-is-better, 0–1
    preferred; columns are re-normalised across methods).  **No accuracy
    columns are used even if present.**  Only scored methods get a real_score.
    """
    w = dict(weights or DEFAULT_REAL_WEIGHTS)
    df = metrics.copy()
    dims = pd.DataFrame(index=df.index)
    for d in w:
        if d in df.columns:
            dims[d] = normalize_metric_direction(df[d], higher_is_better=True)
    score = _weighted_over_present(dims, w)
    scored = _scored_mask(df)
    score[~scored] = np.nan
    out = dims.add_suffix("_score")
    out["real_score"] = score
    if "status" in df.columns:
        out["status"] = df["status"]
    return out


def compute_synthetic_spatial_accuracy_score(metrics: pd.DataFrame, *,
                                             weights: Optional[Mapping[str, float]] = None) -> pd.DataFrame:
    """Synthetic (known-truth) spatial ACCURACY score.

    *metrics* indexed by method with any of ``rmse, mae, jsd`` (lower better),
    ``pearson, spearman, ccc`` (higher better), ``dominant_accuracy`` (higher),
    ``unresolved_aware`` (higher), ``runtime_completeness`` (higher).
    """
    w = dict(weights or DEFAULT_SYNTHETIC_WEIGHTS)
    df = metrics.copy()
    dims = pd.DataFrame(index=df.index)
    # inverse error from rmse/mae/jsd (mean of available, inverted)
    err_cols = [c for c in ("rmse", "mae", "jsd") if c in df.columns]
    if err_cols:
        inv = pd.concat([normalize_metric_direction(df[c], higher_is_better=False)
                         for c in err_cols], axis=1).mean(axis=1)
        dims["inverse_error"] = inv
    corr_cols = [c for c in ("pearson", "spearman", "ccc") if c in df.columns]
    if corr_cols:
        dims["correlation"] = pd.concat(
            [normalize_metric_direction(df[c], higher_is_better=True) for c in corr_cols],
            axis=1).mean(axis=1)
    for d, col in (("dominant_accuracy", "dominant_accuracy"),
                   ("unresolved_aware", "unresolved_aware"),
                   ("runtime_completeness", "runtime_completeness")):
        if col in df.columns:
            dims[d] = normalize_metric_direction(df[col], higher_is_better=True)
    score = _weighted_over_present(dims, w)
    scored = _scored_mask(df)
    score[~scored] = np.nan
    out = dims.add_suffix("_score")
    out["synthetic_accuracy_score"] = score
    if "status" in df.columns:
        out["status"] = df["status"]
    return out


def compute_spatial_overall_scorecard(real: Optional[pd.DataFrame] = None,
                                      synthetic: Optional[pd.DataFrame] = None,
                                      *, real_weight: float = 0.5,
                                      synthetic_weight: float = 0.5) -> pd.DataFrame:
    """Combine real + synthetic evidence into a labelled SCORECARD (not truth).

    If only one scenario is present, the scorecard equals that scenario's score.
    The output is explicitly a multi-criteria scorecard, not an objective ranking.
    """
    frames = {}
    if real is not None and "real_score" in real.columns:
        frames["real_score"] = real["real_score"]
    if synthetic is not None and "synthetic_accuracy_score" in synthetic.columns:
        frames["synthetic_accuracy_score"] = synthetic["synthetic_accuracy_score"]
    if not frames:
        return pd.DataFrame(columns=["scorecard"])
    combined = pd.DataFrame(frames)
    if combined.shape[1] == 1:
        combined["scorecard"] = combined.iloc[:, 0]
    else:
        wsum = real_weight + synthetic_weight or 1.0
        rw, sw = real_weight / wsum, synthetic_weight / wsum
        r = combined["real_score"]
        s = combined["synthetic_accuracy_score"]
        # average available scenarios per method (renormalised)
        num = r.fillna(0) * rw + s.fillna(0) * sw
        den = r.notna() * rw + s.notna() * sw
        combined["scorecard"] = (num / den.replace(0.0, np.nan)).round(6)
    combined["scorecard_note"] = "multi-criteria scorecard (NOT objective truth)"
    return combined


def rank_spatial_methods(scores: pd.DataFrame, *, score_col: str,
                         status_col: str = "status") -> pd.DataFrame:
    """Rank methods by *score_col*, excluding skipped/exported/failed tools.

    Imported tools are ranked but flagged ``imported=True``.  When fewer than two
    rankable methods remain, ``rank`` is left NaN and ``ranking_claim`` records
    that no ranking is claimed (metrics are still returned).
    """
    df = scores.copy()
    if status_col in df.columns:
        rankable = ~df[status_col].astype(str).isin(_NOT_RANKED)
    else:
        rankable = pd.Series(True, index=df.index)
    rankable &= df[score_col].notna()
    df["rankable"] = rankable
    df["imported"] = (df[status_col].astype(str).isin(("imported", "executed_imported"))
                      if status_col in df.columns else False)
    n = int(rankable.sum())
    df["rank"] = np.nan
    if n >= 2:
        df.loc[rankable, "rank"] = (
            df.loc[rankable, score_col].rank(ascending=False, method="min"))
        df["ranking_claim"] = f"ranked {n} methods by {score_col}"
    else:
        df["ranking_claim"] = ("only %d rankable method(s); metrics shown, "
                               "no ranking claimed" % n)
    return df.sort_values(score_col, ascending=False)


def explain_spatial_ranking(ranked: pd.DataFrame, *, score_col: str,
                            scenario: str,
                            weights: Optional[Mapping[str, float]] = None) -> str:
    """Human-readable Markdown explanation of a ranking (fully transparent)."""
    lines = [f"## Spatial ranking — {scenario}", ""]
    if weights:
        lines.append("**Weights (renormalised over present dimensions):** "
                     + ", ".join(f"{k}={v}" for k, v in weights.items()))
    no_acc = "real" in scenario.lower() and "synthetic" not in scenario.lower()
    if no_acc:
        lines.append("These are concordance / structure / usability metrics — "
                     "**NOT accuracy** (no spot-level ground truth).")
    lines.append("")
    claim = ranked["ranking_claim"].iloc[0] if "ranking_claim" in ranked.columns and len(ranked) else "—"
    lines.append(f"_{claim}._")
    lines.append("")
    lines.append("| rank | method | score | status |")
    lines.append("|---|---|---|---|")
    for m, r in ranked.iterrows():
        rk = "" if pd.isna(r.get("rank")) else int(r["rank"])
        st = r.get("status", "")
        sc = r.get(score_col)
        sc = "" if pd.isna(sc) else round(float(sc), 4)
        lines.append(f"| {rk} | {m} | {sc} | {st} |")
    return "\n".join(lines)
