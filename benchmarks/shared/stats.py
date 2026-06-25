"""Statistical comparison helpers for benchmarking (PART 14).

Bootstrap confidence intervals, paired Wilcoxon signed-rank tests, paired effect
size, and bootstrap rank stability.  All deterministic given a seed; no global
RNG state is touched.
"""
from __future__ import annotations

from typing import Callable

import numpy as np
import pandas as pd

__all__ = ["bootstrap_ci", "paired_wilcoxon", "matched_pairs_effect_size",
           "rank_stability"]


def bootstrap_ci(values, *, statistic: Callable = np.mean, n_boot: int = 2000,
                 alpha: float = 0.05, seed: int = 0) -> dict:
    """Percentile bootstrap CI for *statistic* over *values* (resampled with replacement)."""
    v = np.asarray([x for x in values if x is not None and np.isfinite(x)], float)
    if v.size == 0:
        return {"point": float("nan"), "lo": float("nan"), "hi": float("nan"), "n": 0}
    rng = np.random.default_rng(seed)
    idx = rng.integers(0, v.size, size=(n_boot, v.size))
    boot = np.array([statistic(v[i]) for i in idx])
    return {"point": float(statistic(v)),
            "lo": float(np.percentile(boot, 100 * alpha / 2)),
            "hi": float(np.percentile(boot, 100 * (1 - alpha / 2))),
            "n": int(v.size)}


def paired_wilcoxon(a, b) -> dict:
    """Paired Wilcoxon signed-rank test between two equal-length metric vectors.

    Pairs with NaN in either arm are dropped.  Returns statistic, two-sided p,
    the matched-pairs rank-biserial effect size, and n.
    """
    a = np.asarray(a, float)
    b = np.asarray(b, float)
    mask = np.isfinite(a) & np.isfinite(b)
    a, b = a[mask], b[mask]
    if a.size < 1 or np.allclose(a, b):
        return {"statistic": float("nan"), "p_value": float("nan"),
                "effect_size": 0.0, "n": int(a.size)}
    from scipy.stats import wilcoxon
    try:
        stat, p = wilcoxon(a, b, zero_method="wilcox", correction=False,
                           mode="auto")
    except ValueError:
        return {"statistic": float("nan"), "p_value": float("nan"),
                "effect_size": 0.0, "n": int(a.size)}
    return {"statistic": float(stat), "p_value": float(p),
            "effect_size": matched_pairs_effect_size(a, b), "n": int(a.size)}


def matched_pairs_effect_size(a, b) -> float:
    """Rank-biserial correlation for matched pairs (signed; range [-1, 1]).

    Positive ⇒ *a* tends to exceed *b*.  Based on the signed-rank decomposition.
    """
    a = np.asarray(a, float)
    b = np.asarray(b, float)
    d = a - b
    d = d[d != 0]
    if d.size == 0:
        return 0.0
    ranks = pd.Series(np.abs(d)).rank().to_numpy()
    r_plus = ranks[d > 0].sum()
    r_minus = ranks[d < 0].sum()
    total = r_plus + r_minus
    return float((r_plus - r_minus) / total) if total > 0 else 0.0


def rank_stability(per_replicate: pd.DataFrame, *, higher_is_better: bool = True,
                   n_boot: int = 2000, seed: int = 0) -> pd.DataFrame:
    """Bootstrap rank distribution of methods across replicates.

    *per_replicate*: rows = replicates (e.g. split×scenario), columns = methods,
    values = a single metric.  Resamples replicates with replacement, ranks the
    per-resample column means, and reports each method's mean rank and the
    fraction of resamples it is ranked #1.
    """
    methods = list(per_replicate.columns)
    X = per_replicate.to_numpy(float)
    n = X.shape[0]
    rng = np.random.default_rng(seed)
    ranks = np.zeros((n_boot, len(methods)))
    top1 = np.zeros(len(methods))
    for b in range(n_boot):
        idx = rng.integers(0, n, size=n)
        means = np.nanmean(X[idx], axis=0)
        order = (-means if higher_is_better else means)
        r = pd.Series(order).rank(method="average").to_numpy()
        ranks[b] = r
        best = np.nanargmax(means) if higher_is_better else np.nanargmin(means)
        top1[best] += 1
    return pd.DataFrame({
        "method": methods,
        "mean_rank": ranks.mean(axis=0),
        "rank_sd": ranks.std(axis=0),
        "prob_best": top1 / n_boot,
    }).set_index("method").sort_values("mean_rank")
