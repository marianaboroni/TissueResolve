"""
Composite benchmark score across accuracy, robustness, usability,
interpretability, resolution-awareness, and runtime/resource.

Only **executed** or **imported** tools receive a final score; skipped /
exported-only / failed tools do not.  For real spatial data (no ground truth),
the accuracy dimension is excluded from the final score (weights renormalised).
The calculation is transparent — every dimension is written out.
"""
from __future__ import annotations

from pathlib import Path
from typing import Optional

import numpy as np
import pandas as pd

DEFAULT_WEIGHTS = {
    "accuracy": 0.35, "robustness": 0.20, "usability": 0.15,
    "interpretability": 0.15, "resolution_awareness": 0.10, "runtime_resource": 0.05,
}

# qualitative per-method priors for non-accuracy dimensions (0..1)
_PRIORS = {
    "TissueResolve_auto": dict(usability=0.9, interpretability=0.95,
                               resolution_awareness=0.8, robustness=0.85),
    "TissueResolve_hierarchical": dict(usability=0.9, interpretability=1.0,
                                       resolution_awareness=1.0, robustness=0.85),
    "TissueResolve_flat": dict(usability=0.9, interpretability=0.8,
                               resolution_awareness=0.4, robustness=0.7),
    "NNLS_baseline": dict(usability=0.7, interpretability=0.4,
                          resolution_awareness=0.1, robustness=0.5),
}
_DEFAULT_PRIOR = dict(usability=0.6, interpretability=0.5,
                      resolution_awareness=0.3, robustness=0.5)


def load_weights(path: Optional[Path] = None) -> dict:
    if path is None:
        path = Path(__file__).resolve().parents[1] / "configs" / "composite_score_weights.yaml"
    try:
        import yaml
        w = yaml.safe_load(Path(path).read_text())
        return {k: float(v) for k, v in w.items()}
    except Exception:
        return dict(DEFAULT_WEIGHTS)


def _runtime_score(runtimes: pd.Series) -> pd.Series:
    r = runtimes.astype(float).clip(lower=1e-6)
    return (r.min() / r).clip(0, 1)   # fastest = 1.0


def compute_composite_scores(per_method: pd.DataFrame, *, modality: str = "bulk",
                             has_ground_truth: bool = True,
                             weights: Optional[dict] = None) -> pd.DataFrame:
    """*per_method* indexed by method with columns possibly including
    ``status``, ``modality`` (``"bulk"`` / ``"spatial"``), ``accuracy`` (0..1),
    ``runtime_seconds``, ``robustness``.

    Returns a scored table.  Scoring and ranking are **per-modality**: bulk and
    spatial methods are never pooled into a single leaderboard, because their
    final scores are not comparable (bulk has pseudobulk ground truth and a real
    accuracy dimension; real spatial has none, so its accuracy is excluded and
    the remaining weights are renormalised).  ``rank`` is therefore computed
    *within* each modality, and the runtime dimension is normalised within each
    modality (the fastest bulk tool and the fastest spatial tool each score
    1.0).  Only executed/imported tools are scored.

    ``modality`` / ``has_ground_truth`` are fallbacks used only when the frame
    lacks a ``modality`` column (and per-row ``accuracy``).
    """
    base_w = dict(weights or DEFAULT_WEIGHTS)
    df = per_method.copy()
    scored_mask = df.get("status", pd.Series("executed", index=df.index)).isin(
        ["executed", "success", "imported", "executed_imported"])

    # Per-row modality (column wins; else fall back to the scalar argument).
    if "modality" in df.columns:
        row_modality = df["modality"].astype(str)
    else:
        row_modality = pd.Series(modality, index=df.index)

    # Runtime is normalised *within* each modality so the two groups are scored
    # on their own scale rather than the global fastest tool.
    runtime_raw = df.get("runtime_seconds", pd.Series(1.0, index=df.index))
    runtime = pd.Series(np.nan, index=df.index, dtype=float)
    for mod, idx in row_modality.groupby(row_modality).groups.items():
        runtime.loc[idx] = _runtime_score(runtime_raw.loc[idx])

    rows = []
    for m in df.index:
        prior = _PRIORS.get(m, _DEFAULT_PRIOR)
        acc = float(df.loc[m].get("accuracy", np.nan)) if "accuracy" in df.columns else np.nan
        # ``has_ground_truth`` is a global override (False => nobody gets an
        # accuracy dimension).  When ground truth is allowed, a row earns the
        # accuracy dimension only if it actually carries an accuracy value
        # (pseudobulk does; real spatial does not) — which is what de-leaks the
        # spatial methods out of the accuracy-bearing bulk ranking.
        if "accuracy" in df.columns:
            row_has_gt = has_ground_truth and (not np.isnan(acc))
        else:
            row_has_gt = has_ground_truth
        w = dict(base_w)
        if not row_has_gt:
            w.pop("accuracy", None)
        total_w = sum(w.values())
        w = {k: v / total_w for k, v in w.items()}
        dims = {
            "accuracy_score": acc if row_has_gt else np.nan,
            "robustness_score": float(df.loc[m].get("robustness", prior["robustness"])),
            "usability_score": prior["usability"],
            "interpretability_score": prior["interpretability"],
            "resolution_awareness_score": prior["resolution_awareness"],
            "runtime_resource_score": float(runtime.loc[m]),
        }
        if scored_mask.loc[m]:
            fs = 0.0
            for dim, wt in w.items():
                v = dims.get(f"{dim}_score", np.nan)
                if not np.isnan(v):
                    fs += wt * v
            final = round(fs, 4)
        else:
            final = np.nan
        rows.append({"method": m, "modality": row_modality.loc[m],
                     "status": df.loc[m].get("status", "executed"),
                     **{k: round(v, 4) if not (isinstance(v, float) and np.isnan(v)) else np.nan
                        for k, v in dims.items()},
                     "final_score": final})
    out = pd.DataFrame(rows).set_index("method")
    # Rank within modality — bulk and spatial are separate leaderboards.
    out["rank"] = out.groupby("modality")["final_score"].rank(
        ascending=False, method="min")
    return out.sort_values(["modality", "final_score"], ascending=[True, False])
