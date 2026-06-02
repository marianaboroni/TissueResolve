"""
Batch-effect diagnostics for a single-cell/single-nucleus reference.

Default philosophy: **diagnose, don't blindly correct.**  Aggressive batch
correction can remove real cell-type signal, so these functions quantify
batch/donor/library confounding and marker stability rather than altering the
expression matrix.
"""
from __future__ import annotations

from typing import Optional

import numpy as np
import pandas as pd

__all__ = [
    "BATCH_COL_CANDIDATES", "detect_batch_columns",
    "compute_celltype_batch_confounding", "compute_marker_batch_stability",
    "compute_batch_mixing_score", "batch_aware_reference_summary",
]

BATCH_COL_CANDIDATES = (
    "batch", "donor", "donor_id", "sample", "sample_id", "patient", "patient_id",
    "library", "library_type", "dataset", "center", "technology", "assay",
    "suspension_type",
)


def detect_batch_columns(obs: pd.DataFrame) -> list[str]:
    lower = {c.lower(): c for c in obs.columns}
    return [lower[c] for c in BATCH_COL_CANDIDATES if c in lower]


def compute_celltype_batch_confounding(obs: pd.DataFrame, celltype_col: str,
                                       batch_col: str) -> pd.DataFrame:
    """Per-cell-type batch confounding: dominant-batch fraction & n_batches.

    A cell type whose cells nearly all come from one batch/donor is confounded
    and its signature may reflect batch rather than biology.
    """
    df = pd.DataFrame({"ct": obs[celltype_col].astype(str),
                       "batch": obs[batch_col].astype(str)})
    tab = df.groupby(["ct", "batch"]).size().unstack(fill_value=0)
    frac = tab.div(tab.sum(axis=1).replace(0, np.nan), axis=0).fillna(0.0)
    out = pd.DataFrame({
        "n_batches": (tab > 0).sum(axis=1),
        "dominant_batch_fraction": frac.max(axis=1),
        "n_cells": tab.sum(axis=1),
    })
    out["single_batch_only"] = out["n_batches"] <= 1
    out["confounded"] = out["dominant_batch_fraction"] > 0.9
    return out.sort_values("dominant_batch_fraction", ascending=False)


def compute_marker_batch_stability(ref, mapping: Optional[dict] = None,
                                   *, top_n: int = 50) -> pd.DataFrame:
    """Proxy marker-stability score per cell type using cross-donor CV.

    Uses ``ReferenceSignature.donor_cv`` when present (lower CV → more stable
    markers).  Returns per-cell-type mean CV over each type's top markers; NaN
    when donor CV is unavailable (reported, not silently treated as stable).
    """
    rows = []
    cv = getattr(ref, "donor_cv", None)
    R = ref.as_R_cpm()
    for k, ct in enumerate(ref.cell_types):
        top = np.argsort(R[k])[::-1][:top_n]
        if cv is not None:
            stab = float(np.nanmean(cv[top, k]))
        else:
            stab = float("nan")
        rows.append({"cell_type": ct, "mean_marker_donor_cv": stab,
                     "donor_cv_available": cv is not None})
    return pd.DataFrame(rows).set_index("cell_type")


def compute_batch_mixing_score(obs: pd.DataFrame, batch_col: str) -> float:
    """Shannon-evenness of batch sizes (1 = perfectly even, 0 = one batch)."""
    counts = obs[batch_col].astype(str).value_counts().to_numpy(float)
    p = counts / counts.sum()
    if len(p) <= 1:
        return 0.0
    ent = -(p * np.log(p)).sum()
    return float(ent / np.log(len(p)))


def batch_aware_reference_summary(obs: pd.DataFrame, celltype_col: str,
                                  batch_col: str) -> dict:
    conf = compute_celltype_batch_confounding(obs, celltype_col, batch_col)
    return {
        "batch_col": batch_col,
        "n_batches": int(obs[batch_col].nunique()),
        "batch_mixing_score": round(compute_batch_mixing_score(obs, batch_col), 4),
        "n_confounded_cell_types": int(conf["confounded"].sum()),
        "n_single_batch_cell_types": int(conf["single_batch_only"].sum()),
    }
