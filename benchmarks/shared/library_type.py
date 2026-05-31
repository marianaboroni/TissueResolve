"""
Single-cell reference library-type detection and synthetic library stress tests.

Distinguishes scRNA-seq (3'/5'), single-nucleus RNA-seq, SMART-seq, and mixed
references from metadata when available, and can build controlled synthetic
scRNA/snRNA/mixed references for stress testing (clearly labelled as synthetic,
not biological validation).
"""
from __future__ import annotations

from typing import Optional

import numpy as np
import pandas as pd

__all__ = [
    "LIBRARY_TYPE_COLS", "detect_library_type", "summarize_library_composition",
    "create_reference_subset_by_library_type", "simulate_library_shift",
    "cell_type_by_library_type",
]

LIBRARY_TYPE_COLS = (
    "library_type", "assay", "technology", "suspension_type", "protocol",
    "cell_or_nucleus", "method",
)

_SN_KEYWORDS = ("nucleus", "snrna", "sn-rna", "single nucleus", "single-nucleus", "nuclei")
_SMARTSEQ = ("smart-seq", "smartseq", "smart seq", "ss2")
_TENX3 = ("3'", "3 prime", "three prime", "10x 3", "v2", "v3")
_TENX5 = ("5'", "5 prime", "five prime", "10x 5")


def _classify(value: str) -> str:
    v = str(value).lower()
    if any(k in v for k in _SN_KEYWORDS):
        return "single_nucleus"
    if any(k in v for k in _SMARTSEQ):
        return "smartseq"
    if any(k in v for k in _TENX5):
        return "scrna_10x_5p"
    if any(k in v for k in _TENX3) or "10x" in v or "chromium" in v:
        return "scrna_10x_3p"
    if "cell" in v:
        return "scrna"
    return "unknown"


def detect_library_type(obs: pd.DataFrame, col: Optional[str] = None) -> dict:
    """Detect per-cell library type and an overall label.

    Returns ``{"column": col|None, "per_value": {...}, "overall": label,
    "is_mixed": bool, "confidence": float}``.
    """
    use = None
    if col and col != "auto" and col in obs.columns:
        use = col
    else:
        lower = {c.lower(): c for c in obs.columns}
        for cand in LIBRARY_TYPE_COLS:
            if cand in lower:
                use = lower[cand]
                break
    if use is None:
        return {"column": None, "per_value": {}, "overall": "unknown",
                "is_mixed": False, "confidence": 0.0}
    vc = obs[use].astype(str).value_counts()
    classes = {v: _classify(v) for v in vc.index}
    class_counts: dict[str, int] = {}
    for v, n in vc.items():
        class_counts[classes[v]] = class_counts.get(classes[v], 0) + int(n)
    known = {k: v for k, v in class_counts.items() if k != "unknown"}
    distinct = set(known)
    is_mixed = len({c for c in distinct if c.startswith("scrna")} | (
        {"single_nucleus"} if "single_nucleus" in distinct else set())) > 1 or (
        "single_nucleus" in distinct and any(c.startswith("scrna") for c in distinct))
    overall = "mixed" if is_mixed else (max(class_counts, key=class_counts.get)
                                        if class_counts else "unknown")
    total = sum(class_counts.values())
    conf = (max(class_counts.values()) / total) if total else 0.0
    return {"column": use, "per_value": class_counts, "overall": overall,
            "is_mixed": is_mixed, "confidence": round(float(conf), 3)}


def summarize_library_composition(obs: pd.DataFrame, library_col: str,
                                  celltype_col: str) -> pd.DataFrame:
    """cell-type × library-type cell-count table."""
    df = pd.DataFrame({"cell_type": obs[celltype_col].astype(str),
                       "library": obs[library_col].astype(str)})
    return (df.groupby(["cell_type", "library"]).size()
            .unstack(fill_value=0))


def cell_type_by_library_type(obs: pd.DataFrame, library_col: str,
                              celltype_col: str) -> pd.DataFrame:
    """Per-cell-type fraction coming from each library type (confounding view)."""
    tab = summarize_library_composition(obs, library_col, celltype_col)
    frac = tab.div(tab.sum(axis=1).replace(0, np.nan), axis=0).fillna(0.0)
    frac["dominant_library_fraction"] = frac.max(axis=1)
    return frac


def create_reference_subset_by_library_type(adata, library_col: str, value: str):
    """Subset an AnnData to cells of a given library type."""
    mask = adata.obs[library_col].astype(str) == value
    return adata[mask].copy()


def simulate_library_shift(R_cpm: np.ndarray, gene_names: list[str], *,
                           kind: str = "single_nucleus",
                           seed: int = 0) -> np.ndarray:
    """Apply a SYNTHETIC, label-only library-type shift to a CPM matrix.

    This is a *stress test*, not biological validation: snRNA is modelled as
    relative depletion of a random set of genes (proxy for cytoplasmic/3' bias)
    and slight dispersion, then re-normalised to CPM.  Returns a new matrix.
    """
    rng = np.random.default_rng(seed)
    R = np.asarray(R_cpm, float).copy()
    G = R.shape[1]
    if kind == "single_nucleus":
        depleted = rng.choice(G, size=max(1, G // 5), replace=False)
        factor = np.ones(G)
        factor[depleted] = rng.uniform(0.2, 0.6, size=len(depleted))
        R = R * factor[None, :]
    R = R * rng.uniform(0.9, 1.1, size=R.shape)
    R = np.clip(R, 0, None)
    row = R.sum(axis=1, keepdims=True)
    row[row == 0] = 1.0
    return (R / row) * 1e6
