"""
Normalization detection and compatibility for benchmark inputs.

Detects whether a matrix is raw counts, CPM/TPM-like, log1p-normalized, or
scaled, and converts when a method requires a specific scale.  All decisions
are recorded (never silent).
"""
from __future__ import annotations

from typing import Any

import numpy as np
import pandas as pd

__all__ = [
    "detect_matrix_scale", "detect_normalization_status",
    "validate_method_normalization_requirements", "convert_counts_to_cpm",
    "convert_counts_to_log_cpm", "preserve_raw_counts_if_available",
    "normalization_report",
]


def _as_array(x):
    try:
        import scipy.sparse as sp
        if sp.issparse(x):
            return x
    except Exception:
        pass
    return np.asarray(x)


def detect_matrix_scale(x) -> dict:
    """Summarise a matrix's scale: integer-ness, range, zero fraction, max."""
    arr = _as_array(x)
    try:
        import scipy.sparse as sp
        if sp.issparse(arr):
            data = arr.data
            vmax = float(data.max()) if data.size else 0.0
            vmin = float(data.min()) if data.size else 0.0
            nnz = arr.nnz
            total = arr.shape[0] * arr.shape[1]
            zero_frac = 1.0 - (nnz / max(total, 1))
            integerish = bool(np.allclose(data, np.round(data))) if data.size else True
            row_sums = np.asarray(arr.sum(axis=1)).ravel()
        else:
            vmax, vmin = float(arr.max()), float(arr.min())
            zero_frac = float((arr == 0).mean())
            integerish = bool(np.allclose(arr, np.round(arr)))
            row_sums = arr.sum(axis=1)
    except Exception:
        arr = np.asarray(arr, float)
        vmax, vmin = float(arr.max()), float(arr.min())
        zero_frac = float((arr == 0).mean())
        integerish = bool(np.allclose(arr, np.round(arr)))
        row_sums = arr.sum(axis=1)
    return {
        "max": vmax, "min": vmin, "zero_fraction": round(zero_frac, 4),
        "integerish": integerish,
        "row_sum_median": float(np.median(row_sums)) if len(row_sums) else 0.0,
    }


def detect_normalization_status(x) -> str:
    """Classify a matrix scale into a normalization status label."""
    s = detect_matrix_scale(x)
    if s["min"] < 0:
        return "scaled"  # z-scored / centered
    if s["integerish"] and s["max"] > 30:
        return "counts"
    rs = s["row_sum_median"]
    if abs(rs - 1e6) / 1e6 < 0.05:
        return "cpm"
    if s["max"] <= 30 and not s["integerish"]:
        return "log_normalized"
    if s["integerish"]:
        return "counts"
    return "unknown"


def validate_method_normalization_requirements(method, status: str) -> list[str]:
    warns = []
    if getattr(method, "requires_raw_counts", False) and status not in ("counts", "unknown"):
        warns.append(f"{method.name} expects raw counts but input is '{status}'.")
    if getattr(method, "requires_normalized_input", False) and status == "counts":
        warns.append(f"{method.name} expects normalized input but got raw counts.")
    return warns


def convert_counts_to_cpm(df: pd.DataFrame) -> pd.DataFrame:
    """genes × samples → CPM (column-wise)."""
    col_sums = df.sum(axis=0).replace(0, np.nan)
    return (df.div(col_sums, axis=1) * 1e6).fillna(0.0)


def convert_counts_to_log_cpm(df: pd.DataFrame) -> pd.DataFrame:
    return np.log1p(convert_counts_to_cpm(df))


def preserve_raw_counts_if_available(adata) -> tuple[Any, str]:
    """Return (counts_matrix, source) preferring layers['counts']/raw over X."""
    if getattr(adata, "layers", None) is not None and "counts" in adata.layers:
        return adata.layers["counts"], "layers['counts']"
    if getattr(adata, "raw", None) is not None:
        return adata.raw.X, "raw.X"
    return adata.X, "X"


def normalization_report(inputs: dict[str, Any]) -> pd.DataFrame:
    """Build a tidy normalization report from {name: matrix}."""
    rows = []
    for name, mat in inputs.items():
        s = detect_matrix_scale(mat)
        rows.append({"input": name, "status": detect_normalization_status(mat),
                     **s})
    return pd.DataFrame(rows).set_index("input")
