"""TCGA prediction sparsity audit

Scans prediction TSVs (samples × cell types) under
`benchmarks/outputs/bulk/predictions/` and computes per-sample and
method-level sparsity/complexity metrics described in the validation plan.

Writes results to `benchmarks/outputs/tcga_prediction_complexity_audit_by_sample.tsv`
and `benchmarks/outputs/tcga_prediction_complexity_audit_by_method.tsv`.
"""
from __future__ import annotations

import math
from pathlib import Path
from typing import Dict, List

import numpy as np
import pandas as pd


ROOT = Path(__file__).resolve().parents[1]
PRED_DIR = ROOT / "outputs" / "bulk" / "predictions"
OUT_DIR = ROOT / "outputs"
OUT_DIR.mkdir(parents=True, exist_ok=True)


def gini(x: np.ndarray) -> float:
    # Gini coefficient for a distribution vector (non-negative)
    x = np.asarray(x, dtype=float)
    if x.sum() == 0:
        return 0.0
    x = x[x >= 0]
    if x.size == 0:
        return 0.0
    sorted_x = np.sort(x)
    n = x.size
    cumx = np.cumsum(sorted_x)
    g = (2.0 * np.sum((np.arange(1, n + 1) * sorted_x))) / (n * cumx[-1]) - (n + 1) / n
    return float(g)


def shannon_entropy(p: np.ndarray, base: float = math.e) -> float:
    p = np.asarray(p, dtype=float)
    if p.sum() == 0:
        return 0.0
    p = p / p.sum()
    p = p[p > 0]
    return float(-np.sum(p * np.log(p)) / (np.log(base)))


THRESHOLDS = [0.0, 1e-4, 1e-3, 5e-3, 1e-2]


def analyze_df(df: pd.DataFrame) -> pd.DataFrame:
    # df: rows = samples, cols = cell types, values are proportions (rows sum ~1)
    rows = []
    for sid, row in df.iterrows():
        vals = row.values.astype(float)
        total = float(vals.sum()) if vals.size else 0.0
        if total == 0:
            vals_norm = vals
        else:
            vals_norm = vals / total

        counts = {f"n_gt_{t}": int((vals_norm > t).sum()) for t in THRESHOLDS}
        nonzero = int((vals_norm > 0).sum())
        eff_num = math.exp(shannon_entropy(vals_norm))
        ent = shannon_entropy(vals_norm)
        g = gini(vals_norm)
        dominant = float(vals_norm.max()) if vals_norm.size else 0.0
        sorted_vals = np.sort(vals_norm)[::-1] if vals_norm.size else np.array([])
        cum1 = float(sorted_vals[0]) if sorted_vals.size >= 1 else 0.0
        cum3 = float(sorted_vals[:3].sum()) if sorted_vals.size >= 3 else float(sorted_vals.sum())
        cum5 = float(sorted_vals[:5].sum()) if sorted_vals.size >= 5 else float(sorted_vals.sum())
        unresolved_mass = 0.0
        other_mass = 0.0
        # detect columns named 'Other' or unresolved_*
        for cname in df.columns:
            if cname.lower() == "other":
                other_mass = float(row[cname])
            if cname.startswith("unresolved_"):
                unresolved_mass += float(row[cname])

        frac_below_001 = float((vals_norm[(vals_norm > 0) & (vals_norm < 0.001)].sum()))

        rows.append({
            "sample": sid,
            "n_cell_types": int(vals.size),
            "raw_nonzero": nonzero,
            **counts,
            "effective_number": eff_num,
            "entropy": ent,
            "gini": g,
            "dominant_fraction": dominant,
            "top1_cum": cum1,
            "top3_cum": cum3,
            "top5_cum": cum5,
            "unresolved_mass": unresolved_mass,
            "other_mass": other_mass,
            "fraction_removed_by_0.001": frac_below_001,
            "total_mass": total,
        })
    return pd.DataFrame(rows).set_index("sample")


def run():
    pred_files = list(PRED_DIR.glob("*.tsv"))
    if not pred_files:
        print("No prediction TSVs found in:", PRED_DIR)
        return

    summaries: Dict[str, pd.DataFrame] = {}
    per_sample_rows: List[pd.DataFrame] = []
    for p in pred_files:
        method = p.stem
        try:
            df = pd.read_csv(p, sep="\t", index_col=0)
        except Exception:
            df = pd.read_csv(p, sep="\t", index_col=0, engine="python")
        # Ensure numeric
        df = df.fillna(0.0).astype(float)
        # If rows are samples or columns are samples, try to detect orientation.
        row_sums = df.sum(axis=1)
        col_sums = df.sum(axis=0)
        if (col_sums.abs().sum() > row_sums.abs().sum()):
            df = df.T

        # Normalize rows to sum to 1 if not already
        rsum = df.sum(axis=1)
        df = df.div(rsum.replace(0, 1), axis=0)

        summary = analyze_df(df)
        summaries[method] = summary
        per_sample = summary.copy()
        per_sample["method"] = method
        per_sample_rows.append(per_sample.reset_index())

    all_samples = pd.concat(per_sample_rows, ignore_index=True)
    all_samples.to_csv(OUT_DIR / "tcga_prediction_complexity_audit_by_sample.tsv", sep="\t", index=False)

    # aggregate per-method statistics
    method_stats = []
    for m, df in summaries.items():
        method_stats.append({
            "method": m,
            "n_samples": int(df.shape[0]),
            "mean_raw_nonzero": float(df["raw_nonzero"].mean()),
            "median_raw_nonzero": float(df["raw_nonzero"].median()),
            "mean_effective_number": float(df["effective_number"].mean()),
            "median_effective_number": float(df["effective_number"].median()),
            "mean_entropy": float(df["entropy"].mean()),
            "mean_gini": float(df["gini"].mean()),
            "mean_dominant_fraction": float(df["dominant_fraction"].mean()),
            "mean_top5_cum": float(df["top5_cum"].mean()),
        })
    pd.DataFrame(method_stats).to_csv(OUT_DIR / "tcga_prediction_complexity_audit_by_method.tsv", sep="\t", index=False)

    print("Wrote:")
    print(" -", OUT_DIR / "tcga_prediction_complexity_audit_by_sample.tsv")
    print(" -", OUT_DIR / "tcga_prediction_complexity_audit_by_method.tsv")


if __name__ == "__main__":
    run()
