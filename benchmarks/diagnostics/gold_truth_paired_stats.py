"""Paired bootstrap statistics for the gold-truth bulk benchmark.

Aligns per-sample TissueResolve predictions (saved by
``gold_truth_performance_benchmark.py``) with the external-tool predictions
(saved by ``gold_truth_external_benchmark.py``) on the SAME donor-held-out
``full_overlap`` mixtures, then computes paired bootstrap confidence intervals
for metric *differences*, bootstrap rank stability, and per-task ranking.

It refuses to run (and does NOT fabricate CIs) if the internal per-sample
predictions are missing. No estimates are modified here; this is read-only
scoring of already-saved predictions.
"""
from __future__ import annotations

import argparse
import sys
from pathlib import Path

import numpy as np
import pandas as pd

REPO = Path(__file__).resolve().parents[2]
sys.path.insert(0, str(REPO))
sys.path.insert(0, str(REPO / "src"))

from benchmarks.shared.stats import bootstrap_ci, rank_stability

PERF_DIR = REPO / "benchmarks" / "outputs" / "gold_truth_performance"
EXTERNAL_DIR = REPO / "benchmarks" / "outputs" / "gold_truth_external"
INPUT_DIR = REPO / "benchmarks" / "external_tools" / "inputs" / "gold_truth_bulk"
PRED_DIR = PERF_DIR / "predictions"

DATASETS = ["breast", "lung_hlca_subset"]

# Internal canonical method name -> saved prediction file stem.
INTERNAL_METHODS = [
    "TissueResolve_flat",
    "TissueResolve_hierarchical_soft",
    "TissueResolve_auto",
    "TissueResolve_broad_only",
    "TissueResolve_hard_legacy",
    "TissueResolve_ungated_diagnostic",
]
EXTERNAL_METHODS = ["NNLS_external_control", "MuSiC", "BisqueRNA"]

# Pairwise comparisons requested (A vs B); difference reported as A - B.
COMPARISONS = [
    ("TissueResolve_flat", "MuSiC"),
    ("TissueResolve_flat", "BisqueRNA"),
    ("TissueResolve_flat", "NNLS_external_control"),
    ("TissueResolve_hierarchical_soft", "MuSiC"),
    ("TissueResolve_auto", "MuSiC"),
]

# metric -> (higher_is_better)
METRIC_DIRECTION = {
    "broad_rmse": False, "broad_pearson": True,
    "fine_rmse": False, "fine_pearson": True,
    "conditional_rmse": False,
    "absent_subtype_mass": False,
    "rare_subtype_abs_error": False,
}


# ---------------------------------------------------------------------------
# Loading
# ---------------------------------------------------------------------------


def _truth(kind: str, dataset: str) -> pd.DataFrame:
    f = {"fine": "fine_truth.tsv", "broad": "broad_truth.tsv",
         "conditional": "conditional_truth.tsv"}[kind]
    df = pd.read_csv(PERF_DIR / f, sep="\t")
    df = df[(df["dataset"] == dataset) & (df["group"] == "full_overlap")].copy()
    df = df.drop(columns=[c for c in ("dataset", "group") if c in df.columns])
    df = df.set_index("sample")
    # drop label columns that are entirely NaN for this dataset (other dataset's labels)
    df = df.dropna(axis=1, how="all")
    return df.astype(float).fillna(0.0)


def _label_mapping(dataset: str) -> dict[str, str]:
    p = INPUT_DIR / dataset / "label_mapping.tsv"
    if not p.exists():
        return {}
    m = pd.read_csv(p, sep="\t")
    return dict(zip(m["fine_cell_type"].astype(str), m["broad_family"].astype(str)))


def _read_internal(method: str, dataset: str) -> pd.DataFrame | None:
    p = PRED_DIR / f"{method}__{dataset}.tsv"
    if not p.exists():
        return None
    df = pd.read_csv(p, sep="\t", index_col=0)
    df.index = df.index.map(str)
    df.columns = df.columns.map(str)
    return df


def _read_external(method: str, dataset: str) -> pd.DataFrame | None:
    p = EXTERNAL_DIR / "raw_predictions" / f"{method}__{dataset}.tsv"
    if not p.exists():
        return None
    df = pd.read_csv(p, sep="\t", index_col=0)
    df.index = df.index.map(str)
    df.columns = df.columns.map(str)
    df = df.clip(lower=0.0)
    rs = df.sum(axis=1).replace(0, np.nan)
    df = df.div(rs, axis=0).fillna(0.0)
    return df


def _fine_aligned(pred: pd.DataFrame, truth_fine: pd.DataFrame) -> pd.DataFrame:
    # drop unresolved_* columns, reindex to truth fine labels (matches internal scoring)
    cols = [c for c in pred.columns if not str(c).startswith("unresolved_")]
    return pred[cols].reindex(index=truth_fine.index,
                              columns=truth_fine.columns).fillna(0.0)


def _broad_from_fine(fine_pred: pd.DataFrame, mapping: dict[str, str],
                     truth_broad: pd.DataFrame) -> pd.DataFrame:
    out = pd.DataFrame(0.0, index=fine_pred.index, columns=truth_broad.columns)
    for fine_lbl in fine_pred.columns:
        fam = mapping.get(fine_lbl)
        if fam in out.columns:
            out[fam] = out[fam] + fine_pred[fine_lbl]
    return out


def _conditional_from_fine(fine_pred: pd.DataFrame, mapping: dict[str, str],
                           truth_cond: pd.DataFrame) -> pd.DataFrame:
    # within-family conditional proportions: fine / family-sum
    out = pd.DataFrame(0.0, index=fine_pred.index, columns=truth_cond.columns)
    fam_members: dict[str, list[str]] = {}
    for lbl in fine_pred.columns:
        fam_members.setdefault(mapping.get(lbl, lbl), []).append(lbl)
    for fam, members in fam_members.items():
        members = [m for m in members if m in fine_pred.columns]
        fam_sum = fine_pred[members].sum(axis=1).replace(0, np.nan)
        for m in members:
            if m in out.columns:
                out[m] = (fine_pred[m] / fam_sum).fillna(0.0)
    return out.reindex(columns=truth_cond.columns).fillna(0.0)


# ---------------------------------------------------------------------------
# Per-sample metrics
# ---------------------------------------------------------------------------


def _per_sample_rmse(pred: pd.DataFrame, truth: pd.DataFrame) -> np.ndarray:
    p, t = pred.to_numpy(float), truth.to_numpy(float)
    return np.sqrt(np.mean((p - t) ** 2, axis=1))


def _per_sample_pearson(pred: pd.DataFrame, truth: pd.DataFrame) -> np.ndarray:
    p, t = pred.to_numpy(float), truth.to_numpy(float)
    out = np.full(p.shape[0], np.nan)
    for i in range(p.shape[0]):
        if np.std(p[i]) > 1e-12 and np.std(t[i]) > 1e-12:
            out[i] = np.corrcoef(p[i], t[i])[0, 1]
    return out


def _per_sample_absent_mass(fine_pred: pd.DataFrame, fine_truth: pd.DataFrame) -> np.ndarray:
    absent = fine_truth.to_numpy(float) < 1e-9
    return np.where(absent, fine_pred.to_numpy(float), 0.0).sum(axis=1)


def _per_sample_rare_abs_error(fine_pred: pd.DataFrame, fine_truth: pd.DataFrame,
                               rare: str) -> np.ndarray:
    if rare not in fine_pred.columns or rare not in fine_truth.columns:
        return np.full(fine_pred.shape[0], np.nan)
    return np.abs(fine_pred[rare].to_numpy(float) - fine_truth[rare].to_numpy(float))


# ---------------------------------------------------------------------------
# Assemble per-sample metric matrices per dataset
# ---------------------------------------------------------------------------


def _collect(dataset: str, n_boot: int, seed: int):
    fine_t = _truth("fine", dataset)
    broad_t = _truth("broad", dataset)
    cond_t = _truth("conditional", dataset)
    mapping = _label_mapping(dataset)
    rare = list(fine_t.columns)[-1]  # internal benchmark uses the last fine label as rare

    methods = {}
    for m in INTERNAL_METHODS:
        raw = _read_internal(m, dataset)
        if raw is not None:
            methods[m] = raw
    for m in EXTERNAL_METHODS:
        raw = _read_external(m, dataset)
        if raw is not None:
            methods[m] = raw

    per_sample = {}  # method -> {metric: np.ndarray over samples}
    for m, raw in methods.items():
        fine_p = _fine_aligned(raw, fine_t)
        broad_p = _broad_from_fine(fine_p, mapping, broad_t)
        cond_p = _conditional_from_fine(fine_p, mapping, cond_t)
        per_sample[m] = {
            "fine_rmse": _per_sample_rmse(fine_p, fine_t),
            "fine_pearson": _per_sample_pearson(fine_p, fine_t),
            "broad_rmse": _per_sample_rmse(broad_p, broad_t),
            "broad_pearson": _per_sample_pearson(broad_p, broad_t),
            "conditional_rmse": _per_sample_rmse(cond_p, cond_t),
            "absent_subtype_mass": _per_sample_absent_mass(fine_p, fine_t),
            "rare_subtype_abs_error": _per_sample_rare_abs_error(fine_p, fine_t, rare),
        }
    return per_sample, list(fine_t.index)


def run(n_boot: int, seed: int) -> int:
    EXTERNAL_DIR.mkdir(parents=True, exist_ok=True)

    # refuse to run if internal predictions are missing
    missing = [m for m in ("TissueResolve_flat",)
               if not (PRED_DIR / f"{m}__{DATASETS[0]}.tsv").exists()]
    if missing:
        print(f"[abort] internal per-sample predictions missing ({missing}). "
              f"Re-run gold_truth_performance_benchmark.py with --save-predictions. "
              f"Refusing to fabricate CIs.")
        return 2

    ci_rows, rank_task_rows = [], []
    rank_frames = []

    for dataset in DATASETS:
        per_sample, samples = _collect(dataset, n_boot, seed)
        if not per_sample:
            print(f"[skip] {dataset}: no predictions found")
            continue

        # ---- paired bootstrap CIs of metric differences ----
        for a, b in COMPARISONS:
            if a not in per_sample or b not in per_sample:
                continue
            for metric in METRIC_DIRECTION:
                va = per_sample[a][metric]
                vb = per_sample[b][metric]
                mask = np.isfinite(va) & np.isfinite(vb)
                diff = (va - vb)[mask]
                if diff.size == 0:
                    continue
                ci = bootstrap_ci(diff, statistic=np.mean, n_boot=n_boot,
                                  seed=seed)
                higher_better = METRIC_DIRECTION[metric]
                # favors A if difference sign is in A's favor and CI excludes 0
                excludes_zero = (ci["lo"] > 0) or (ci["hi"] < 0)
                if not excludes_zero:
                    favors = "tie"
                else:
                    a_better = (ci["point"] > 0) == higher_better
                    favors = a if a_better else b
                ci_rows.append({
                    "dataset": dataset, "comparison": f"{a}_vs_{b}",
                    "method_a": a, "method_b": b, "metric": metric,
                    "mean_diff_a_minus_b": ci["point"], "ci_lo": ci["lo"],
                    "ci_hi": ci["hi"], "n_samples": ci["n"],
                    "higher_is_better": higher_better,
                    "ci_excludes_zero": excludes_zero, "favors": favors,
                })

        # ---- rank stability over per-sample fine RMSE ----
        fine_rmse = pd.DataFrame({m: per_sample[m]["fine_rmse"] for m in per_sample},
                                 index=samples)
        rs = rank_stability(fine_rmse, higher_is_better=False, n_boot=n_boot, seed=seed)
        rs = rs.reset_index().assign(dataset=dataset, metric="fine_rmse")
        rank_frames.append(rs)

        # ---- per-task ranking (aggregate mean per method) ----
        for metric, higher in METRIC_DIRECTION.items():
            def _safe_mean(arr):
                arr = np.asarray(arr, float)
                return float(np.nanmean(arr)) if np.isfinite(arr).any() else float("nan")
            agg = {m: _safe_mean(per_sample[m][metric]) for m in per_sample}
            ser = pd.Series(agg).sort_values(ascending=not higher)
            for rank, (m, val) in enumerate(ser.items(), start=1):
                rank_task_rows.append({
                    "dataset": dataset, "task": metric, "method": m,
                    "value": float(val), "rank": rank,
                    "higher_is_better": higher,
                })

    pd.DataFrame(ci_rows).to_csv(EXTERNAL_DIR / "paired_bootstrap_ci.tsv", sep="\t", index=False)
    if rank_frames:
        pd.concat(rank_frames, ignore_index=True).to_csv(
            EXTERNAL_DIR / "rank_stability.tsv", sep="\t", index=False)
    pd.DataFrame(rank_task_rows).to_csv(EXTERNAL_DIR / "rank_by_task.tsv", sep="\t", index=False)
    print(f"[done] wrote paired_bootstrap_ci.tsv, rank_stability.tsv, rank_by_task.tsv "
          f"under {EXTERNAL_DIR}")
    return 0


def parse_args(argv=None):
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument("--n-boot", type=int, default=5000)
    ap.add_argument("--seed", type=int, default=2026)
    return ap.parse_args(argv)


def main(argv=None) -> int:
    args = parse_args(argv)
    return run(args.n_boot, args.seed)


if __name__ == "__main__":
    raise SystemExit(main())
