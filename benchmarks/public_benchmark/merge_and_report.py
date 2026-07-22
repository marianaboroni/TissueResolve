#!/usr/bin/env python
"""Public benchmark — merge + re-score + aggregate + figures (bulk).

Re-scores EVERY saved prediction (in-process TissueResolve/baselines + external MuSiC/Bisque/
BayesPrism) with the SAME fair scoring (broad credits unresolved mass; fine/conditional/swap/rare/
absent/identifiability-stratified). Produces the prescribed output tree + figures with source data +
paired comparisons vs TissueResolve_PoissonGLM. Honest: methods without predictions are recorded in
method_status, never fabricated.

Usage: PYTHONPATH=src:. python benchmarks/public_benchmark/merge_and_report.py --run-real-data
"""
from __future__ import annotations

import argparse
import json
import sys
import warnings
from pathlib import Path

import numpy as np
import pandas as pd

REPO = Path(__file__).resolve().parents[2]
sys.path.insert(0, str(REPO / "examples" / "real_breast_cancer" / "scripts"))
sys.path.insert(0, str(REPO)); sys.path.insert(0, str(REPO / "src"))

from benchmarks.public_benchmark import scoring as SC  # noqa: E402

BULK = REPO / "benchmarks" / "results" / "public_benchmark" / "bulk"
EXPORT = BULK / "external_inputs"
COMBINED = REPO / "benchmarks" / "results" / "public_benchmark" / "combined"
FIGSRC = COMBINED / "figure_source_data"
FIGDIR = REPO / "docs" / "figures" / "public_benchmark"
SCEN = ["base", "low_depth", "reduced_overlap"]
REF_METHOD = "TissueResolve_PoissonGLM"


def _mapping_and_cert(dataset):
    """Fine→broad mapping + saved per-scenario certificate (class + clusters) — no recompute."""
    from benchmarks.spatial.run_weak_smoothing_grid import _load_dataset
    D = _load_dataset(dataset)
    cert = pd.read_csv(BULK / "identifiability_certificate.tsv", sep="\t")
    cert = cert[cert.dataset == dataset]
    types = sorted(cert.cell_type.astype(str).unique().tolist())
    mp = {t: str(D["mapping"].get(t, t)) for t in types}
    return mp, types, cert


def score_bulk(datasets):
    rows, ct_rows = [], []
    for ds in datasets:
        mp, types, cert_tbl = _mapping_and_cert(ds)
        for scen in SCEN:
            sd = EXPORT / ds / scen
            tf = sd / "truth_mrna.tsv"
            if not tf.exists():
                continue
            truth = pd.read_csv(tf, sep="\t", index_col=0).reindex(columns=types).fillna(0.0)
            sub = cert_tbl[cert_tbl.scenario == scen]
            class_by = dict(zip(sub.cell_type.astype(str), sub.recoverability.astype(str)))
            clusters = [str(c).split(";") for c in sub.get("cluster", pd.Series(dtype=str)).dropna().unique()
                        if str(c) and str(c) != "nan"]
            for pf in sorted(sd.glob("pred_*.tsv")):
                method = pf.name[len("pred_"):-len(".tsv")]
                pred = pd.read_csv(pf, sep="\t", index_col=0)
                pred.index = [str(i) for i in pred.index]
                # external tools may return samples×types already aligned to truth.index
                if not set(truth.index) <= set(pred.index):
                    if pred.shape[0] == truth.shape[0]:
                        pred.index = truth.index
                pred = pred.reindex(index=truth.index)
                keep = [c for c in pred.columns if str(c).startswith("unresolved_") or c in types]
                pred = pred[keep].fillna(0.0)
                try:
                    m = SC.score_all(truth, pred, mp, clusters, class_by)
                except Exception as e:
                    print(f"score fail {ds}/{scen}/{method}: {e}"); continue
                rows.append(dict(dataset=ds, scenario=scen, method=method, n_samples=len(truth), **m))
                for cr in SC.per_celltype_rows(truth, pred.reindex(columns=types).fillna(0.0), class_by, mp):
                    ct_rows.append(dict(dataset=ds, scenario=scen, method=method, **cr))
    return pd.DataFrame(rows), pd.DataFrame(ct_rows)


def paired(df, metrics):
    rng = np.random.default_rng(0); out = []
    piv = {mtr: df.pivot_table(index=["dataset", "scenario"], columns="method", values=mtr) for mtr in metrics}
    for mtr in metrics:
        P = piv[mtr]
        if REF_METHOD not in P.columns:
            continue
        for method in P.columns:
            if method == REF_METHOD:
                continue
            a = P[REF_METHOD].to_numpy(float); b = P[method].to_numpy(float)
            ok = np.isfinite(a) & np.isfinite(b); a, b = a[ok], b[ok]
            if len(a) < 2:
                continue
            diff = b - a  # >0 => method worse than TR Poisson on an error metric
            boot = [np.mean(diff[rng.integers(0, len(diff), len(diff))]) for _ in range(2000)]
            out.append(dict(metric=mtr, method=method, vs=REF_METHOD, mean_diff=float(np.mean(diff)),
                            ci_low=float(np.percentile(boot, 2.5)), ci_high=float(np.percentile(boot, 97.5)),
                            tr_poisson_better_rate=float(np.mean(diff > 0)), n=int(len(diff))))
    return pd.DataFrame(out)


def _bar(df, value, title, path, ascending=True):
    import matplotlib; matplotlib.use("Agg"); import matplotlib.pyplot as plt
    s = df.groupby("method")[value].mean().sort_values(ascending=ascending)
    s.to_csv(str(path) + ".data.tsv", sep="\t")
    fig, ax = plt.subplots(figsize=(7, 0.45 * len(s) + 1))
    colors = ["#2c7fb8" if "TissueResolve" in m else "#bdbdbd" for m in s.index]
    ax.barh(range(len(s)), s.values, color=colors); ax.set_yticks(range(len(s))); ax.set_yticklabels(s.index)
    ax.set_xlabel(value); ax.set_title(title); fig.tight_layout()
    for ext in ("png", "pdf", "svg"):
        fig.savefig(f"{path}.{ext}", dpi=140)
    plt.close(fig)


def main(argv=None):
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument("--run-real-data", action="store_true")
    ap.add_argument("--datasets", nargs="+", default=["breast", "lung"])
    ap.parse_args(argv)
    for d in (COMBINED, FIGSRC, FIGDIR):
        d.mkdir(parents=True, exist_ok=True)
    ns = argparse.Namespace()
    df, ctdf = score_bulk(["breast", "lung"])
    df.to_csv(BULK / "metrics_by_scenario.tsv", sep="\t", index=False)
    ctdf.to_csv(BULK / "metrics_by_celltype.tsv", sep="\t", index=False)

    def _sub(cols, name):
        c = [x for x in cols if x in df.columns]
        df.groupby("method")[c].mean().round(4).to_csv(BULK / name, sep="\t")
    _sub(["broad_rmse", "broad_mae", "broad_pearson", "broad_spearman", "broad_bias",
          "broad_jsd", "broad_dominant_acc"], "broad_metrics.tsv")
    _sub(["fine_rmse", "fine_mae", "fine_pearson", "fine_spearman", "fine_bias", "fine_jsd",
          "fine_dominant_acc", "unresolved_mass"], "fine_metrics.tsv")
    _sub(["condfam_rmse", "condfam_pearson", "swap_rmse", "swap_pearson"],
         "conditional_within_family_metrics.tsv")
    _sub(["rare_recall", "rare_precision", "rare_fpr", "detection_recall"], "rare_detection_metrics.tsv")
    _sub(["absent_mass", "false_positive_subtype_rate", "eff_n_pred", "eff_n_true"], "absent_mass_metrics.tsv")
    _sub(["err_RESOLVABLE", "err_WEAKLY_RESOLVABLE", "err_GROUP_ONLY", "err_UNRESOLVABLE",
          "overconfident_error_rate", "unresolved_mass"], "identifiability_stratified_metrics.tsv")

    pc = paired(df, ["broad_rmse", "fine_rmse", "condfam_rmse", "swap_rmse"])
    pc.to_csv(BULK / "paired_comparisons.tsv", sep="\t", index=False)

    # method status = in-process + external
    st = []
    ip = BULK / "method_status.tsv"
    if ip.exists():
        st.append(pd.read_csv(ip, sep="\t"))
    ext = EXPORT / "external_method_status.tsv"
    if ext.exists():
        e = pd.read_csv(ext, sep="\t"); st.append(e.rename(columns={"runtime_seconds": "runtime_s"}))
    status = pd.concat(st, ignore_index=True) if st else pd.DataFrame()
    status.to_csv(BULK / "method_status_full.tsv", sep="\t", index=False)

    # figures
    _bar(df, "broad_pearson", "Bulk broad-family accuracy (Pearson, higher=better)",
         FIGDIR / "bulk_broad_accuracy", ascending=True)
    _bar(df, "fine_pearson", "Bulk fine-subpopulation accuracy (Pearson, higher=better)",
         FIGDIR / "bulk_fine_accuracy", ascending=True)
    _bar(df, "condfam_rmse", "Bulk conditional within-family RMSE (lower=better)",
         FIGDIR / "bulk_conditional_rmse", ascending=False)
    _bar(df, "rare_recall", "Bulk rare-population recall (higher=better)",
         FIGDIR / "bulk_rare_detection", ascending=True)
    _bar(df, "overconfident_error_rate", "Overconfident fine-state error in GROUP_ONLY/UNRESOLVABLE (lower=better)",
         FIGDIR / "bulk_overconfident_error", ascending=False)
    _bar(df, "runtime_s", "Bulk runtime per scenario (s)", FIGDIR / "runtime_comparison", ascending=True) if "runtime_s" in df.columns else None

    summ = df.groupby("method")[[c for c in ("broad_pearson", "broad_rmse", "fine_pearson", "fine_rmse",
        "condfam_rmse", "swap_rmse", "rare_recall", "absent_mass", "overconfident_error_rate",
        "unresolved_mass") if c in df.columns]].mean().round(4)
    summ.to_csv(COMBINED / "benchmark_summary.tsv", sep="\t")
    (COMBINED / "run_manifest.json").write_text(json.dumps(
        {"methods": sorted(df.method.unique().tolist()), "datasets": ["breast", "lung"],
         "scenarios": SCEN, "n_rows": int(len(df))}, indent=2))
    print("=== BENCHMARK SUMMARY (mean over datasets×scenarios) ===\n", summ.to_string())
    print("\n=== PAIRED vs TR_PoissonGLM (mean_diff>0 => TR better on error metric) ===\n",
          pc.to_string(index=False) if len(pc) else "(none)")
    print("\n=== METHOD STATUS ===\n", status.to_string(index=False) if len(status) else "(none)")
    print(f"\nWrote {COMBINED} and {FIGDIR}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
