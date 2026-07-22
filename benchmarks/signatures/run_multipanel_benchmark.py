#!/usr/bin/env python
"""Permanent, validated multipanel signature benchmark (Stage 1.5, real-data).

Replaces the /tmp scripts. Donor-disjoint, leakage-guarded (validation.py aborts on
hierarchy/dataset mismatch, donor overlap, unrecorded budget, mis-orientation), equal
STRATIFIED gene budget across strategies, Poisson solver, paired by mixture. Writes a
manifest + per-mixture / per-cell-type / summary metrics under a versioned output dir.

Strategies this runner compares (extend as Strategy C/D/E mature):
  A_current       GeneSelector markers (baseline) at the same budget
  B_fine_global   optimizer fine_global panel (stratified per-type budget)

Usage:
  PYTHONPATH=src:. python benchmarks/signatures/run_multipanel_benchmark.py \
      --run-real-data --datasets breast lung --seeds 0 1 2 3 4 --budget 300
"""
from __future__ import annotations

import argparse
import json
import sys
import time
import warnings
from pathlib import Path

import numpy as np
import pandas as pd

REPO = Path(__file__).resolve().parents[2]
sys.path.insert(0, str(REPO / "examples" / "real_breast_cancer" / "scripts"))
sys.path.insert(0, str(REPO))
sys.path.insert(0, str(REPO / "src"))
import _harness as H  # noqa: E402
from benchmarks.spatial.run_weak_smoothing_grid import _load_dataset  # noqa: E402
from benchmarks.signatures import validation as V  # noqa: E402

OUT = REPO / "benchmarks" / "results" / "signatures"


def _git(*a):
    import subprocess
    try:
        return subprocess.run(["git", *a], cwd=str(REPO), capture_output=True,
                              text=True, timeout=15).stdout.strip()
    except Exception:
        return ""


def _cond_rmse(truth, pred, mp):
    fam = {}
    for c in truth.columns:
        fam.setdefault(str(mp.get(c, c)), []).append(c)
    multi = [c for m in fam.values() if len(m) > 1 for c in m]
    if not multi:
        return np.nan

    def cd(df):
        o = pd.DataFrame(0.0, index=df.index, columns=truth.columns)
        for m in fam.values():
            m = [x for x in m if x in df.columns]
            if len(m) < 2:
                continue
            s = df[m].sum(1).replace(0, np.nan)
            for x in m:
                o[x] = (df[x] / s).fillna(0)
        return o
    tc, pc = cd(truth), cd(pred.reindex(columns=truth.columns).fillna(0))
    return float(np.sqrt(np.mean((tc[multi].to_numpy(float) - pc[multi].to_numpy(float)) ** 2)))


def _score_mixture(truth_row, pred_row, cols):
    t = truth_row[cols].to_numpy(float); p = pred_row.reindex(cols).fillna(0).to_numpy(float)
    absent = t < 1e-9; present = ~absent; rare = present & (t < 0.05); pos = p > 0.01
    fp = int((pos & absent).sum()); tn = int((~pos & absent).sum())
    return dict(RMSE=float(np.sqrt(np.mean((p - t) ** 2))), MAE=float(np.mean(np.abs(p - t))),
                rare_fpr=(fp / (fp + tn) if fp + tn else np.nan),
                rare_recall=(int((pos & rare).sum()) / int(rare.sum()) if rare.sum() else np.nan))


def _strategies(train, ref, ctc, dc, mp, budget):
    from tissueresolve.reference.markers import GeneSelector
    from tissueresolve.reference.signature_optimizer import ReferenceSignatureOptimizer
    from tissueresolve.config import GeneConfig
    with warnings.catch_warnings():
        warnings.simplefilter("ignore")
        cur = GeneSelector(config=GeneConfig(n_genes=budget)).select(ref).selected_genes
        m = ReferenceSignatureOptimizer(ctc, dc, mapping=mp, gene_budget=budget,
                                        minimal_sizes=(50, 100, 150), seed=0).optimize(train)
    return {"A_current": (cur, {"strategy": "A_current", "n_genes": len(cur)}),
            "B_fine_global": (m.fine_global_genes,
                              {"strategy": "B_fine_global", "n_genes": len(m.fine_global_genes),
                               "status": m.status})}, m


def run_dataset(dataset, seeds, budget, rows_mix, rows_ct, strat_meta):
    from tissueresolve.solver import PoissonGLMSolver
    D = _load_dataset(dataset)
    ad, ctc, dc, ref, types = (D["adata"], D["celltype_col"], D["donor_col"],
                               D["ref"], D["ref_types"])
    mp = {t: D["mapping"].get(t, t) for t in types}
    # --- leakage / consistency guards (abort on violation) ---
    V.validate_hierarchy_covers_dataset(mp, types, dataset=dataset)
    ref_donors = set(D["ref_donors"]); query_donors = set(D["query_donors"])
    V.validate_donor_disjoint(ref_donors, query_donors)
    V.validate_selection_donors(ref_donors, query_donors)

    train = ad[ad.obs[dc].astype(str).isin(ref_donors).to_numpy()].copy()
    strategies, _model = _strategies(train, ref, ctc, dc, mp, budget)
    for _, meta in strategies.values():
        V.validate_budget_recorded(meta); meta["dataset"] = dataset; strat_meta.append(meta)

    for seed in seeds:
        counts, truth, _ = H.generate_pseudobulk(
            ad[ad.obs[dc].astype(str).isin(query_donors).to_numpy()].copy(),
            ctc, n_per_regime=2, n_cells=200, seed=seed)
        V.validate_matrix_orientation(counts, ref.gene_names)
        cols = list(truth.columns)
        preds = {}
        for name, (genes, _) in strategies.items():
            with warnings.catch_warnings():
                warnings.simplefilter("ignore")
                preds[name] = PoissonGLMSolver(genes=genes).solve(counts, ref).proportions
        for name, pred in preds.items():
            V.validate_truth_pred_alignment(truth, pred)
            for s in truth.index:                      # per-mixture (paired unit)
                sc = _score_mixture(truth.loc[s], pred.loc[s] if s in pred.index else
                                    pd.Series(0.0, index=cols), cols)
                sc.update(dataset=dataset, seed=seed, mixture=s, strategy=name,
                          cond_RMSE=_cond_rmse(truth.loc[[s]], pred.reindex([s]).fillna(0), mp))
                rows_mix.append(sc)
            # per cell type (mean abs error across mixtures)
            pa = pred.reindex(index=truth.index, columns=cols).fillna(0)
            for c in cols:
                rows_ct.append(dict(dataset=dataset, seed=seed, strategy=name, cell_type=c,
                                    mae=float(np.mean(np.abs(pa[c].to_numpy(float) - truth[c].to_numpy(float))))))


def _paired(df_mix):
    """Paired-by-mixture bootstrap: B vs A per dataset for RMSE and cond_RMSE."""
    out = []
    rng = np.random.default_rng(0)
    for ds in df_mix["dataset"].unique():
        d = df_mix[df_mix.dataset == ds]
        piv = d.pivot_table(index=["seed", "mixture"], columns="strategy",
                            values=["RMSE", "cond_RMSE"])
        for metric in ("RMSE", "cond_RMSE"):
            if ("A_current" not in piv[metric]) or ("B_fine_global" not in piv[metric]):
                continue
            a = piv[metric]["A_current"].to_numpy(float); b = piv[metric]["B_fine_global"].to_numpy(float)
            ok = np.isfinite(a) & np.isfinite(b); a, b = a[ok], b[ok]
            if len(a) < 3:
                continue
            diff = b - a                      # negative => B better (lower is better)
            boot = [np.mean(diff[rng.integers(0, len(diff), len(diff))]) for _ in range(2000)]
            out.append(dict(dataset=ds, metric=metric, mean_diff_BminusA=float(np.mean(diff)),
                            ci_low=float(np.percentile(boot, 2.5)), ci_high=float(np.percentile(boot, 97.5)),
                            frac_B_better=float(np.mean(diff < 0)), n_mixtures=int(len(diff))))
    return pd.DataFrame(out)


def main(argv=None):
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument("--run-real-data", action="store_true")
    ap.add_argument("--datasets", nargs="+", default=["breast", "lung"])
    ap.add_argument("--seeds", type=int, nargs="+", default=[0, 1, 2, 3, 4])
    ap.add_argument("--budget", type=int, default=300)
    args = ap.parse_args(argv)
    if not H.real_data_enabled(args.run_real_data):
        print(f"Refusing without --run-real-data (or {H.REAL_DATA_ENV}=1).", file=sys.stderr)
        return 2
    OUT.mkdir(parents=True, exist_ok=True)
    t0 = time.perf_counter()
    rows_mix, rows_ct, strat_meta = [], [], []
    for ds in args.datasets:
        run_dataset(ds, args.seeds, args.budget, rows_mix, rows_ct, strat_meta)
    mix = pd.DataFrame(rows_mix); ct = pd.DataFrame(rows_ct)
    mix.to_csv(OUT / "metrics_by_mixture.tsv", sep="\t", index=False)
    ct.to_csv(OUT / "metrics_by_celltype.tsv", sep="\t", index=False)
    pd.DataFrame(strat_meta).to_csv(OUT / "strategy_manifest.tsv", sep="\t", index=False)
    summ = mix.groupby(["dataset", "strategy"])[["RMSE", "MAE", "cond_RMSE", "rare_fpr", "rare_recall"]].mean().round(4)
    summ.to_csv(OUT / "summary_metrics.tsv", sep="\t")
    paired = _paired(mix); paired.to_csv(OUT / "paired_comparisons.tsv", sep="\t", index=False)
    (OUT / "benchmark_manifest.json").write_text(json.dumps({
        "datasets": args.datasets, "seeds": args.seeds, "budget": args.budget,
        "budget_strategy": "stratified_round_robin", "solver": "poisson",
        "git_commit": _git("rev-parse", "HEAD"), "git_dirty": bool(_git("status", "--porcelain")),
        "runtime_s": round(time.perf_counter() - t0, 1)}, indent=2))
    print("\n=== SUMMARY (stratified budget) ==="); print(summ.to_string())
    print("\n=== PAIRED (B_fine_global − A_current; negative => B better) ===")
    print(paired.to_string(index=False))
    print(f"\nWrote {OUT}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
