#!/usr/bin/env python
"""Stage 1.6 benchmark: separability-aware adaptive budget (D) vs fixed-per-type (C) vs
current GeneSelector (A), donor-disjoint, Poisson, paired-by-mixture. Answers §30.

Strategies (same mixtures, same solver):
  A_current        GeneSelector (n_genes=500)
  C_fixed_per_type fine_global at top_n_per_type=15 (recommended 1.5 baseline)
  D_adaptive       SeparabilityAwareBudgetAllocator (per-type, difficulty-driven)
  F_all_eligible   union of all donor-aware DE candidates (control: "just use many genes")

Usage: PYTHONPATH=src:. python benchmarks/signatures/run_stage1_6.py --run-real-data \
       --datasets breast lung --seeds 0 1 2 3 4
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

OUT = REPO / "benchmarks" / "results" / "signatures" / "stage1_6"


def _cond(truth, pred, mp):
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


def _panels(train, ref, ctc, dc, mp, adaptive_kw):
    from tissueresolve.reference.markers import GeneSelector
    from tissueresolve.reference import gene_selection as GS
    from tissueresolve.reference.adaptive_signature import SeparabilityAwareBudgetAllocator
    from tissueresolve.config import GeneConfig
    t0 = time.perf_counter()
    with warnings.catch_warnings():
        warnings.simplefilter("ignore")
        A = GeneSelector(config=GeneConfig(n_genes=500)).select(ref).selected_genes
        C = GS.select_donor_aware_genes(train, ctc, dc, top_n_per_type=15)
        adapt = SeparabilityAwareBudgetAllocator(ctc, dc, mp, **adaptive_kw).allocate(train)
        Dg = adapt.genes
        types = [t for t in pd.unique(train.obs[ctc].astype(str))
                 if (train.obs[ctc].astype(str) == t).sum() >= 10]
        F = list(dict.fromkeys(g for t in types
                 for g in GS.donor_aware_de(train, ctc, dc, t, mode="one_vs_rest",
                                            top_n=100, min_cells=10).genes))
    build_s = time.perf_counter() - t0
    return {"A_current": A, "C_fixed_per_type": C, "D_adaptive": Dg, "F_all_eligible": F}, adapt, build_s


def run(dataset, seeds, adaptive_kw, rows_mix, rows_type, meta):
    from tissueresolve.solver import PoissonGLMSolver
    D = _load_dataset(dataset)
    ad, ctc, dc, ref, types = (D["adata"], D["celltype_col"], D["donor_col"], D["ref"], D["ref_types"])
    mp = {t: D["mapping"].get(t, t) for t in types}
    V.validate_hierarchy_covers_dataset(mp, types, dataset=dataset)
    V.validate_donor_disjoint(set(D["ref_donors"]), set(D["query_donors"]))
    V.validate_selection_donors(set(D["ref_donors"]), set(D["query_donors"]))
    train = ad[ad.obs[dc].astype(str).isin(set(D["ref_donors"])).to_numpy()].copy()
    panels, adapt, build_s = _panels(train, ref, ctc, dc, mp, adaptive_kw)
    for strat, g in panels.items():
        V.validate_budget_recorded({"strategy": strat, "n_genes": len(g)})
    meta.append(dict(dataset=dataset, **{f"n_{k}": len(v) for k, v in panels.items()},
                     adaptive_build_s=round(build_s, 1),
                     adaptive_status_counts=json.dumps(adapt.manifest.get("status_counts", {}))))
    # per-type adaptive gene counts + status
    for t, v in adapt.per_type.items():
        rows_type.append(dict(dataset=dataset, cell_type=t, adaptive_genes=v["gene_count"],
                              status=v["status"], best_rmse=v.get("best_rmse"),
                              separability_bc=v.get("separability_bc"),
                              closest_confounder=v.get("closest_confounder")))
    for seed in seeds:
        counts, truth, _ = H.generate_pseudobulk(
            ad[ad.obs[dc].astype(str).isin(set(D["query_donors"])).to_numpy()].copy(),
            ctc, n_per_regime=2, n_cells=200, seed=seed)
        V.validate_matrix_orientation(counts, ref.gene_names)
        for strat, g in panels.items():
            with warnings.catch_warnings():
                warnings.simplefilter("ignore")
                props = PoissonGLMSolver(genes=g).solve(counts, ref).proportions
            V.validate_truth_pred_alignment(truth, props)
            cols = list(truth.columns)
            for s in truth.index:
                t = truth.loc[s]; p = props.reindex(columns=cols).fillna(0).loc[s] if s in props.index else pd.Series(0.0, index=cols)
                tv = t.to_numpy(float); pv = p.to_numpy(float); absent = tv < 1e-9
                rows_mix.append(dict(dataset=dataset, seed=seed, mixture=f"{seed}:{s}", strategy=strat,
                                     n_genes=len(g), RMSE=float(np.sqrt(np.mean((pv - tv) ** 2))),
                                     Pearson=float(np.corrcoef(pv, tv)[0, 1]) if pv.std() > 0 else np.nan,
                                     cond_RMSE=_cond(truth.loc[[s]], props.reindex([s]).fillna(0), mp),
                                     spillover_absent=float(pv[absent].sum())))


def _paired(mix):
    rng = np.random.default_rng(0); out = []
    for ds in mix.dataset.unique():
        d = mix[mix.dataset == ds]
        piv = d.pivot_table(index=["seed", "mixture"], columns="strategy", values=["RMSE", "cond_RMSE"])
        for base in ("C_fixed_per_type", "A_current"):
            for metric in ("RMSE", "cond_RMSE"):
                if "D_adaptive" not in piv[metric] or base not in piv[metric]:
                    continue
                a = piv[metric][base].to_numpy(float); b = piv[metric]["D_adaptive"].to_numpy(float)
                ok = np.isfinite(a) & np.isfinite(b); a, b = a[ok], b[ok]
                if len(a) < 3:
                    continue
                diff = b - a
                boot = [np.mean(diff[rng.integers(0, len(diff), len(diff))]) for _ in range(2000)]
                out.append(dict(dataset=ds, comparison=f"D_adaptive-vs-{base}", metric=metric,
                                mean_diff=float(np.mean(diff)), ci_low=float(np.percentile(boot, 2.5)),
                                ci_high=float(np.percentile(boot, 97.5)),
                                win_rate=float(np.mean(diff < 0)), n=int(len(diff))))
    return pd.DataFrame(out)


def main(argv=None):
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument("--run-real-data", action="store_true")
    ap.add_argument("--datasets", nargs="+", default=["breast", "lung"])
    ap.add_argument("--seeds", type=int, nargs="+", default=[0, 1, 2, 3, 4])
    args = ap.parse_args(argv)
    if not H.real_data_enabled(args.run_real_data):
        print(f"Refusing without --run-real-data (or {H.REAL_DATA_ENV}=1).", file=sys.stderr)
        return 2
    OUT.mkdir(parents=True, exist_ok=True)
    adaptive_kw = dict(sizes=(5, 10, 15, 20, 30, 40, 60), min_cells=10, perf_tol=0.05,
                       easy_rmse=0.10, unresolvable_rmse=0.30, max_genes=60, seed=0)
    rows_mix, rows_type, meta = [], [], []
    t0 = time.perf_counter()
    for ds in args.datasets:
        run(ds, args.seeds, adaptive_kw, rows_mix, rows_type, meta)
    mix = pd.DataFrame(rows_mix)
    mix.to_csv(OUT / "metrics_by_mixture.tsv", sep="\t", index=False)
    pd.DataFrame(rows_type).to_csv(OUT / "adaptive_gene_budget_by_celltype.tsv", sep="\t", index=False)
    pd.DataFrame(meta).to_csv(OUT / "budget_strategy_comparison.tsv", sep="\t", index=False)
    summ = mix.groupby(["dataset", "strategy"])[["RMSE", "Pearson", "cond_RMSE", "spillover_absent", "n_genes"]].mean().round(4)
    summ.to_csv(OUT / "summary_metrics.tsv", sep="\t")
    paired = _paired(mix); paired.to_csv(OUT / "paired_comparisons.tsv", sep="\t", index=False)
    (OUT / "benchmark_manifest.json").write_text(json.dumps(
        {"datasets": args.datasets, "seeds": args.seeds, "adaptive_kw": adaptive_kw,
         "runtime_s": round(time.perf_counter() - t0, 1)}, indent=2, default=str))
    print("\n=== SUMMARY ===\n", summ.to_string())
    print("\n=== PER-TYPE ADAPTIVE STATUS COUNTS ===")
    print(pd.DataFrame(rows_type).groupby(["dataset", "status"]).size().to_string())
    print("\n=== PAIRED (D_adaptive − baseline; negative ⇒ adaptive better) ===\n", paired.to_string(index=False))
    print("\n=== gene counts per strategy ===\n", pd.DataFrame(meta).to_string(index=False))
    print(f"\nWrote {OUT}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
