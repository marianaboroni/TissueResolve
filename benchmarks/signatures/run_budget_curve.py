#!/usr/bin/env python
"""Gene-budget calibration curves (Stage 1.5 §7). Tests the 'lung needs ~550 genes'
hypothesis instead of assuming it. Donor-disjoint, Poisson, real-data.

Efficient: per-type one-vs-rest donor-aware DE is computed ONCE (large top_n), then the
top_n_per_type sweep just slices the cached per-type rankings (no repeated DE). Current
GeneSelector is evaluated at matched total budgets. Outputs budget_curve*.tsv + pareto_front.tsv.
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
sys.path.insert(0, str(REPO))
sys.path.insert(0, str(REPO / "src"))
import _harness as H  # noqa: E402
from benchmarks.spatial.run_weak_smoothing_grid import _load_dataset  # noqa: E402

OUT = REPO / "benchmarks" / "results" / "signatures"
TOP_N = [5, 10, 15, 20, 30, 40]
TOTAL_BUDGETS = [100, 200, 300, 500, 750, 1000]


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


def _score(truth, props, mp):
    p = props.reindex(index=truth.index, columns=truth.columns).fillna(0)
    d = p.to_numpy(float) - truth.to_numpy(float)
    return dict(RMSE=float(np.sqrt((d ** 2).mean())),
                Pearson=float(np.corrcoef(p.to_numpy(float).ravel(), truth.to_numpy(float).ravel())[0, 1]),
                cond_RMSE=_cond_rmse(truth, props, mp))


def run(dataset, seeds, rows):
    from tissueresolve.solver import PoissonGLMSolver
    from tissueresolve.reference.markers import GeneSelector
    from tissueresolve.reference import gene_selection as GS
    from tissueresolve.config import GeneConfig
    D = _load_dataset(dataset)
    ad, ctc, dc, ref, types = (D["adata"], D["celltype_col"], D["donor_col"], D["ref"], D["ref_types"])
    mp = {t: D["mapping"].get(t, t) for t in types}
    train = ad[ad.obs[dc].astype(str).isin(set(D["ref_donors"])).to_numpy()].copy()
    fine_types = [t for t in pd.unique(train.obs[ctc].astype(str))
                  if (train.obs[ctc].astype(str) == t).sum() >= 10]
    # per-type DE ONCE (large), cache ranked genes
    with warnings.catch_warnings():
        warnings.simplefilter("ignore")
        per_type = {t: GS.donor_aware_de(train, ctc, dc, t, mode="one_vs_rest",
                                         top_n=max(TOP_N), min_cells=10).genes for t in fine_types}
    query_mask = ad.obs[dc].astype(str).isin(set(D["query_donors"])).to_numpy()
    queries = []
    for seed in seeds:
        counts, truth, _ = H.generate_pseudobulk(ad[query_mask].copy(), ctc,
                                                 n_per_regime=2, n_cells=200, seed=seed)
        queries.append((seed, counts, truth))

    def _solve_score(genes, method, budget_kind, n_req):
        for seed, counts, truth in queries:
            with warnings.catch_warnings():
                warnings.simplefilter("ignore")
                props = PoissonGLMSolver(genes=genes).solve(counts, ref).proportions
            sc = _score(truth, props, mp)
            sc.update(dataset=dataset, seed=seed, method=method, budget_kind=budget_kind,
                      requested=n_req, n_genes=len(genes), n_fine_types=len(fine_types))
            rows.append(sc)

    # fine_global at each top_n (natural, scales with #types)
    for tn in TOP_N:
        genes = list(dict.fromkeys(g for t in fine_types for g in per_type[t][:tn]))
        _solve_score(genes, "fine_global", "per_type", tn)
    # current GeneSelector at total budgets
    for B in TOTAL_BUDGETS:
        with warnings.catch_warnings():
            warnings.simplefilter("ignore")
            cur = GeneSelector(config=GeneConfig(n_genes=B)).select(ref).selected_genes
        _solve_score(cur, "current", "total", B)


def main(argv=None):
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument("--run-real-data", action="store_true")
    ap.add_argument("--datasets", nargs="+", default=["breast", "lung"])
    ap.add_argument("--seeds", type=int, nargs="+", default=[0, 1, 2])
    args = ap.parse_args(argv)
    if not H.real_data_enabled(args.run_real_data):
        print(f"Refusing without --run-real-data (or {H.REAL_DATA_ENV}=1).", file=sys.stderr)
        return 2
    OUT.mkdir(parents=True, exist_ok=True)
    rows = []
    for ds in args.datasets:
        run(ds, args.seeds, rows)
    df = pd.DataFrame(rows)
    df.to_csv(OUT / "budget_curve.tsv", sep="\t", index=False)
    by = df.groupby(["dataset", "method", "n_genes", "n_fine_types"])[["RMSE", "Pearson", "cond_RMSE"]].mean().round(4).reset_index()
    by.to_csv(OUT / "budget_curve_by_tissue.tsv", sep="\t", index=False)
    # Pareto: for each tissue, non-dominated (min n_genes, min RMSE)
    pareto = []
    for ds in df.dataset.unique():
        d = by[by.dataset == ds].sort_values("n_genes")
        best = np.inf
        for _, r in d.iterrows():
            if r.RMSE < best:            # fewer genes AND lower RMSE than all smaller panels
                best = r.RMSE; pareto.append(r.to_dict())
    pd.DataFrame(pareto).to_csv(OUT / "pareto_front.tsv", sep="\t", index=False)
    # best fine_global budget per tissue (min RMSE), to test the ~550 hypothesis
    print("\n=== BUDGET CURVE (fine_global per_type sweep + current) ===")
    print(by.to_string(index=False))
    for ds in by.dataset.unique():
        fg = by[(by.dataset == ds) & (by.method == "fine_global")]
        if not fg.empty:
            b = fg.loc[fg.RMSE.idxmin()]
            print(f"[{ds}] best fine_global RMSE={b.RMSE} at n_genes={int(b.n_genes)} "
                  f"(n_fine_types={int(b.n_fine_types)})")
    (OUT / "budget_curve_manifest.json").write_text(json.dumps(
        {"datasets": args.datasets, "seeds": args.seeds, "top_n": TOP_N,
         "total_budgets": TOTAL_BUDGETS, "solver": "poisson"}, indent=2))
    print(f"\nWrote {OUT}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
