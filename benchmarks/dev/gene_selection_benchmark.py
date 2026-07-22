#!/usr/bin/env python
"""Gene-selection benchmark (Etapa 3) — donor-held-out, §28-gated (dev, real-data).

Tests whether donor-aware gene-selection strategies beat the current markers on
donor-held-out **deconvolution** (not classification), using the Poisson solver on a
shared reference. Answers Q11 across the promotion gates: breast + lung, 5 donor-disjoint
seeds, cross-platform (breast chemistries), with rare recall/FPR, spillover, conditional
within-family RMSE, and signature condition number.

Gene sets (built on TRAIN donors only, per scenario):
  current_markers   GeneSelector (baseline)
  donor_de_ovr      donor-aware pseudobulk DE, one-vs-rest union
  ml_minimal        L1-logistic donor-held-out minimal set (variance-prefiltered)
  hybrid_DE_ML      DE candidates -> ML-minimal
  all_shared        every shared gene (no selection)

Experimental; outputs gitignored under benchmarks/results/gene_selection/; changes NO
defaults. Usage:
  PYTHONPATH=src:. python benchmarks/dev/gene_selection_benchmark.py \
      --run-real-data --datasets breast lung --seeds 0 1 2 3 4
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

OUT = REPO / "benchmarks" / "results" / "gene_selection"
MIN_CELLS = 10
PRESENCE_THR = 0.01
RARE_MAX = 0.05


def _build_ref(adata_sub, ctc):
    with warnings.catch_warnings():
        warnings.simplefilter("ignore")
        return H.prepare_reference(adata_sub.copy(), min_cells=30, cell_type_col=ctc,
                                   estimate_overdispersion=True).reference


def _variance_prefilter(adata, ctc, dc, n=3000):
    from tissueresolve.reference.gene_selection import donor_pseudobulk
    X, _, _, genes = donor_pseudobulk(adata, ctc, dc, min_cells=MIN_CELLS)
    var = X.var(axis=0)
    return [genes[i] for i in np.argsort(-var)[:n]]


def _gene_sets(train_adata, ref, ctc, dc, types, seed):
    from tissueresolve.reference.markers import GeneSelector
    from tissueresolve.reference import gene_selection as GS
    sets = {}
    with warnings.catch_warnings():
        warnings.simplefilter("ignore")
        sets["current_markers"] = GeneSelector().select(ref).selected_genes
        sets["donor_de_ovr"] = list(dict.fromkeys(
            g for t in types
            for g in GS.donor_aware_de(train_adata, ctc, dc, t, top_n=15,
                                       min_cells=MIN_CELLS).genes))
        cand = _variance_prefilter(train_adata, ctc, dc, n=3000)
        sets["ml_minimal"] = GS.ml_minimal_genes(
            train_adata, ctc, dc, candidate_genes=cand,
            sizes=(20, 50, 100, 150, 200), seed=seed).genes
        sets["hybrid_DE_ML"] = GS.hybrid_gene_set(
            train_adata, ctc, dc, types, de_top_n=40,
            ml_sizes=(50, 100, 150), min_cells=MIN_CELLS, seed=seed).genes
    sets["all_shared"] = None
    return sets


def _score(truth, props, mapping):
    from tissueresolve.solver import SolverResult  # noqa: F401 (keep import cheap)
    p = props.reindex(index=truth.index, columns=truth.columns).fillna(0.0)
    t, pv = truth.to_numpy(float), p.to_numpy(float)
    rmse = float(np.sqrt(np.mean((pv - t) ** 2)))
    pear = float(np.corrcoef(pv.ravel(), t.ravel())[0, 1])
    # conditional within-family RMSE
    fam = {}
    for c in truth.columns:
        fam.setdefault(str(mapping.get(c, c)), []).append(c)
    multi = [c for m in fam.values() if len([x for x in m if x in truth.columns]) > 1
             for c in m if c in truth.columns]

    def cond(df):
        o = pd.DataFrame(0.0, index=df.index, columns=truth.columns)
        for m in fam.values():
            m = [x for x in m if x in df.columns]
            if len(m) < 2:
                continue
            s = df[m].sum(1).replace(0, np.nan)
            for x in m:
                o[x] = (df[x] / s).fillna(0.0)
        return o
    cr = (float(np.sqrt(np.mean((cond(truth)[multi].to_numpy(float)
                                 - cond(p)[multi].to_numpy(float)) ** 2)))
          if multi else float("nan"))
    absent = t < 1e-9
    present = ~absent
    rare_present = present & (t < RARE_MAX)
    pos = pv > PRESENCE_THR
    fp = int((pos & absent).sum()); tn = int((~pos & absent).sum())
    return {
        "RMSE": rmse, "Pearson": pear, "cond_RMSE": cr,
        "rare_recall": (int((pos & rare_present).sum()) / int(rare_present.sum()))
        if rare_present.sum() else np.nan,
        "rare_fpr": fp / (fp + tn) if (fp + tn) else np.nan,
        "spillover_absent_mass": float(np.where(absent, pv, 0.0).sum(1).mean()),
    }


def _cond_number(ref, genes):
    from tissueresolve.solver.base import align_query_to_reference  # noqa: F401
    g = [x for x in (genes or list(ref.gene_names)) if x in set(map(str, ref.gene_names))]
    sub = ref.subset_genes(g)
    R = sub.as_R_cpm().T
    try:
        return float(np.linalg.cond(R))
    except Exception:
        return float("nan")


def _run(dataset, seeds, scenarios):
    from tissueresolve.solver import PoissonGLMSolver
    D = _load_dataset(dataset)
    ad, ctc, dc, mp, types = (D["adata"], D["celltype_col"], D["donor_col"],
                              D["mapping"], D["ref_types"])
    rows = []
    for scen in scenarios:
        if scen == "complete":
            ref = D["ref"]
            train = ad[ad.obs[dc].astype(str).isin(set(D["ref_donors"])).to_numpy()].copy()
            qmask = ad.obs[dc].astype(str).isin(set(D["query_donors"])).to_numpy()
        elif scen == "cross_platform":
            if "assay" not in ad.obs.columns or ad.obs["assay"].nunique() < 2:
                print(f"[{dataset}/{scen}] N/A"); continue
            assay = ad.obs["assay"].astype(str)
            pa, pb = assay.value_counts().index.tolist()[:2]
            train = ad[(assay == pa).to_numpy()].copy()
            ref = _build_ref(train, ctc)
            qmask = (assay == pb).to_numpy()
        else:
            continue
        t0 = time.time()
        sets = _gene_sets(train, ref, ctc, dc, [t for t in types if t in
                          set(train.obs[ctc].astype(str))], seed=seeds[0])
        print(f"[{dataset}/{scen}] gene sets built in {time.time()-t0:.0f}s: "
              + ", ".join(f"{k}={len(v) if v else 'all'}" for k, v in sets.items()))
        cond = {k: _cond_number(ref, v) for k, v in sets.items()}
        for seed in seeds:
            counts, truth, _ = H.generate_pseudobulk(ad[qmask].copy(), ctc,
                                                     n_per_regime=2, n_cells=200, seed=seed)
            for name, gset in sets.items():
                with warnings.catch_warnings():
                    warnings.simplefilter("ignore")
                    props = PoissonGLMSolver(genes=gset).solve(counts, ref).proportions
                sc = _score(truth, props, mp)
                sc.update(dataset=dataset, scenario=scen, seed=seed, gene_set=name,
                          n_genes=(len(gset) if gset else int(props.shape[1] and
                                   len(ref.gene_names))), condition_number=cond[name])
                rows.append(sc)
    return rows


def main(argv=None):
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument("--run-real-data", action="store_true")
    ap.add_argument("--datasets", nargs="+", default=["breast", "lung"])
    ap.add_argument("--seeds", type=int, nargs="+", default=[0, 1, 2, 3, 4])
    ap.add_argument("--scenarios", nargs="+", default=["complete", "cross_platform"])
    args = ap.parse_args(argv)
    if not H.real_data_enabled(args.run_real_data):
        print(f"Refusing without --run-real-data (or {H.REAL_DATA_ENV}=1).", file=sys.stderr)
        return 2
    OUT.mkdir(parents=True, exist_ok=True)
    rows = []
    for ds in args.datasets:
        rows += _run(ds, args.seeds, args.scenarios)
    df = pd.DataFrame(rows)
    df.to_csv(OUT / "gene_selection_results.tsv", sep="\t", index=False)
    summ = (df.groupby(["dataset", "scenario", "gene_set"])
            [["RMSE", "Pearson", "cond_RMSE", "rare_recall", "rare_fpr",
              "spillover_absent_mass", "n_genes", "condition_number"]]
            .mean(numeric_only=True).round(4))
    summ.to_csv(OUT / "gene_selection_summary.tsv", sep="\t")
    (OUT / "manifest.json").write_text(json.dumps(
        {"datasets": args.datasets, "seeds": args.seeds, "scenarios": args.scenarios,
         "solver": "poisson", "presence_thr": PRESENCE_THR}, indent=2))
    print("\n=== GENE-SELECTION SUMMARY (seed means, Poisson, donor-held-out) ===")
    print(summ.to_string())
    print(f"\nWrote {OUT}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
