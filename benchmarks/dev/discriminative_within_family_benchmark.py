#!/usr/bin/env python
"""Test the discriminative-marker-reweighted within-family split on COLLINEAR T/NK
states (breast pseudobulk gold-truth), vs the plain (uniform-weight) mean NNLS split.

Setup mirrors the pilot Probe B: T/NK-restricted pseudobulk from held-out query
donors; reference = mean profile per T/NK state from ref donors (top HVG). The
within-family split is solved by NNLS with (a) uniform gene weights [baseline] vs
(b) discriminative gene weights ∝ between-sibling variance, with a power sweep.

Question: does discriminative reweighting lower the conditional within-T/NK RMSE
(and improve the collinear-pair recovery) vs the mean baseline?

Experimental; outputs gitignored; no defaults changed.

Usage:
  PYTHONPATH=src:. python benchmarks/dev/discriminative_within_family_benchmark.py \
      --run-real-data --seeds 0 1 2 3 4
"""
from __future__ import annotations

import argparse
import json
import sys
from pathlib import Path

import numpy as np
import pandas as pd

REPO = Path(__file__).resolve().parents[2]
sys.path.insert(0, str(REPO / "examples" / "real_breast_cancer" / "scripts"))
sys.path.insert(0, str(REPO))
sys.path.insert(0, str(REPO / "src"))
import _harness as H  # noqa: E402
from benchmarks.spatial.run_weak_smoothing_grid import _load_dataset  # noqa: E402

OUT = REPO / "benchmarks" / "outputs" / "discriminative_within_family"
N_HVG = 800
MIN_CELLS_PER_CLASS = 60
N_SAMPLES = 16
CELLS_PER_SAMPLE = 400
POWERS = [0.0, 1.0, 2.0, 4.0]   # 0.0 = uniform baseline


def _dense(adata):
    import scipy.sparse as sp
    X = adata.X
    return np.asarray(X.todense()) if sp.issparse(X) else np.asarray(X)


def _lognorm(X):
    lib = np.maximum(X.sum(1, keepdims=True), 1.0)
    return np.log1p(X / lib * 1e4)


def main(argv=None) -> int:
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument("--run-real-data", action="store_true")
    ap.add_argument("--seeds", type=int, nargs="+", default=[0, 1, 2, 3, 4])
    args = ap.parse_args(argv)
    if not H.real_data_enabled(args.run_real_data):
        print(f"Refusing without --run-real-data (or {H.REAL_DATA_ENV}=1).", file=sys.stderr)
        return 2

    from tissueresolve.experimental.discriminative_within_family import (
        discriminative_gene_weights, solve_conditional_split)

    D = _load_dataset("breast")
    adata, mapping, ctc, dc = D["adata"], D["mapping"], D["celltype_col"], D["donor_col"]
    tnk = [t for t, f in mapping.items() if f == "T/NK"]
    obs = adata.obs
    ct = obs[ctc].astype(str).to_numpy(); dn = obs[dc].astype(str).to_numpy()
    ref_d, qry_d = set(D["ref_donors"]), set(D["query_donors"])
    refmask, qmask = np.isin(dn, list(ref_d)), np.isin(dn, list(qry_d))

    # HVG within T/NK
    tnk_mask = np.isin(ct, tnk)
    hvg = np.argsort(_lognorm(_dense(adata[tnk_mask])).var(0))[::-1][:N_HVG]

    # usable T/NK states + mean profiles (ref donors) on HVG
    states, Phi = [], []
    for s in tnk:
        ridx = np.where(refmask & (ct == s))[0]
        if len(ridx) < MIN_CELLS_PER_CLASS:
            continue
        X = _dense(adata[ridx])[:, hvg].astype(float)
        prof = X.sum(0); prof = prof / (prof.sum() + 1e-9)
        states.append(s); Phi.append(prof)
    Phi = np.array(Phi)                  # (S, G) proportion scale
    S = len(states)
    print(f"[breast T/NK] {S} states on {N_HVG} HVG; states={states}")

    # discriminative weights (per power) from the family profiles
    weights = {p: discriminative_gene_weights(Phi, power=p) for p in POWERS}

    OUT.mkdir(parents=True, exist_ok=True)
    rng = np.random.default_rng(0)
    rows = []
    for seed in args.seeds:
        rs = np.random.default_rng(seed)
        for k in range(N_SAMPLES):
            true = rs.dirichlet(np.ones(S))
            bulk = np.zeros(N_HVG)
            for si, s in enumerate(states):
                qidx = np.where(qmask & (ct == s))[0]
                if len(qidx) == 0:
                    continue
                ncell = max(1, int(round(true[si] * CELLS_PER_SAMPLE)))
                pick = rs.choice(qidx, ncell, replace=len(qidx) < ncell)
                bulk += _dense(adata[pick])[:, hvg].astype(float).sum(0)
            bn = bulk / (bulk.sum() + 1e-9)
            for p in POWERS:
                theta = solve_conditional_split(bn, Phi.T, weights[p], total_mass=1.0)
                rows.append({"seed": seed, "sample": k, "power": p,
                             "rmse": float(np.sqrt(np.mean((theta - true) ** 2))),
                             "pearson": float(np.corrcoef(theta, true)[0, 1])
                             if np.std(theta) > 0 else 0.0,
                             **{f"err[{states[i]}]": float(theta[i] - true[i]) for i in range(S)}})
    df = pd.DataFrame(rows)
    df.to_csv(OUT / "per_sample_metrics.tsv", sep="\t", index=False)
    summ = df.groupby("power").agg(rmse_mean=("rmse", "mean"),
                                   pearson_mean=("pearson", "mean")).reset_index()
    summ.to_csv(OUT / "by_power_summary.tsv", sep="\t", index=False)

    base = summ[summ.power == 0.0].iloc[0]
    best = summ.sort_values("rmse_mean").iloc[0]
    verdict = {
        "baseline_uniform_rmse": round(float(base.rmse_mean), 4),
        "baseline_uniform_pearson": round(float(base.pearson_mean), 3),
        "best_power": float(best.power),
        "best_rmse": round(float(best.rmse_mean), 4),
        "best_pearson": round(float(best.pearson_mean), 3),
        "discriminative_helps": bool(best.power > 0 and best.rmse_mean < base.rmse_mean - 1e-6),
        "rmse_reduction_pct": round(float((base.rmse_mean - best.rmse_mean) / base.rmse_mean * 100), 1),
    }
    (OUT / "verdict.json").write_text(json.dumps(verdict, indent=2), encoding="utf-8")
    (OUT / "manifest.json").write_text(json.dumps({
        "dataset": "breast", "family": "T/NK", "n_states": S, "states": states,
        "n_hvg": N_HVG, "powers": POWERS, "seeds": args.seeds,
        "note": "experimental; discriminative within-family reweighting; defaults unchanged",
    }, indent=2), encoding="utf-8")
    print("\n=== conditional within-T/NK split: discriminative-weight power sweep ===")
    print(summ.round(4).to_string(index=False))
    print("\nVERDICT:", json.dumps(verdict))
    print(f"Wrote -> {OUT}/")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
