#!/usr/bin/env python
"""Stage 2A — Experiment B: does cohort cross-sample covariance recover information that
per-sample deconvolution cannot, when the reference R is KNOWN?

Two parts:
  1. SYNTHETIC FACTOR GRID (always runs, offline, deterministic): sweep collinearity ×
     proportion-variation × depth × cohort-size, comparing single / first_moment /
     cross_covariance / oracle. Measures overall RMSE and error along R's null direction.
  2. REAL COHORT (breast + lung, behind --run-real-data): build a donor-disjoint reference,
     form a cohort of pseudobulk samples with known truth, and compare the same four estimators
     overall and specifically on the certificate's confounded pairs, across cohort sizes.

Answers §5 decision D/E. Changes no default; the estimators are experimental (not registered).

Usage:
  PYTHONPATH=src:. python benchmarks/signatures/run_stage2a_cohort.py            # synthetic only
  PYTHONPATH=src:. python benchmarks/signatures/run_stage2a_cohort.py --run-real-data
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

from tissueresolve.experimental.cohort_covariance import cohort_deconvolve, null_direction  # noqa: E402
from tissueresolve.reference.identifiability import identifiability_certificate  # noqa: E402

OUT = REPO / "benchmarks" / "results" / "signatures" / "stage2a" / "cohort"
METHODS = ("single", "first_moment", "cross_covariance", "oracle")


def _rms(x):
    return float(np.sqrt(np.mean(np.square(x))))


# ----------------------------- synthetic factor grid -----------------------------
def _make_R(kind, G=120):
    """well_conditioned: 3 orthogonal blocks; confounded: types 0,1 share a profile."""
    R = np.zeros((G, 3))
    R[:40, 0] = 1.0
    R[40:80, 1] = 1.0 if kind == "well_conditioned" else 0.0
    if kind == "confounded":
        R[:40, 1] = 1.0            # col 1 == col 0  → null direction (1,-1,0)
    R[80:120, 2] = 1.0
    return R


def _make_theta(variation, K, N, rng):
    if variation == "constant":
        return np.tile(np.array([0.5, 0.3, 0.2])[:, None], (1, N))
    if variation == "clustered":
        return np.array([0.5, 0.3, 0.2])[:, None] + rng.normal(0, 0.03, (K, N))
    return np.abs(rng.normal(0.4, 0.15, (K, N)))          # varying (independent)


def synthetic_grid():
    rows = []
    for kind in ("well_conditioned", "confounded"):
        for variation in ("constant", "clustered", "varying"):
            for sigma in (0.05, 0.5):
                for N in (10, 50, 250):
                    rng = np.random.default_rng(hash((kind, variation, int(sigma * 100), N)) % (2**32))
                    R = _make_R(kind); K = R.shape[1]
                    Th = _make_theta(variation, K, N, rng)
                    Y = R @ Th + rng.normal(0, sigma, (R.shape[0], N))
                    nd = null_direction(R)
                    est = {m: cohort_deconvolve(Y, R, m, theta_true=Th, sigma2=sigma**2).theta
                           for m in METHODS}
                    row = dict(kind=kind, variation=variation, sigma=sigma, N=N)
                    for m in METHODS:
                        row[f"rmse_{m}"] = _rms(est[m] - Th)
                        if nd is not None:
                            row[f"nullrmse_{m}"] = _rms((est[m] - Th).T @ nd)
                    row["gain_crosscov_vs_single"] = row["rmse_single"] - row["rmse_cross_covariance"]
                    if nd is not None:
                        row["nullgain_crosscov_vs_single"] = (
                            row["nullrmse_single"] - row["nullrmse_cross_covariance"])
                    rows.append(row)
    return pd.DataFrame(rows)


# ----------------------------- real cohort -----------------------------
def _confounded_pairs(cert, types):
    pairs = []
    for g in cert.confounded_groups:
        members = [m for m in g["group"].split(";") if m in types]
        for i in range(len(members)):
            for j in range(i + 1, len(members)):
                pairs.append((members[i], members[j]))
    return pairs


def real_cohort(dataset, sizes, seeds, rows, paired):
    import _harness as H
    from benchmarks.spatial.run_weak_smoothing_grid import _load_dataset
    from benchmarks.signatures import validation as V
    from tissueresolve.solver import PoissonGLMSolver

    D = _load_dataset(dataset)
    ad, ctc, dc, ref, types = (D["adata"], D["celltype_col"], D["donor_col"], D["ref"], D["ref_types"])
    V.validate_donor_disjoint(set(D["ref_donors"]), set(D["query_donors"]))
    R_cpm = ref.as_R_cpm()                                  # (K, G)
    genes = list(ref.gene_names)
    gidx = {g: i for i, g in enumerate(genes)}

    query = ad[ad.obs[dc].astype(str).isin(set(D["query_donors"])).to_numpy()].copy()
    # a large pool of pseudobulk samples (known truth) to subsample cohorts from
    pool_counts, pool_truth = [], []
    for seed in seeds:
        c, t, _ = H.generate_pseudobulk(query, ctc, n_per_regime=4, n_cells=200, seed=seed)
        pool_counts.append(c); pool_truth.append(t)
    counts = pd.concat(pool_counts, axis=1)                 # genes × samples
    truth = pd.concat(pool_truth, axis=0)                   # samples × types
    truth = truth.reindex(columns=[str(c) for c in ref.cell_types]).fillna(0.0)
    counts.columns = [f"s{i}" for i in range(counts.shape[1])]
    truth.index = counts.columns

    # Restrict to a small gene panel so structural confounding actually exists (at full depth /
    # all genes these references are well-conditioned and no confounded pairs appear — see the
    # certificate benchmark). This makes the "does cohort covariance help on confounded pairs?" test
    # non-vacuous; the overall-RMSE comparison is still valid on the panel.
    prng = np.random.default_rng(12345)
    panel = list(prng.choice(genes, size=min(1000, len(genes)), replace=False))
    det_frac = (counts.loc[panel] > 0).mean(axis=1)
    detectable = [g for g in panel if det_frac.get(g, 0) >= 0.10]
    lib = float(counts.loc[detectable].sum(axis=0).median())
    n_don = {str(t_): int(query.obs.loc[query.obs[ctc].astype(str) == str(t_), dc].nunique())
             for t_ in ref.cell_types}
    cert = identifiability_certificate(R_cpm, [str(c) for c in ref.cell_types], genes,
                                       query_detectable_genes=detectable, library_size=lib,
                                       n_donors_by_type=n_don)
    pairs = _confounded_pairs(cert, [str(c) for c in ref.cell_types])

    # align to panel genes, CPM-normalise both sides to the same scale
    ov = [g for g in panel if g in counts.index]
    Rov = R_cpm[:, [gidx[g] for g in ov]].T                # (G, K)
    Ymat = counts.loc[ov].to_numpy(float)
    Ycpm = Ymat / np.clip(Ymat.sum(0, keepdims=True), 1, None) * 1e6
    Th_true = truth.to_numpy(float).T                       # (K, Nall)
    Nall = Ycpm.shape[1]
    tcols = list(truth.columns)
    pidx = {t: i for i, t in enumerate(tcols)}

    rng = np.random.default_rng(0)
    for N in sizes:
        if N > Nall:
            continue
        for rep in range(5):
            sel = rng.choice(Nall, size=N, replace=False)
            Y = Ycpm[:, sel]; Tt = Th_true[:, sel]
            ests = {}
            for m in METHODS:
                ests[m] = cohort_deconvolve(Y, Rov, m, theta_true=Tt,
                                            nonneg=(m != "oracle")).theta
            # production single-sample Poisson GLM as an external reference point
            with warnings.catch_warnings():
                warnings.simplefilter("ignore")
                pois = PoissonGLMSolver(genes=ov).solve(
                    counts.iloc[:, sel], ref).proportions.reindex(columns=tcols).fillna(0)
            pois = pois.to_numpy(float).T
            for m, est in list(ests.items()) + [("single_poisson", pois)]:
                overall = _rms(est - Tt)
                pair_err = np.nan
                if pairs:
                    diffs = []
                    for a, b in pairs:
                        if a in pidx and b in pidx:
                            da = (est[pidx[a]] - Tt[pidx[a]]) - (est[pidx[b]] - Tt[pidx[b]])
                            diffs.append(da)
                    if diffs:
                        pair_err = _rms(np.concatenate(diffs))
                rows.append(dict(dataset=dataset, N=N, rep=rep, method=m,
                                 overall_rmse=overall, confounded_pair_rmse=pair_err,
                                 n_confounded_pairs=len(pairs)))
    # paired: cross_covariance vs single (linear, apples-to-apples)
    df = pd.DataFrame([r for r in rows if r["dataset"] == dataset])
    for N in sorted(df.N.unique()):
        sub = df[df.N == N]
        piv = sub.pivot_table(index="rep", columns="method", values="overall_rmse")
        if {"single", "cross_covariance"} <= set(piv.columns):
            d = (piv["cross_covariance"] - piv["single"]).to_numpy()
            paired.append(dict(dataset=dataset, N=int(N), comparison="cross_covariance-vs-single",
                               metric="overall_rmse", mean_diff=float(np.mean(d)),
                               win_rate=float(np.mean(d < 0)), n=int(len(d))))


def main(argv=None):
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument("--run-real-data", action="store_true")
    ap.add_argument("--datasets", nargs="+", default=["breast", "lung"])
    ap.add_argument("--sizes", type=int, nargs="+", default=[10, 25, 50, 100])
    ap.add_argument("--seeds", type=int, nargs="+", default=[0, 1, 2, 3, 4, 5])
    args = ap.parse_args(argv)
    OUT.mkdir(parents=True, exist_ok=True)
    t0 = time.perf_counter()

    grid = synthetic_grid()
    grid.to_csv(OUT / "cohort_synthetic_grid.tsv", sep="\t", index=False)
    print("=== SYNTHETIC GRID (gain>0 ⇒ cohort helps) ===")
    show = grid.groupby(["kind", "variation", "sigma"])[
        [c for c in ("gain_crosscov_vs_single", "nullgain_crosscov_vs_single") if c in grid]].mean().round(4)
    print(show.to_string())

    rows, paired = [], []
    if args.run_real_data:
        import _harness as H
        if not H.real_data_enabled(True):
            print("real-data env not enabled", file=sys.stderr); return 2
        for ds in args.datasets:
            print(f"\n[cohort] {ds} ...")
            real_cohort(ds, args.sizes, args.seeds, rows, paired)
        rdf = pd.DataFrame(rows); rdf.to_csv(OUT / "cohort_covariance_results.tsv", sep="\t", index=False)
        pdf = pd.DataFrame(paired); pdf.to_csv(OUT / "paired_comparisons.tsv", sep="\t", index=False)
        summ = rdf.groupby(["dataset", "N", "method"])[
            ["overall_rmse", "confounded_pair_rmse"]].mean().round(4)
        print("\n=== REAL COHORT SUMMARY ===\n", summ.to_string())
        print("\n=== PAIRED cross_covariance − single (negative ⇒ cohort helps) ===\n",
              pdf.to_string(index=False))
    else:
        print("\n[real cohort skipped — pass --run-real-data]")

    (OUT / "benchmark_manifest.json").write_text(json.dumps(
        {"datasets": args.datasets if args.run_real_data else [], "sizes": args.sizes,
         "seeds": args.seeds, "methods": list(METHODS) + ["single_poisson"],
         "ran_real_data": bool(args.run_real_data),
         "runtime_s": round(time.perf_counter() - t0, 1)}, indent=2, default=str))
    print(f"\nWrote {OUT}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
