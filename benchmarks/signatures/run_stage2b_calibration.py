#!/usr/bin/env python
"""Stage 2B — calibration of the identifiability certificate (breast + lung).

Validates the calibrated (nonnegative, swap-based, donor-uncertainty-aware) certificate against
donor-disjoint held-out deconvolution. Donor-disjoint TRAIN / CALIBRATION / TEST split (reference =
train; class->probability map fit on calibration; every metric measured on test). Deconvolver =
production Poisson GLM, the same solver the certificate simulates.

Reports, across a realistic depth/panel sweep, per dataset:
  * PREDICTED vs EMPIRICAL within-cluster swap RMSE (does the certificate predict confounding?);
  * swap-based FALSE-MERGE and FALSE-RESOLUTION rates + Brier/ECE (vs the abundance-dominated 2A
    criterion, which gave false-merge ~1.0);
  * donor-aware PREDICTED-interval coverage vs the 2A-style gene-bootstrap coverage (which
    under-covered 3-4x) at nominal 90%.

No default changed; no new solver (the certificate reuses the production solver in a forward
simulation for analysis only). Third tissue: NOT EXECUTED (no donor-annotated cell-typed reference
offline; see run_stage2b_third_tissue.py).

Usage:
  PYTHONPATH=src:. python benchmarks/signatures/run_stage2b_calibration.py --run-real-data \
      --datasets breast lung --seeds 0 1 2 3 4
"""
from __future__ import annotations

import argparse
import json
import platform
import sys
import time
import warnings
from collections import defaultdict
from pathlib import Path

import numpy as np
import pandas as pd

REPO = Path(__file__).resolve().parents[2]
sys.path.insert(0, str(REPO / "examples" / "real_breast_cancer" / "scripts"))
sys.path.insert(0, str(REPO))
sys.path.insert(0, str(REPO / "src"))

from tissueresolve.reference.identifiability_calibration import (  # noqa: E402
    calibrated_identifiability_certificate as ccert)

OUT = REPO / "benchmarks" / "results" / "signatures" / "stage2b"

SWAP_CONFOUNDED = 0.15      # empirical within-cluster swap RMSE >= this => not recoverable (pre-specified)
RMSE_RESOLVABLE = 0.10      # standalone-type held-out RMSE threshold
NOMINAL = 0.90
Z = 1.645
BOOT_B, BOOT_FRAC, BOOT_MAX = 20, 0.8, 24
N_SIM = 200
CLASS_PROB_PRIOR = {"RESOLVABLE": 0.90, "WEAKLY_RESOLVABLE": 0.60, "GROUP_ONLY": 0.20, "UNRESOLVABLE": 0.10}
# depth factor of the (~1e6-count) pseudobulk; panel size (None = all genes)
CONDITIONS = [("full_depth", 1.0, None), ("shallow_0p1", 0.1, None), ("panel_1000", 1.0, 1000)]


def _thin(counts, factor, rng):
    if factor >= 1.0:
        return counts
    return pd.DataFrame(rng.binomial(counts.to_numpy().astype(np.int64), factor),
                        index=counts.index, columns=counts.columns)


def _pool(query, ctc, seeds, n_per_regime, H):
    cl, tl = [], []
    for s in seeds:
        c, t, _ = H.generate_pseudobulk(query, ctc, n_per_regime=n_per_regime, n_cells=200, seed=s)
        cl.append(c); tl.append(t)
    counts = pd.concat(cl, axis=1)
    counts.columns = [f"m{i}" for i in range(counts.shape[1])]
    truth = pd.concat(tl, axis=0).reset_index(drop=True)
    truth.index = counts.columns
    return counts, truth


def _empirical_swap(truth, pred, groups, types):
    """Empirical within-cluster swap RMSE for each certificate cluster (member/cluster-sum)."""
    out = {}
    for g in groups:
        members = [m for m in g if m in types]
        if len(members) < 2:
            continue
        T = truth.reindex(columns=members).fillna(0).to_numpy(float)
        P = pred.reindex(columns=members).fillna(0).to_numpy(float)
        ts, ps = T.sum(1), P.sum(1)
        present = ts > 0.05
        if present.sum() < 3:
            out[";".join(members)] = (np.nan, int(present.sum()))
            continue
        tc = T[present] / ts[present, None]
        pc = P[present] / np.clip(ps[present, None], 1e-9, None)
        out[";".join(members)] = (float(np.sqrt(np.mean((tc - pc) ** 2))), int(present.sum()))
    return out


def _gene_bootstrap_coverage(counts, ref, truth, class_by_type, types, panel, rng, Solver):
    avail = [g for g in (panel or list(ref.gene_names)) if g in counts.index]
    if len(avail) < 10:
        return {}, np.nan
    cols = list(counts.columns)[:BOOT_MAX]
    c_sub = counts[cols]; Tv = truth.reindex(index=cols, columns=types).fillna(0).to_numpy(float)
    m = max(5, int(len(avail) * BOOT_FRAC))
    preds = []
    for _ in range(BOOT_B):
        gsel = list(rng.choice(avail, size=m, replace=False))
        with warnings.catch_warnings():
            warnings.simplefilter("ignore")
            p = Solver(genes=gsel).solve(c_sub, ref).proportions.reindex(index=cols, columns=types).fillna(0.0)
        preds.append(p.to_numpy(float))
    stack = np.stack(preds, 0)
    lo = np.percentile(stack, (1 - NOMINAL) / 2 * 100, 0); hi = np.percentile(stack, (1 + NOMINAL) / 2 * 100, 0)
    by = defaultdict(list)
    for j, k in enumerate(types):
        inside = (Tv[:, j] >= lo[:, j] - 1e-9) & (Tv[:, j] <= hi[:, j] + 1e-9)
        by[class_by_type.get(k, "?")].extend(inside.tolist())
    return ({c: round(float(np.mean(v)), 4) for c, v in by.items() if v},
            round(float(np.mean([x for v in by.values() for x in v])), 4))


def _predicted_coverage(truth, pred, pred_sd, class_by_type, types):
    by = defaultdict(list)
    for k in types:
        tv = truth[k].to_numpy(float) if k in truth else np.zeros(len(truth))
        pv = pred[k].to_numpy(float) if k in pred else np.zeros(len(truth))
        sd = float(pred_sd.get(k, np.nan))
        if not (sd == sd):
            continue
        inside = np.abs(pv - tv) <= Z * sd + 1e-9
        by[class_by_type.get(k, "?")].extend(inside.tolist())
    return ({c: round(float(np.mean(v)), 4) for c, v in by.items() if v},
            round(float(np.mean([x for v in by.values() for x in v])), 4) if by else np.nan)


def run(dataset, seeds, outdir):
    import _harness as H
    from benchmarks.spatial.run_weak_smoothing_grid import _load_dataset
    from benchmarks.signatures import validation as V
    from tissueresolve.solver import PoissonGLMSolver
    from tissueresolve.reference.build import ReferenceBuilder
    from tissueresolve.config import ReferenceConfig

    D = _load_dataset(dataset)
    ad, ctc, dc = D["adata"], D["celltype_col"], D["donor_col"]
    train = set(map(str, D["ref_donors"]))
    # rebuild the reference on TRAIN donors WITH donor_col so donor_cv (cross-donor CV) is populated
    # — required for the donor-level uncertainty of the calibrated certificate (focus 3).
    train_ad = ad[ad.obs[dc].astype(str).isin(train).to_numpy()].copy()
    with warnings.catch_warnings():
        warnings.simplefilter("ignore")
        ref = ReferenceBuilder(ReferenceConfig(celltype_col=ctc, donor_col=dc, min_cells=25)).build_from_adata(train_ad)
    assert ref.donor_cv is not None, "donor_cv not populated — donor-level uncertainty unavailable"
    types = [str(c) for c in ref.cell_types]
    q = sorted(map(str, D["query_donors"]))
    perm = np.random.default_rng(2024).permutation(q); half = max(1, len(perm) // 2)
    cal_donors, test_donors = sorted(perm[:half]), sorted(perm[half:])
    assert train.isdisjoint(cal_donors) and train.isdisjoint(test_donors) and set(cal_donors).isdisjoint(test_donors)
    V.validate_donor_disjoint(train, set(cal_donors) | set(test_donors))

    genes = list(ref.gene_names)
    q_cal = ad[ad.obs[dc].astype(str).isin(set(cal_donors)).to_numpy()].copy()
    q_test = ad[ad.obs[dc].astype(str).isin(set(test_donors)).to_numpy()].copy()
    n_don = {t_: int(q_test.obs.loc[q_test.obs[ctc].astype(str) == t_, dc].nunique()) for t_ in types}
    counts_test_full, truth_test = _pool(q_test, ctc, seeds, 3, H)
    counts_cal_full, truth_cal = _pool(q_cal, ctc, seeds[:3], 2, H)
    truth_test = truth_test.reindex(columns=types).fillna(0.0)
    truth_cal = truth_cal.reindex(columns=types).fillna(0.0)
    panel_1000 = list(np.random.default_rng(12345).choice(genes, size=min(1000, len(genes)), replace=False))

    swap_rows, cov_rows, type_rows, cal_rows = [], [], [], []
    for cname, depth, n_panel in CONDITIONS:
        rng = np.random.default_rng(abs(hash((dataset, cname))) % (2**32))
        panel = panel_1000 if n_panel else genes
        counts_t = _thin(counts_test_full, depth, rng)
        det_frac = (counts_t.loc[panel] > 0).mean(axis=1)
        detectable = [g for g in panel if det_frac.get(g, 0) >= 0.10]
        if len(detectable) < 2:
            continue
        lib = float(counts_t.loc[detectable].sum(axis=0).median())
        # ---- solve CALIBRATION held-out once; fit donor-shift inflation on it (donor-disjoint) ----
        counts_c = _thin(counts_cal_full, depth, np.random.default_rng(7))
        with warnings.catch_warnings():
            warnings.simplefilter("ignore")
            pred_c = PoissonGLMSolver(genes=(panel if n_panel else None)).solve(
                counts_c, ref).proportions.reindex(columns=types).fillna(0.0)
        best_s, best_gap, sd_s1 = 1.0, np.inf, None
        for s in (1.0, 2.0, 3.0, 5.0, 8.0):
            cs = ccert(ref, query_detectable_genes=detectable, library_size=lib, n_sim=100,
                       n_donors_by_type=n_don, shift_scale=s, seed=0)
            sd_s = dict(zip(cs.per_type.cell_type, cs.per_type.pred_sd))
            cby = dict(zip(cs.per_type.cell_type, cs.per_type.recoverability))
            if s == 1.0:
                sd_s1, cby_s1 = sd_s, cby
            _, cov_cal_s = _predicted_coverage(truth_cal, pred_c, sd_s, cby, types)
            if cov_cal_s == cov_cal_s and abs(cov_cal_s - NOMINAL) < best_gap:
                best_gap, best_s = abs(cov_cal_s - NOMINAL), s
        cert = ccert(ref, query_detectable_genes=detectable, library_size=lib, n_sim=N_SIM,
                     n_donors_by_type=n_don, shift_scale=best_s, seed=0)
        if cname == "panel_1000":
            cert.save(outdir)
        class_by = dict(zip(cert.per_type.cell_type, cert.per_type.recoverability))
        pred_sd = dict(zip(cert.per_type.cell_type, cert.per_type.pred_sd))
        groups = [g["group"].split(";") for _, g in cert.clusters.iterrows()] if len(cert.clusters) else []
        pred_swap = {g["group"]: g["swap_rmse"] for _, g in cert.clusters.iterrows()} if len(cert.clusters) else {}

        with warnings.catch_warnings():
            warnings.simplefilter("ignore")
            pred_t = PoissonGLMSolver(genes=(panel if n_panel else None)).solve(
                counts_t, ref).proportions.reindex(columns=types).fillna(0.0)
        V.validate_truth_pred_alignment(truth_test, pred_t)

        emp_swap = _empirical_swap(truth_test, pred_t, groups, types)
        for gkey, (esw, npres) in emp_swap.items():
            swap_rows.append(dict(dataset=dataset, condition=cname, group=gkey,
                                  predicted_swap=pred_swap.get(gkey, np.nan),
                                  empirical_swap=esw, n_present=npres,
                                  n_members=len(gkey.split(";"))))
        # per-type empirical held-out RMSE + swap-based recoverable outcome
        cluster_of = {t: gkey for gkey in emp_swap for t in gkey.split(";")}
        for k in types:
            tv = truth_test[k].to_numpy(float); pv = pred_t[k].to_numpy(float)
            hr = float(np.sqrt(np.mean((pv - tv) ** 2)))
            gk = cluster_of.get(k)
            esw = emp_swap.get(gk, (np.nan, 0))[0] if gk else np.nan
            rec_emp = (esw < SWAP_CONFOUNDED) if (esw == esw) else (hr < RMSE_RESOLVABLE)
            type_rows.append(dict(dataset=dataset, condition=cname, cell_type=k,
                                  recoverability=class_by[k], held_out_rmse=round(hr, 4),
                                  empirical_swap=round(esw, 4) if esw == esw else np.nan,
                                  recoverable_empirical=bool(rec_emp)))
        # coverage (predicted donor-aware vs gene-bootstrap baseline)
        pcov, pov = _predicted_coverage(truth_test, pred_t, pred_sd, class_by, types)
        gcov, gov = _gene_bootstrap_coverage(counts_t, ref, truth_test, class_by, types,
                                             (panel if n_panel else None), np.random.default_rng(9), PoissonGLMSolver)
        _, pov_s1 = _predicted_coverage(truth_test, pred_t, sd_s1, cby_s1, types)
        cov_rows.append(dict(dataset=dataset, condition=cname, method="predicted_uncalibrated_s1",
                             overall=pov_s1, shift_scale=1.0, by_class={}))
        cov_rows.append(dict(dataset=dataset, condition=cname, method="predicted_donor_aware_calibrated",
                             overall=pov, shift_scale=best_s, by_class=pcov))
        cov_rows.append(dict(dataset=dataset, condition=cname, method="gene_bootstrap_baseline",
                             overall=gov, shift_scale=np.nan, by_class=gcov))
        # calibration donors (pred_c solved above): swap-based recoverable outcome for class_prob fit
        emp_swap_c = _empirical_swap(truth_cal, pred_c, groups, types)
        cluster_of_c = {t: gk for gk in emp_swap_c for t in gk.split(";")}
        for k in types:
            hr = float(np.sqrt(np.mean((pred_c[k].to_numpy(float) - truth_cal[k].to_numpy(float)) ** 2)))
            gk = cluster_of_c.get(k); esw = emp_swap_c.get(gk, (np.nan, 0))[0] if gk else np.nan
            rec = (esw < SWAP_CONFOUNDED) if esw == esw else (hr < RMSE_RESOLVABLE)
            cal_rows.append(dict(cell_type=k, recoverability=class_by[k], recoverable_empirical=bool(rec)))

    type_df = pd.DataFrame(type_rows); swap_df = pd.DataFrame(swap_rows)
    cal_df = pd.DataFrame(cal_rows); cov_df = pd.DataFrame(cov_rows)
    type_df.to_csv(outdir / "recoverability_vs_empirical.tsv", sep="\t", index=False)
    swap_df.to_csv(outdir / "predicted_vs_empirical_swap.tsv", sep="\t", index=False)
    cov_df.to_csv(outdir / "uncertainty_coverage.tsv", sep="\t", index=False)

    # class->prob fit on CAL (swap-based), evaluate Brier/ECE + false rates on TEST
    class_prob = dict(CLASS_PROB_PRIOR)
    for cls, grp in cal_df.groupby("recoverability"):
        if len(grp) >= 5:
            class_prob[cls] = float(np.clip(grp.recoverable_empirical.mean(), 0.02, 0.98))
    d = type_df[type_df.recoverability.isin(class_prob)].copy()
    d["p"] = d.recoverability.map(class_prob); d["y"] = d.recoverable_empirical.astype(float)
    brier = float(np.mean((d.p - d.y) ** 2)) if len(d) else np.nan
    ece = float(sum((len(g) / len(d)) * abs(p - g.y.mean()) for p, g in d.groupby("p"))) if len(d) else np.nan
    resolved = d[d.recoverability == "RESOLVABLE"]; merged = d[d.recoverability.isin(["GROUP_ONLY", "UNRESOLVABLE"])]
    false_res = float((~resolved.y.astype(bool)).mean()) if len(resolved) else np.nan
    false_merge = float(merged.y.mean()) if len(merged) else np.nan

    # predicted vs empirical swap correlation
    sd = swap_df.dropna(subset=["predicted_swap", "empirical_swap"])
    swap_corr = float(sd.predicted_swap.corr(sd.empirical_swap)) if len(sd) >= 3 else np.nan
    swap_rank = float(sd.predicted_swap.corr(sd.empirical_swap, method="spearman")) if len(sd) >= 3 else np.nan

    summary = dict(dataset=dataset, swap_pearson=round(swap_corr, 4) if swap_corr == swap_corr else np.nan,
                   swap_spearman=round(swap_rank, 4) if swap_rank == swap_rank else np.nan,
                   n_clusters=int(len(sd)), brier_swapbased=round(brier, 4) if brier == brier else np.nan,
                   ece_swapbased=round(ece, 4) if ece == ece else np.nan,
                   false_resolution_rate=round(false_res, 4) if false_res == false_res else np.nan,
                   false_merge_rate=round(false_merge, 4) if false_merge == false_merge else np.nan,
                   class_prob={k: round(v, 3) for k, v in class_prob.items()})
    pd.DataFrame([summary]).to_csv(outdir / "calibration_summary.tsv", sep="\t", index=False)

    print(f"\n=== {dataset}: PREDICTED vs EMPIRICAL swap (per cluster) ===")
    if len(sd):
        print(sd[["condition", "group", "predicted_swap", "empirical_swap", "n_members", "n_present"]]
              .to_string(index=False, max_colwidth=44))
    print(f"swap Pearson={summary['swap_pearson']} Spearman={summary['swap_spearman']} (n={len(sd)} clusters)")
    print(f"swap-based: false_merge={summary['false_merge_rate']} false_resolution={summary['false_resolution_rate']} "
          f"Brier={summary['brier_swapbased']} ECE={summary['ece_swapbased']}")
    print("=== uncertainty coverage @90% (uncalibrated s=1 vs CAL-calibrated donor-aware vs gene-bootstrap) ===")
    print(cov_df[["condition", "method", "overall", "shift_scale"]].to_string(index=False))
    return summary


def main(argv=None):
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument("--run-real-data", action="store_true")
    ap.add_argument("--datasets", nargs="+", default=["breast", "lung"])
    ap.add_argument("--seeds", type=int, nargs="+", default=[0, 1, 2, 3, 4])
    args = ap.parse_args(argv)
    import _harness as H
    if not H.real_data_enabled(args.run_real_data):
        print(f"Refusing without --run-real-data (or {H.REAL_DATA_ENV}=1).", file=sys.stderr)
        return 2
    OUT.mkdir(parents=True, exist_ok=True)
    versions = {m: __import__(m).__version__ for m in ("numpy", "pandas", "scipy", "anndata")}
    versions["python"] = platform.python_version()
    t0 = time.perf_counter()
    summ = []
    for ds in args.datasets:
        od = OUT / ds; od.mkdir(parents=True, exist_ok=True)
        summ.append(run(ds, args.seeds, od))
    pd.DataFrame([{"package": k, "version": v} for k, v in versions.items()]).to_csv(
        OUT / "software_versions.tsv", sep="\t", index=False)
    pd.DataFrame(summ).to_csv(OUT / "stage2b_summary.tsv", sep="\t", index=False)
    (OUT / "benchmark_manifest.json").write_text(json.dumps(
        {"datasets": args.datasets, "seeds": args.seeds, "conditions": [c[0] for c in CONDITIONS],
         "n_sim": N_SIM, "swap_confounded": SWAP_CONFOUNDED, "nominal_coverage": NOMINAL,
         "third_tissue": "NOT EXECUTED (no donor-annotated cell-typed reference offline)",
         "runtime_s": round(time.perf_counter() - t0, 1)}, indent=2, default=str))
    print(f"\nWrote {OUT}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
