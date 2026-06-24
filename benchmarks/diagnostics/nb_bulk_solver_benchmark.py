#!/usr/bin/env python
"""Benchmark the experimental Poisson / NB GLM bulk solver (P1) vs wNNLS.

Gold-truth pseudobulk (donor-held-out) on breast + HLCA/lung, multiple scenarios
(imbalanced / rare / similar_subtypes) × seeds. Compares the default weighted-NNLS
against the opt-in count-likelihood GLM solvers and a plain-NNLS external control.
Outputs RNA-derived (mRNA) proportions; truth = true mRNA proportions. CARD is
spatial and is correctly excluded.

Usage:
  PYTHONPATH=src:. python benchmarks/diagnostics/nb_bulk_solver_benchmark.py \
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
from benchmarks.shared import synthetic_holdout as SH  # noqa: E402
from benchmarks.shared import metrics as M  # noqa: E402
from benchmarks.spatial.run_weak_smoothing_grid import _load_dataset  # noqa: E402

OUT = REPO / "benchmarks" / "outputs" / "nb_bulk_solver"
SCENARIOS = ["imbalanced", "rare", "similar_subtypes"]
N_SAMPLES = 12
CELLS_PER_SAMPLE = 500
METHODS = ["wNNLS", "poisson_glm_experimental", "nb_glm_experimental", "external_nnls_control"]


def _conditional_rmse(truth, pred, mapping):
    fam_members = {}
    for c in truth.columns:
        fam_members.setdefault(str(mapping.get(c, c)), []).append(c)
    multi, sq = [], []
    for fam, members in fam_members.items():
        members = [m for m in members if m in pred.columns]
        if len(members) < 2:
            continue
        ts = truth[members].sum(axis=1).replace(0, np.nan)
        ps = pred[members].sum(axis=1).replace(0, np.nan)
        for m in members:
            tc = (truth[m] / ts).fillna(0.0).to_numpy()
            pc = (pred[m] / ps).fillna(0.0).to_numpy()
            sq.append((tc - pc) ** 2)
            multi.append(m)
    return float(np.sqrt(np.mean(np.concatenate(sq)))) if sq else float("nan")


def _rare_metrics(truth, pred, rare, thr=0.01):
    if rare not in truth.columns or rare not in pred.columns:
        return float("nan"), float("nan"), float("nan")
    t = truth[rare].to_numpy() > thr
    p = pred[rare].to_numpy() > thr
    tp = int((t & p).sum()); fp = int((~t & p).sum())
    fn = int((t & ~p).sum()); tn = int((~t & ~p).sum())
    sens = tp / (tp + fn) if (tp + fn) else float("nan")
    prec = tp / (tp + fp) if (tp + fp) else float("nan")
    fpr = fp / (fp + tn) if (fp + tn) else float("nan")
    return sens, prec, fpr


def _score(truth, pred, mapping, rare):
    cols = [c for c in truth.columns if not str(c).startswith("unresolved_")]
    t = truth[cols]
    p = pred.reindex(index=t.index, columns=cols).fillna(0.0)
    acc = M.accuracy_metrics(t, p)
    tfam, pfam = M.aggregate_to_families(t, mapping), M.aggregate_to_families(p, mapping)
    bacc = M.accuracy_metrics(tfam, pfam)
    sens, prec, fpr = _rare_metrics(t, p, rare)
    absent = t.to_numpy(float) < 1e-9
    absent_mass = float(np.where(absent, p.to_numpy(float), 0.0).sum(axis=1).mean())
    fp_subtype = float(((p.to_numpy(float) > 0.01) & absent).mean())
    effn_p = float(np.mean([M.effective_n_populations(p.iloc[i].to_numpy()) for i in range(p.shape[0])]))
    effn_t = float(np.mean([M.effective_n_populations(t.iloc[i].to_numpy()) for i in range(t.shape[0])]))
    return {
        "broad_pearson": bacc["pearson"], "broad_rmse": bacc["rmse"],
        "fine_pearson": acc["pearson"], "fine_rmse": acc["rmse"],
        "conditional_rmse": _conditional_rmse(t, p, mapping),
        "rare_sensitivity": sens, "rare_precision": prec, "rare_fpr": fpr,
        "absent_subtype_mass": absent_mass, "false_positive_subtype_rate": fp_subtype,
        "pairwise_spillover": absent_mass,
        "effective_n_pred": effn_p, "effective_n_truth": effn_t,
    }


def _run_method(method, counts, ref, genes):
    """Return predicted fine proportions (samples × types) for a method."""
    import tissueresolve as tr
    from tissueresolve.config import TissueResolveConfig
    if method == "external_nnls_control":
        from scipy.optimize import nnls
        Phi = (ref.as_R_cpm()[:, [list(ref.gene_names).index(g) for g in genes]].T / 1e6)
        B = counts.loc[genes].to_numpy(float)
        Bn = B / np.where(B.sum(0, keepdims=True) > 0, B.sum(0, keepdims=True), 1.0)
        Pn = Phi / np.where(Phi.sum(0, keepdims=True) > 0, Phi.sum(0, keepdims=True), 1.0)
        out = np.zeros((counts.shape[1], len(ref.cell_types)))
        for n in range(counts.shape[1]):
            theta, _ = nnls(Pn, Bn[:, n])
            s = theta.sum()
            out[n] = theta / s if s > 0 else np.full(len(ref.cell_types), 1.0 / len(ref.cell_types))
        return pd.DataFrame(out, index=counts.columns, columns=list(ref.cell_types))
    cfg = TissueResolveConfig()
    cfg.bulk_solver.method = method
    res = tr.deconv_bulk(counts, ref, config=cfg, resolution_mode="none")
    return res.deconv.proportions


def main(argv=None) -> int:
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument("--run-real-data", action="store_true")
    ap.add_argument("--datasets", nargs="+", default=["breast", "lung"])
    ap.add_argument("--seeds", type=int, nargs="+", default=[0, 1, 2, 3, 4])
    args = ap.parse_args(argv)
    if not H.real_data_enabled(args.run_real_data):
        print(f"Refusing without --run-real-data (or {H.REAL_DATA_ENV}=1).", file=sys.stderr)
        return 2

    OUT.mkdir(parents=True, exist_ok=True)
    rows, runtime_rows = [], []

    for ds in args.datasets:
        D = _load_dataset(ds)
        adata, ref, mapping = D["adata"], D["ref"], D["mapping"]
        rare, q_donors = D["rare_type"], D["query_donors"]
        ct_col, donor_col = D["celltype_col"], D["donor_col"]
        ref_genes = set(ref.gene_names)
        ref_types = [str(c) for c in ref.cell_types]
        for seed in args.seeds:
            for scen in SCENARIOS:
                kw = {"rare_type": rare} if scen == "rare" else {}
                if scen == "similar_subtypes":
                    # pick a collinear within-family pair if any
                    fam = {}
                    for t in ref_types:
                        fam.setdefault(mapping.get(t, t), []).append(t)
                    pair = next((m[:2] for m in fam.values() if len(m) >= 2), None)
                    if pair:
                        kw["similar_pair"] = tuple(pair)
                targets = SH.build_target_proportions(ref_types, N_SAMPLES, scen,
                                                      seed=seed, **kw)
                ds_pb = SH.realize_pseudobulk(adata, targets, celltype_col=ct_col,
                                              donor_col=donor_col, query_donors=q_donors,
                                              seed=seed, cells_per_sample=CELLS_PER_SAMPLE)
                counts = ds_pb.counts
                counts = counts.loc[counts.index.intersection(ref_genes)]
                truth = ds_pb.true_mrna_proportions
                genes = sorted(set(counts.index) & ref_genes)
                for method in METHODS:
                    t0 = time.perf_counter()
                    try:
                        with warnings.catch_warnings():
                            warnings.simplefilter("ignore")
                            pred = _run_method(method, counts, ref, genes)
                    except Exception as exc:  # noqa: BLE001
                        runtime_rows.append({"dataset": ds, "seed": seed, "scenario": scen,
                                             "method": method, "runtime_seconds":
                                             round(time.perf_counter() - t0, 2),
                                             "status": f"failed: {exc}"})
                        print(f"  {ds} seed={seed} {scen} {method} FAILED: {exc}")
                        continue
                    rt = float(time.perf_counter() - t0)
                    m = _score(truth, pred, mapping, rare)
                    rows.append({"dataset": ds, "seed": seed, "scenario": scen,
                                 "method": method, **m})
                    runtime_rows.append({"dataset": ds, "seed": seed, "scenario": scen,
                                         "method": method, "runtime_seconds": round(rt, 2),
                                         "status": "ok"})
                    print(f"  {ds} s={seed} {scen:<16} {method:<26} {rt:5.1f}s "
                          f"fine={m['fine_pearson']:.3f} cond={m['conditional_rmse']:.3f} "
                          f"rareP={m['rare_precision']} effN={m['effective_n_pred']:.1f}"
                          f"/{m['effective_n_truth']:.1f}")

    per_scen = pd.DataFrame(rows)
    per_scen.to_csv(OUT / "per_scenario_metrics.tsv", sep="\t", index=False)
    pd.DataFrame(runtime_rows).to_csv(OUT / "runtime_metrics.tsv", sep="\t", index=False)
    mcols = [c for c in per_scen.columns if c not in ("dataset", "seed", "scenario", "method")]
    per_ds = per_scen.groupby(["dataset", "method"])[mcols].mean(numeric_only=True).reset_index()
    per_ds.to_csv(OUT / "per_dataset_metrics.tsv", sep="\t", index=False)
    per_scen.groupby("method")[mcols].mean(numeric_only=True).reset_index().to_csv(
        OUT / "bulk_metrics.tsv", sep="\t", index=False)
    per_scen.groupby(["dataset", "scenario", "method"])[["conditional_rmse"]].mean().reset_index().to_csv(
        OUT / "per_family_metrics.tsv", sep="\t", index=False)

    # promotion gates: each GLM method vs wNNLS, per dataset
    gate_rows = []
    for dsn in per_ds["dataset"].unique():
        d = per_ds[per_ds.dataset == dsn].set_index("method")
        if "wNNLS" not in d.index:
            continue
        base = d.loc["wNNLS"]
        base_effn_err = abs(base.effective_n_pred - base.effective_n_truth)
        for method in ("poisson_glm_experimental", "nb_glm_experimental"):
            if method not in d.index:
                continue
            v = d.loc[method]
            v_effn_err = abs(v.effective_n_pred - v.effective_n_truth)
            def _ge(a, b, tol):  # noqa: E306
                return a >= b - tol * max(abs(b), 1e-9)
            checks = {
                "broad_preserved": _ge(v.broad_pearson, base.broad_pearson, 0.02),
                "fine_preserved": _ge(v.fine_pearson, base.fine_pearson, 0.02),
                "conditional_rmse_not_worse": v.conditional_rmse <= base.conditional_rmse * 1.02
                    if np.isfinite(base.conditional_rmse) else True,
                "rare_precision_not_worse": (not np.isfinite(base.rare_precision)) or
                    (np.isfinite(v.rare_precision) and v.rare_precision >= base.rare_precision - 0.05),
                "rare_sensitivity_not_worse": (not np.isfinite(base.rare_sensitivity)) or
                    (np.isfinite(v.rare_sensitivity) and v.rare_sensitivity >= base.rare_sensitivity - 0.05),
                "fp_subtype_not_increased": v.false_positive_subtype_rate <= base.false_positive_subtype_rate + 1e-9,
                "absent_mass_not_increased": v.absent_subtype_mass <= base.absent_subtype_mass + 1e-9,
                "effn_closer_or_equal": v_effn_err <= base_effn_err + 1e-9,
            }
            gate_rows.append({"dataset": dsn, "method": method,
                              "n_passed": int(sum(checks.values())), "n_gates": len(checks),
                              "all_passed": bool(all(checks.values())),
                              **{f"gate_{k}": bool(x) for k, x in checks.items()}})
    pd.DataFrame(gate_rows).to_csv(OUT / "promotion_gates.tsv", sep="\t", index=False)

    statuses = [{"method": m, "n_runs": int((per_scen["method"] == m).sum()) if not per_scen.empty else 0,
                 "kind": "default" if m == "wNNLS" else "external" if "control" in m else "experimental"}
                for m in METHODS]
    statuses.append({"method": "MuSiC/Bisque/DWLS/SCDC/BayesPrism", "n_runs": 0,
                     "kind": "deferred (not installed in this env)"})
    pd.DataFrame(statuses).to_csv(OUT / "method_status.tsv", sep="\t", index=False)

    (OUT / "manifest.json").write_text(json.dumps({
        "datasets": args.datasets, "seeds": args.seeds, "scenarios": SCENARIOS,
        "methods": METHODS, "default_solver": "wNNLS", "default_changed": False,
        "truth": "true_mrna_proportions", "note": "CARD excluded (spatial)",
    }, indent=2), encoding="utf-8")
    print(f"\nWrote -> {OUT}/")
    if not per_ds.empty:
        show = ["dataset", "method", "fine_pearson", "broad_pearson", "conditional_rmse",
                "rare_precision", "rare_sensitivity", "false_positive_subtype_rate",
                "effective_n_pred", "effective_n_truth"]
        print(per_ds[show].round(3).to_string(index=False))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
