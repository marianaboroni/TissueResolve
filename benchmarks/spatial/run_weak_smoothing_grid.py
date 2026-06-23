#!/usr/bin/env python
"""Experimental spatial weak-smoothing lambda grid (gold-truth synthetic).

Sweeps the CAR spatial smoothing strength lambda_spatial over a grid on the SAME
synthetic spatial scenario used by run_synthetic_spatial.py, replicated across
seeds (each seed = a distinct realisation of the structured scenario: sharp
border + gradient + rare niche + mixed/collinear spots). Scores accuracy,
spatial-structure, composition, and runtime for each (lambda, seed).

Default behaviour is NOT changed: lambda=0.1 is the package default and is one of
the grid points; lower lambdas are evaluated only here, on synthetic scenarios
(never on a final/real test set). No solver added, soft gating/unresolved mass
untouched.

Usage:  PYTHONPATH=src:. python benchmarks/spatial/run_weak_smoothing_grid.py --run-real-data
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
from benchmarks.shared import synthetic_spatial as SSP  # noqa: E402
from benchmarks.shared import metrics as M  # noqa: E402
from benchmarks.shared import spatial_metrics as SM  # noqa: E402

OUT = REPO / "benchmarks" / "outputs" / "spatial_weak_smoothing"
HIERARCHY_TSV = REPO / "examples" / "real_breast_cancer" / "config" / \
    "breast_cancer_cell_type_hierarchy.tsv"

LAMBDA_GRID = [0.0, 0.01, 0.02, 0.03, 0.05, 0.1]   # 0.1 = package default
DEFAULT_LAMBDA = 0.1
SEEDS = [0, 1, 2]
N_SIDE = 20
CELLS_PER_SPOT = 40
MIN_CELLS = 30


def _conditional(df, truth_cols, mapping):
    """Within-family conditional proportions over the given fine columns."""
    fam_members: dict[str, list[str]] = {}
    for c in truth_cols:
        fam_members.setdefault(str(mapping.get(c, c)), []).append(c)
    out = pd.DataFrame(0.0, index=df.index, columns=truth_cols)
    for fam, members in fam_members.items():
        members = [m for m in members if m in df.columns]
        if len(members) < 2:
            continue
        s = df[members].sum(axis=1).replace(0, np.nan)
        for m in members:
            out[m] = (df[m] / s).fillna(0.0)
    multi = [c for fam, mem in fam_members.items() if len(mem) > 1 for c in mem]
    return out, multi


def _score(truth, pred, coords, domains, mapping, rare_type):
    cols = [c for c in truth.columns if not str(c).startswith("unresolved_")]
    t = truth[cols]
    p = pred.reindex(index=t.index, columns=cols).fillna(0.0)
    acc = M.accuracy_metrics(t, p)
    comp = M.compositional_metrics(t, p)
    tfam, pfam = M.aggregate_to_families(t, mapping), M.aggregate_to_families(p, mapping)
    bacc = M.accuracy_metrics(tfam, pfam)
    fid = SM.spatial_fidelity_metrics(t, p, coords, domain_labels=domains, k=6)
    bnd = SM.boundary_metrics(t, p, coords, domains, mapping=mapping, k=6)

    # conditional within-family RMSE (fine subtypes that share a family)
    tc, multi = _conditional(t, cols, mapping)
    pc, _ = _conditional(p, cols, mapping)
    cond_rmse = float(np.sqrt(np.mean(
        (tc[multi].to_numpy(float) - pc[multi].to_numpy(float)) ** 2))) if multi else float("nan")

    # rare niche detection in niche spots
    dom = np.asarray(domains)
    niche = dom == "niche"
    rare_sens = SM.rare_niche_sensitivity(t, p, domains, rare_type)
    if rare_type in p.columns and niche.any():
        pred_pos = p[rare_type].to_numpy(float) > 0.01
        tp = int((pred_pos & niche).sum()); fp = int((pred_pos & ~niche).sum())
        fn = int((~pred_pos & niche).sum()); tn = int((~pred_pos & ~niche).sum())
        rare_prec = tp / (tp + fp) if (tp + fp) else float("nan")
        rare_fpr = fp / (fp + tn) if (fp + tn) else float("nan")
    else:
        rare_prec = rare_fpr = float("nan")

    # composition behaviour
    ent_pred = float(M.proportion_entropy(p).mean())
    ent_truth = float(M.proportion_entropy(t).mean())
    effn_pred = float(np.mean([M.effective_n_populations(p.iloc[i].to_numpy()) for i in range(p.shape[0])]))
    effn_truth = float(np.mean([M.effective_n_populations(t.iloc[i].to_numpy()) for i in range(t.shape[0])]))
    sparsity = M.near_zero_fraction(p, eps=1e-3)
    dominant_frac = float(p.max(axis=1).mean())
    absent_mask = t.to_numpy(float) < 1e-9
    absent_mass = float(np.where(absent_mask, p.to_numpy(float), 0.0).sum(axis=1).mean())
    fp_subtype_rate = float(((p.to_numpy(float) > 0.01) & absent_mask).mean())
    # pairwise spillover: mean predicted mass placed on truly-absent types per spot
    pairwise_spillover = absent_mass

    return {
        "broad_pearson": bacc["pearson"], "broad_rmse": bacc["rmse"],
        "fine_pearson": acc["pearson"], "fine_rmse": acc["rmse"],
        "conditional_rmse": cond_rmse, "local_rmse": fid["local_rmse"],
        "dominant_accuracy": M.dominant_accuracy(t, p),
        "oversmoothing_score": fid["oversmoothing_score"],
        "boundary_f1": bnd["boundary_f1"], "edge_blurring": bnd["edge_blurring"],
        "morans_corr_true_pred": fid.get("morans_i_corr_true_pred", float("nan")),
        "morans_mae": fid.get("morans_i_mae", float("nan")),
        "domain_recovery_ari": fid.get("domain_recovery_ari", float("nan")),
        "rare_niche_sensitivity": rare_sens, "rare_niche_precision": rare_prec,
        "rare_false_positive_rate": rare_fpr,
        "entropy_pred": ent_pred, "entropy_truth": ent_truth,
        "effective_n_pred": effn_pred, "effective_n_truth": effn_truth,
        "sparsity_near_zero_frac": sparsity, "dominant_fraction": dominant_frac,
        "absent_subtype_mass": absent_mass,
        "false_positive_subtype_rate": fp_subtype_rate,
        "pairwise_spillover": pairwise_spillover,
    }


def main(argv=None) -> int:
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument("--run-real-data", action="store_true")
    ap.add_argument("--n-side", type=int, default=N_SIDE)
    ap.add_argument("--seeds", type=int, nargs="+", default=SEEDS)
    ap.add_argument("--lambdas", type=float, nargs="+", default=LAMBDA_GRID)
    args = ap.parse_args(argv)
    if not H.real_data_enabled(args.run_real_data):
        print(f"Refusing without --run-real-data (or {H.REAL_DATA_ENV}=1).", file=sys.stderr)
        return 2

    import anndata as ad
    import tissueresolve as tr
    from tissueresolve.config import TissueResolveConfig
    from tissueresolve.reference.hierarchy import load_hierarchy_mapping, build_cell_type_hierarchy

    OUT.mkdir(parents=True, exist_ok=True)
    adata = ad.read_h5ad(H.REFERENCE_H5AD)
    if "feature_name" in adata.var.columns:
        fn = adata.var["feature_name"].astype(str)
        if fn.nunique() == adata.n_vars:
            adata.var_names = fn.to_numpy()
    # Reference built ONCE from a fixed donor split; only the spatial scenario
    # realisation varies across seeds (isolates spatial-structure variation).
    ref_donors, query_donors = SH.split_donors(adata, "donor_id", ref_frac=0.5, seed=0)
    ref_mask = adata.obs["donor_id"].astype(str).isin(set(ref_donors)).to_numpy()
    with warnings.catch_warnings():
        warnings.simplefilter("ignore")
        ref = H.prepare_reference(adata[ref_mask].copy(), min_cells=MIN_CELLS,
                                  estimate_overdispersion=True).reference
    ref_types = [str(c) for c in ref.cell_types]
    mapping = build_cell_type_hierarchy(ref_types, load_hierarchy_mapping(HIERARCHY_TSV))
    rare_type = next((t for t in ["regulatory T cell", "natural killer cell"]
                      if t in ref_types), ref_types[-1])
    print(f"reference: {len(ref_types)} types; rare={rare_type}; "
          f"lambdas={args.lambdas}; seeds={args.seeds}")

    per_scenario, runtime_rows, per_family_rows = [], [], []
    for seed in args.seeds:
        sc = SSP.generate_spatial_scenario(adata, ref_types, query_donors, mapping,
                                           n_side=args.n_side, cells_per_spot=CELLS_PER_SPOT,
                                           seed=seed, rare_type=rare_type)
        truth, coords, domains = sc["truth"], sc["coords"], sc["domain_labels"]
        for lam in args.lambdas:
            cfg = TissueResolveConfig()
            cfg.spatial_solver.lambda_spatial = float(lam)
            t0 = time.perf_counter()
            status, n_warn = "ok", 0
            try:
                with warnings.catch_warnings(record=True) as caught:
                    warnings.simplefilter("always")
                    res = tr.deconv_spatial(
                        sc["Y"], ref, sc["array_row"], sc["array_col"], sc["lib_sizes"],
                        sc["gene_names"], spot_ids=sc["spot_ids"], config=cfg,
                        resolution_mode="none", run_neighbourhood=False)
                    n_warn = len(caught)
                pred = res.deconv.proportions
            except Exception as exc:  # noqa: BLE001
                print(f"  seed={seed} lambda={lam} FAILED: {exc}")
                runtime_rows.append({"seed": seed, "lambda_spatial": lam,
                                     "runtime_seconds": float(time.perf_counter() - t0),
                                     "n_warnings": n_warn, "status": "failed"})
                continue
            rt = float(time.perf_counter() - t0)
            preset = ("default" if lam == DEFAULT_LAMBDA else
                      "weak_smoothing" if lam == 0.02 else
                      "no_smoothing" if lam == 0.0 else f"grid({lam})")
            m = _score(truth, pred, coords, domains, mapping, rare_type)
            row = {"seed": seed, "lambda_spatial": lam, "preset": preset,
                   "is_default": lam == DEFAULT_LAMBDA, "runtime_seconds": round(rt, 2),
                   "n_warnings": n_warn, **m}
            per_scenario.append(row)
            runtime_rows.append({"seed": seed, "lambda_spatial": lam,
                                 "runtime_seconds": round(rt, 2), "n_warnings": n_warn,
                                 "status": "ok"})
            print(f"  seed={seed} lambda={lam:<4} {rt:5.1f}s fine_r={m['fine_pearson']:.3f} "
                  f"oversmooth={m['oversmoothing_score']:.2f} boundary_f1={m['boundary_f1']:.2f} "
                  f"rare_sens={m['rare_niche_sensitivity']:.2f}")

    ps = pd.DataFrame(per_scenario)
    ps.to_csv(OUT / "per_scenario_metrics.tsv", sep="\t", index=False)
    pd.DataFrame(runtime_rows).to_csv(OUT / "runtime_metrics.tsv", sep="\t", index=False)

    # lambda grid = mean over seeds
    metric_cols = [c for c in ps.columns if c not in
                   ("seed", "lambda_spatial", "preset", "is_default")]
    grid = ps.groupby("lambda_spatial")[metric_cols].mean().reset_index()
    grid["preset"] = grid["lambda_spatial"].map(
        lambda l: "default" if l == DEFAULT_LAMBDA else "weak_smoothing" if l == 0.02
        else "no_smoothing" if l == 0.0 else f"grid({l})")
    grid["n_seeds"] = ps.groupby("lambda_spatial").size().values
    grid.to_csv(OUT / "lambda_grid_metrics.tsv", sep="\t", index=False)

    # per-family conditional RMSE per lambda (uses the mapping families with >1 member)
    pf = (ps.groupby("lambda_spatial")[["conditional_rmse"]].mean()
          .reset_index().rename(columns={"conditional_rmse": "mean_conditional_rmse"}))
    pf.to_csv(OUT / "per_family_metrics.tsv", sep="\t", index=False)

    (OUT / "manifest.json").write_text(json.dumps({
        "lambda_grid": args.lambdas, "default_lambda": DEFAULT_LAMBDA,
        "seeds": args.seeds, "n_side": args.n_side, "cells_per_spot": CELLS_PER_SPOT,
        "rare_type": rare_type, "dataset": "breast (real_breast_cancer reference)",
        "scenario": "sharp border + gradient + rare niche + mixed/collinear spots",
        "ground_truth": "mRNA-proportion per spot", "default_changed": False,
        "experimental_preset": "weak_smoothing (lambda=0.02)",
    }, indent=2), encoding="utf-8")
    print(f"\nWrote -> {OUT}/")
    show = ["lambda_spatial", "preset", "fine_pearson", "broad_pearson", "local_rmse",
            "oversmoothing_score", "boundary_f1", "rare_niche_sensitivity",
            "false_positive_subtype_rate", "runtime_seconds"]
    print(grid[show].round(3).to_string(index=False))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
