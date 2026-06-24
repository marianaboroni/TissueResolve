#!/usr/bin/env python
"""Benchmark the experimental in-solver state-similarity regularization (Option A).

Compares spatial-smoothing baselines + the previous post-fit state refinement
against the new SEPARATE projected-gradient solver (competition / laplacian state
penalties, with/without sparsity, with/without weak spatial) on breast + HLCA/lung
synthetic gold-truth scenarios (multi-seed). The production NB-CAR solver is used
only as the warm start; it is not modified.

Primary question: does putting the state penalty INSIDE the objective improve what
post-fit redistribution could not — conditional within-family RMSE, effective-N,
spillover, false-positive subtype mass — without hurting rare niches?

Usage:
  PYTHONPATH=src:. python benchmarks/diagnostics/in_solver_state_regularization_benchmark.py \
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
from benchmarks.shared import synthetic_spatial as SSP  # noqa: E402
from benchmarks.spatial.run_weak_smoothing_grid import _load_dataset, _score  # noqa: E402

OUT = REPO / "benchmarks" / "outputs" / "in_solver_state_regularization"

# label -> knobs. spatial baselines + post-fit refinement + in-solver solver modes.
#   lam      : production lambda_spatial
#   edge     : edge-aware spot graph
#   postfit  : post-fit state_regularization mode (None|state_regularized|adaptive_resolution)
#   solver   : in-solver dict(state_penalty, lambda_spatial, lambda_state, lambda_sparse) or None
CONFIGS = {
    "default":               {"lam": 0.1,  "edge": False, "postfit": None, "solver": None},
    "no_smoothing":          {"lam": 0.0,  "edge": False, "postfit": None, "solver": None},
    "weak_smoothing":        {"lam": 0.02, "edge": False, "postfit": None, "solver": None},
    "edge_aware_smoothing":  {"lam": 0.05, "edge": True,  "postfit": None, "solver": None},
    "combined_weak_edge":    {"lam": 0.02, "edge": True,  "postfit": None, "solver": None},
    "postfit_state_regularized": {"lam": 0.1, "edge": False, "postfit": "state_regularized", "solver": None},
    "adaptive_resolution":   {"lam": 0.1,  "edge": False, "postfit": "adaptive_resolution", "solver": None},
    "solver_competition_low":    {"lam": 0.1, "edge": False, "postfit": None,
                                  "solver": {"state_penalty": "competition", "lambda_spatial": 0.02,
                                             "lambda_state": 0.001, "lambda_sparse": 0.0}},
    "solver_competition_medium": {"lam": 0.1, "edge": False, "postfit": None,
                                  "solver": {"state_penalty": "competition", "lambda_spatial": 0.02,
                                             "lambda_state": 0.01, "lambda_sparse": 0.0}},
    "solver_laplacian_low":      {"lam": 0.1, "edge": False, "postfit": None,
                                  "solver": {"state_penalty": "laplacian", "lambda_spatial": 0.02,
                                             "lambda_state": 0.001, "lambda_sparse": 0.0}},
    "solver_laplacian_medium":   {"lam": 0.1, "edge": False, "postfit": None,
                                  "solver": {"state_penalty": "laplacian", "lambda_spatial": 0.02,
                                             "lambda_state": 0.01, "lambda_sparse": 0.0}},
    "solver_competition_sparse": {"lam": 0.1, "edge": False, "postfit": None,
                                  "solver": {"state_penalty": "competition", "lambda_spatial": 0.02,
                                             "lambda_state": 0.005, "lambda_sparse": 0.005}},
    "solver_competition_weak_spatial": {"lam": 0.1, "edge": False, "postfit": None,
                                  "solver": {"state_penalty": "competition", "lambda_spatial": 0.0,
                                             "lambda_state": 0.01, "lambda_sparse": 0.0}},
}


def _build_cfg(TissueResolveConfig, spec):
    cfg = TissueResolveConfig()
    cfg.spatial_solver.lambda_spatial = float(spec["lam"])
    cfg.spatial_solver.edge_aware = bool(spec["edge"])
    if spec["postfit"] is not None:
        cfg.state_regularization.enabled = True
        cfg.state_regularization.mode = spec["postfit"]
    if spec["solver"] is not None:
        s = spec["solver"]
        cfg.state_regularized_solver.enabled = True
        cfg.state_regularized_solver.state_penalty = s["state_penalty"]
        cfg.state_regularized_solver.lambda_spatial = float(s["lambda_spatial"])
        cfg.state_regularized_solver.lambda_state = float(s["lambda_state"])
        cfg.state_regularized_solver.lambda_sparse = float(s["lambda_sparse"])
        cfg.state_regularized_solver.within_family_only = True
        cfg.state_regularized_solver.max_iter = 300
    return cfg


def main(argv=None) -> int:
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument("--run-real-data", action="store_true")
    ap.add_argument("--datasets", nargs="+", default=["breast", "lung"])
    ap.add_argument("--seeds", type=int, nargs="+", default=[0, 1, 2, 3, 4])
    ap.add_argument("--n-side", type=int, default=20)
    args = ap.parse_args(argv)
    if not H.real_data_enabled(args.run_real_data):
        print(f"Refusing without --run-real-data (or {H.REAL_DATA_ENV}=1).", file=sys.stderr)
        return 2

    import tissueresolve as tr
    from tissueresolve.config import TissueResolveConfig
    from tissueresolve.experimental.state_similarity_regularization import (
        compute_state_similarity_graph)

    OUT.mkdir(parents=True, exist_ok=True)
    rows, runtime_rows, graph_rows, loss_rows = [], [], [], []

    for ds in args.datasets:
        D = _load_dataset(ds)
        mapping, rare = D["mapping"], D["rare_type"]
        try:
            g = compute_state_similarity_graph(D["ref"].as_R_cpm(),
                                               [str(c) for c in D["ref"].cell_types],
                                               family_map=mapping, within_family_only=True)
            graph_rows.append({"dataset": ds, "n_states": g.n_states,
                               "within_family_edges": g.n_edges,
                               "n_components": g.metadata["n_connected_components"],
                               "weight_mean": g.metadata["weight_mean"],
                               "weight_median": g.metadata["weight_median"]})
        except Exception as exc:  # noqa: BLE001
            print(f"[{ds}] graph diag failed: {exc}")

        for seed in args.seeds:
            sc = SSP.generate_spatial_scenario(
                D["adata"], D["ref_types"], D["query_donors"], mapping,
                celltype_col=D["celltype_col"], donor_col=D["donor_col"],
                n_side=args.n_side, cells_per_spot=40, seed=seed,
                rare_type=rare, domain_families=D["domain_families"])
            truth, coords, domains = sc["truth"], sc["coords"], sc["domain_labels"]
            for label, spec in CONFIGS.items():
                cfg = _build_cfg(TissueResolveConfig, spec)
                t0 = time.perf_counter()
                try:
                    with warnings.catch_warnings():
                        warnings.simplefilter("ignore")
                        res = tr.deconv_spatial(
                            sc["Y"], D["ref"], sc["array_row"], sc["array_col"], sc["lib_sizes"],
                            sc["gene_names"], spot_ids=sc["spot_ids"], config=cfg,
                            resolution_mode="none", run_neighbourhood=False,
                            family_map=mapping, rare_protection=[rare])
                    pred = res.deconv.proportions
                except Exception as exc:  # noqa: BLE001
                    runtime_rows.append({"dataset": ds, "seed": seed, "method": label,
                                         "runtime_seconds": round(time.perf_counter() - t0, 2),
                                         "status": f"failed: {exc}"})
                    print(f"  {ds} seed={seed} {label} FAILED: {exc}")
                    continue
                rt = float(time.perf_counter() - t0)
                m = _score(truth, pred, coords, domains, mapping, rare)
                md = res.deconv.run_metadata
                rows.append({"dataset": ds, "seed": seed, "method": label,
                             "solver_used": bool(spec["solver"]),
                             "n_iter": md.get("state_solver_n_iter"),
                             "converged": md.get("state_solver_converged"),
                             "final_loss": md.get("state_solver_final_loss"),
                             "recon_loss": md.get("state_solver_recon_loss"),
                             "loss_monotonic": md.get("state_solver_loss_monotonic"),
                             "broad_mass_dev": md.get("state_reg_broad_mass_max_deviation", 0.0),
                             **m})
                runtime_rows.append({"dataset": ds, "seed": seed, "method": label,
                                     "runtime_seconds": round(rt, 2), "status": "ok"})
                if spec["solver"]:
                    loss_rows.append({"dataset": ds, "seed": seed, "method": label,
                                      "n_iter": md.get("state_solver_n_iter"),
                                      "converged": md.get("state_solver_converged"),
                                      "loss_monotonic": md.get("state_solver_loss_monotonic"),
                                      "final_loss": md.get("state_solver_final_loss"),
                                      "recon_loss": md.get("state_solver_recon_loss"),
                                      "spatial_penalty": md.get("state_solver_spatial_penalty"),
                                      "state_penalty_value": md.get("state_solver_state_penalty_value"),
                                      "sparsity_penalty": md.get("state_solver_sparsity_penalty")})
                print(f"  {ds} seed={seed} {label:<30} {rt:5.1f}s fine={m['fine_pearson']:.3f} "
                      f"condRMSE={m['conditional_rmse']:.3f} effN={m['effective_n_pred']:.1f}"
                      f"/{m['effective_n_truth']:.1f} oversm={m['oversmoothing_score']:.2f} "
                      f"rare={m['rare_niche_sensitivity']:.2f} spill={m['pairwise_spillover']:.3f}")

    per_scenario = pd.DataFrame(rows)
    per_scenario.to_csv(OUT / "per_scenario_metrics.tsv", sep="\t", index=False)
    pd.DataFrame(runtime_rows).to_csv(OUT / "runtime_metrics.tsv", sep="\t", index=False)
    pd.DataFrame(graph_rows).to_csv(OUT / "state_graph_summary.tsv", sep="\t", index=False)
    pd.DataFrame(loss_rows).to_csv(OUT / "loss_curves.tsv", sep="\t", index=False)

    mcols = [c for c in per_scenario.columns if c not in
             ("dataset", "seed", "method", "solver_used", "converged", "loss_monotonic")]
    per_ds = per_scenario.groupby(["dataset", "method"])[mcols].mean(numeric_only=True).reset_index()
    per_ds.to_csv(OUT / "per_dataset_metrics.tsv", sep="\t", index=False)
    per_scenario.groupby("method")[mcols].mean(numeric_only=True).reset_index().to_csv(
        OUT / "spatial_metrics.tsv", sep="\t", index=False)
    per_scenario.groupby("method")[["conditional_rmse"]].mean().reset_index().to_csv(
        OUT / "per_family_metrics.tsv", sep="\t", index=False)

    # promotion gates: each solver config vs default and weak, per dataset
    gate_rows = []
    solver_methods = [k for k, v in CONFIGS.items() if v["solver"] is not None]
    for ds in per_ds["dataset"].unique():
        d = per_ds[per_ds.dataset == ds].set_index("method")
        if not {"default", "weak_smoothing"} <= set(d.index):
            continue
        base, weak = d.loc["default"], d.loc["weak_smoothing"]
        base_effn_err = abs(base.effective_n_pred - base.effective_n_truth)
        for method in [m for m in solver_methods if m in d.index]:
            v = d.loc[method]
            v_effn_err = abs(v.effective_n_pred - v.effective_n_truth)
            conv = per_scenario[(per_scenario.dataset == ds) & (per_scenario.method == method)]["converged"]
            checks = {
                "cond_rmse_improved_5pct": v.conditional_rmse <= base.conditional_rmse * 0.95,
                "cond_rmse_not_worse_2pct": v.conditional_rmse <= base.conditional_rmse * 1.02,
                "effn_closer_to_truth": v_effn_err <= base_effn_err + 1e-9,
                "fp_subtype_not_increased": v.false_positive_subtype_rate <= base.false_positive_subtype_rate + 1e-9,
                "spillover_not_increased": v.pairwise_spillover <= base.pairwise_spillover + 1e-9,
                "rare_sens_preserved": v.rare_niche_sensitivity >= base.rare_niche_sensitivity - 1e-9,
                "fine_pearson_within1pct": v.fine_pearson >= base.fine_pearson - 0.01 * max(base.fine_pearson, 1e-9),
                "broad_pearson_within2pct": v.broad_pearson >= base.broad_pearson - 0.02 * max(base.broad_pearson, 1e-9),
                "local_rmse_preserved": v.local_rmse <= base.local_rmse + 0.002,
                "oversmoothing_not_worse_weak_10pct": v.oversmoothing_score <= weak.oversmoothing_score * 1.10,
                "boundary_f1_not_worse_weak_5pct": v.boundary_f1 >= weak.boundary_f1 - 0.05 * max(weak.boundary_f1, 1e-9),
                "broad_mass_dev_lt_1pct": v.get("broad_mass_dev", 0.0) < 0.01,
                "optimizer_converged": bool(conv.fillna(False).mean() >= 0.5),
            }
            gate_rows.append({"dataset": ds, "method": method,
                              "n_passed": int(sum(checks.values())), "n_gates": len(checks),
                              "all_passed": bool(all(checks.values())),
                              **{f"gate_{k}": bool(x) for k, x in checks.items()}})
    pd.DataFrame(gate_rows).to_csv(OUT / "promotion_gates.tsv", sep="\t", index=False)

    statuses = []
    for lbl, spec in CONFIGS.items():
        n = int((per_scenario["method"] == lbl).sum()) if not per_scenario.empty else 0
        kind = ("in_solver" if spec["solver"] else "post_fit" if spec["postfit"] else "spatial_baseline")
        statuses.append({"method": lbl, "kind": kind, "n_runs": n,
                         "status": "executed" if n else "no_runs"})
    statuses.append({"method": "external (CARD/RCTD/cell2location)", "kind": "external", "n_runs": 0,
                     "status": "reused prior reports qualitatively; no new reruns"})
    pd.DataFrame(statuses).to_csv(OUT / "method_status.tsv", sep="\t", index=False)

    (OUT / "manifest.json").write_text(json.dumps({
        "datasets": args.datasets, "seeds": args.seeds, "configs": list(CONFIGS),
        "default_lambda": 0.1, "default_changed": False,
        "module": "in_solver state_regularized_solver_experimental (Option A; separate solver)",
        "note": "production NB-CAR solver untouched; warm-started joint projected-gradient; not Redeconve",
    }, indent=2), encoding="utf-8")
    print(f"\nWrote -> {OUT}/")
    if not per_ds.empty:
        show = ["dataset", "method", "fine_pearson", "broad_pearson", "conditional_rmse",
                "effective_n_pred", "effective_n_truth", "oversmoothing_score",
                "rare_niche_sensitivity", "false_positive_subtype_rate", "pairwise_spillover"]
        print(per_ds[show].round(3).to_string(index=False))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
