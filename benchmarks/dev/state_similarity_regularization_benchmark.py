#!/usr/bin/env python
"""Benchmark experimental state-similarity / Redeconve-*inspired* regularization.

Real NB-CAR solver fits on breast + HLCA/lung synthetic gold-truth spatial
scenarios (multi-seed), comparing spatial-smoothing baselines (default /
no_smoothing / weak / edge_aware / combined) against the new opt-in, within-family,
post-fit state-similarity refinement (state_regularized / sparsity /
adaptive_resolution) and a state+weak combination. The refinement is NOT a solver
change; it concentrates redundant within-family mass and applies sparsity while
conserving each spot's broad-family mass.

Primary question: does coupling *similar states* (not similar spots) improve the
limits that spatial smoothing could not — conditional within-family RMSE, inflated
effective-N, diffuse false-positive subtype mass — without hurting rare niches?

Usage:
  PYTHONPATH=src:. python benchmarks/diagnostics/state_similarity_regularization_benchmark.py \
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

OUT = REPO / "benchmarks" / "outputs" / "state_similarity_regularization"

# Each config: label -> dict of knobs.
#   lam              : lambda_spatial
#   edge_aware       : in-solver edge-weighted graph
#   sr               : None, or dict(mode, lambda_state, lambda_sparse)
# The shipped presets use small lambdas (state 0.01 / sparse 0.001); we also test
# stronger benchmark-only variants to characterise the mechanism's ceiling. These
# do NOT change the shipped preset defaults.
CONFIGS = {
    "default":               {"lam": 0.1,  "edge_aware": False, "sr": None},
    "no_smoothing":          {"lam": 0.0,  "edge_aware": False, "sr": None},
    "weak_smoothing":        {"lam": 0.02, "edge_aware": False, "sr": None},
    "edge_aware_smoothing":  {"lam": 0.05, "edge_aware": True,  "sr": None},
    "combined_weak_edge":    {"lam": 0.02, "edge_aware": True,  "sr": None},
    "state_regularized":     {"lam": 0.1,  "edge_aware": False,
                              "sr": {"mode": "state_regularized", "lambda_state": 0.01, "lambda_sparse": 0.001}},
    "state_regularized_strong": {"lam": 0.1, "edge_aware": False,
                              "sr": {"mode": "state_regularized", "lambda_state": 0.1, "lambda_sparse": 0.05}},
    "sparsity_state":        {"lam": 0.1,  "edge_aware": False,
                              "sr": {"mode": "sparsity", "lambda_state": 0.0, "lambda_sparse": 0.01}},
    "sparsity_state_strong": {"lam": 0.1,  "edge_aware": False,
                              "sr": {"mode": "sparsity", "lambda_state": 0.0, "lambda_sparse": 0.1}},
    "adaptive_resolution":   {"lam": 0.1,  "edge_aware": False,
                              "sr": {"mode": "adaptive_resolution", "lambda_state": 0.01, "lambda_sparse": 0.001}},
    "state_regularized_weak": {"lam": 0.02, "edge_aware": False,
                              "sr": {"mode": "state_regularized", "lambda_state": 0.1, "lambda_sparse": 0.05}},
}


def _build_cfg(TissueResolveConfig, spec):
    cfg = TissueResolveConfig()
    cfg.spatial_solver.lambda_spatial = float(spec["lam"])
    cfg.spatial_solver.edge_aware = bool(spec["edge_aware"])
    if spec["sr"] is not None:
        cfg.state_regularization.enabled = True
        cfg.state_regularization.mode = spec["sr"]["mode"]
        cfg.state_regularization.lambda_state = float(spec["sr"]["lambda_state"])
        cfg.state_regularization.lambda_sparse = float(spec["sr"]["lambda_sparse"])
        cfg.state_regularization.within_family_only = True
        cfg.state_regularization.preserve_broad_mass = True
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
        compute_state_similarity_graph, recommend_state_groups)

    OUT.mkdir(parents=True, exist_ok=True)
    rows, runtime_rows, graph_rows, group_rows = [], [], [], []

    for ds in args.datasets:
        D = _load_dataset(ds)
        mapping, rare = D["mapping"], D["rare_type"]

        # one state-graph diagnostic per dataset (reference is fixed across seeds)
        try:
            profiles = D["ref"].as_R_cpm()
            cts = [str(c) for c in D["ref"].cell_types]
            g_wf = compute_state_similarity_graph(profiles, cts, family_map=mapping,
                                                  within_family_only=True)
            g_all = compute_state_similarity_graph(profiles, cts, family_map=mapping,
                                                   within_family_only=False)
            cross = sum(1 for e in range(g_all.n_edges)
                        if mapping.get(cts[g_all.state_i[e]]) != mapping.get(cts[g_all.state_j[e]]))
            graph_rows.append({
                "dataset": ds, "n_states": g_wf.n_states,
                "within_family_edges": g_wf.n_edges,
                "cross_family_edges_if_unrestricted": cross,
                "n_components_within_family": g_wf.metadata["n_connected_components"],
                "weight_mean": g_wf.metadata["weight_mean"],
                "weight_median": g_wf.metadata["weight_median"],
                "weight_min": g_wf.metadata["weight_min"],
                "weight_max": g_wf.metadata["weight_max"],
                "n_isolated_states": g_wf.metadata["n_isolated_states"]})
            rec = recommend_state_groups(g_all, family_map=mapping, threshold=0.9,
                                         rare_protection=[rare])
            for grp in rec["recommended_groups"]:
                group_rows.append({"dataset": ds, "family": mapping.get(grp[0], grp[0]),
                                   "group": " | ".join(grp), "size": len(grp)})
        except Exception as exc:  # noqa: BLE001
            print(f"[{ds}] state-graph diagnostic failed: {exc}")

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
                             "lambda_spatial": spec["lam"], "edge_aware": spec["edge_aware"],
                             "state_reg": spec["sr"]["mode"] if spec["sr"] else "none",
                             "broad_mass_dev": md.get("state_reg_broad_mass_max_deviation", 0.0),
                             "mass_shifted_state": md.get("state_reg_mass_shifted_state", 0.0),
                             "mass_shifted_sparse": md.get("state_reg_mass_shifted_sparse", 0.0),
                             **m})
                runtime_rows.append({"dataset": ds, "seed": seed, "method": label,
                                     "runtime_seconds": round(rt, 2), "status": "ok"})
                print(f"  {ds} seed={seed} {label:<26} {rt:5.1f}s fine={m['fine_pearson']:.3f} "
                      f"condRMSE={m['conditional_rmse']:.3f} effN={m['effective_n_pred']:.1f}"
                      f"/{m['effective_n_truth']:.1f} oversm={m['oversmoothing_score']:.2f} "
                      f"rare={m['rare_niche_sensitivity']:.2f} spill={m['pairwise_spillover']:.3f}")

    per_scenario = pd.DataFrame(rows)
    per_scenario.to_csv(OUT / "per_scenario_metrics.tsv", sep="\t", index=False)
    pd.DataFrame(runtime_rows).to_csv(OUT / "runtime_metrics.tsv", sep="\t", index=False)
    pd.DataFrame(graph_rows).to_csv(OUT / "state_graph_summary.tsv", sep="\t", index=False)
    pd.DataFrame(group_rows).to_csv(OUT / "state_group_recommendations.tsv", sep="\t", index=False)

    mcols = [c for c in per_scenario.columns if c not in
             ("dataset", "seed", "method", "lambda_spatial", "edge_aware", "state_reg")]
    per_ds = per_scenario.groupby(["dataset", "method"])[mcols].mean().reset_index()
    per_ds.to_csv(OUT / "per_dataset_metrics.tsv", sep="\t", index=False)
    per_scenario.groupby("method")[mcols].mean().reset_index().to_csv(
        OUT / "spatial_metrics.tsv", sep="\t", index=False)
    per_scenario.groupby("method")[["conditional_rmse"]].mean().reset_index().to_csv(
        OUT / "per_family_metrics.tsv", sep="\t", index=False)

    # ---- promotion gates: each state-reg config vs default and weak, per dataset ----
    gate_rows = []
    sr_methods = [k for k, v in CONFIGS.items() if v["sr"] is not None]
    for ds in per_ds["dataset"].unique():
        d = per_ds[per_ds.dataset == ds].set_index("method")
        if not {"default", "weak_smoothing"} <= set(d.index):
            continue
        base, weak = d.loc["default"], d.loc["weak_smoothing"]
        base_effn_err = abs(base.effective_n_pred - base.effective_n_truth)
        for method in [m for m in sr_methods if m in d.index]:
            v = d.loc[method]
            v_effn_err = abs(v.effective_n_pred - v.effective_n_truth)
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
                "oversmoothing_not_worse_than_default": v.oversmoothing_score <= base.oversmoothing_score + 1e-9,
                "boundary_f1_not_worse_than_weak_5pct": v.boundary_f1 >= weak.boundary_f1 - 0.05 * max(weak.boundary_f1, 1e-9),
                "broad_mass_error_lt_1pct": v.get("broad_mass_dev", 0.0) < 0.01,
            }
            gate_rows.append({"dataset": ds, "method": method,
                              "n_passed": int(sum(checks.values())), "n_gates": len(checks),
                              "all_passed": bool(all(checks.values())),
                              **{f"gate_{k}": bool(x) for k, x in checks.items()}})
    pd.DataFrame(gate_rows).to_csv(OUT / "promotion_gates.tsv", sep="\t", index=False)

    statuses = []
    for lbl, spec in CONFIGS.items():
        n = int((per_scenario["method"] == lbl).sum()) if not per_scenario.empty else 0
        kind = ("diagnostic_only" if (spec["sr"] and spec["sr"]["mode"] == "adaptive_resolution")
                else "state_regularized" if spec["sr"] else "spatial_baseline")
        statuses.append({"method": lbl, "kind": kind, "lambda_spatial": spec["lam"],
                         "n_runs": n, "status": "executed" if n else "no_runs"})
    statuses.append({"method": "external (CARD/RCTD/cell2location)", "kind": "external",
                     "lambda_spatial": "-", "n_runs": 0,
                     "status": "reused prior reports qualitatively; no new reruns"})
    pd.DataFrame(statuses).to_csv(OUT / "method_status.tsv", sep="\t", index=False)

    (OUT / "manifest.json").write_text(json.dumps({
        "datasets": args.datasets, "seeds": args.seeds, "configs": list(CONFIGS),
        "default_lambda": 0.1, "default_changed": False,
        "module": "state_similarity_regularization (experimental, opt-in)",
        "note": "post-fit within-family refinement; not Redeconve; no solver change",
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
