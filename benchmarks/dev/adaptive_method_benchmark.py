#!/usr/bin/env python
"""Screening benchmark for the experimental adaptive methods (gold-truth synthetic).

Runtime-feasible SCREEN (not the full matrix): on the breast synthetic spatial
scenario (validation seeds), it compares the package default (lambda=0.1),
weak_smoothing (0.02), no_smoothing (0.0), and the two post-fit experimental
smoothers (edge_aware, family_adaptive) applied to the default fit. It scores the
same metric panel as the weak-smoothing grid and computes promotion gates vs
default.

The full 2-dataset x 8-scenario x 5-seed x 7-config matrix and external-tool
reruns are DEFERRED for runtime (documented in method_status.tsv, never silently
skipped). Bulk donor-stable weighting / unresolved calibration / auto mode are
pure modules validated by unit tests; their end-to-end pipeline wiring is
deferred and recorded in method_status.tsv.

Usage:
  PYTHONPATH=src:. python benchmarks/diagnostics/adaptive_method_benchmark.py \
      --run-real-data --datasets breast --seeds 0 1
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

OUT = REPO / "benchmarks" / "outputs" / "adaptive_methods"
DEFAULT_LAMBDA = 0.1

# Promotion-gate directions for spatial screening (weakness-targeted).
GATES = {
    "oversmoothing_score": ("lower", 0.10),     # improve target weakness >=10%
    "fine_pearson": ("higher_preserve", 0.0),
    "broad_pearson": ("higher_preserve", 0.0),
    "local_rmse": ("lower_preserve", 0.0),
    "boundary_f1": ("higher_preserve", 0.0),
    "rare_niche_sensitivity": ("higher_preserve", 0.0),
    "false_positive_subtype_rate": ("lower_preserve", 0.0),
}


def _gate_pass(metric, base, val):
    direction, thresh = GATES[metric]
    if direction == "lower":
        return (base - val) / base >= thresh if base > 0 else val <= base
    if direction == "lower_preserve":
        return val <= base + 0.002
    if direction in ("higher_preserve",):
        return val >= base - 0.01
    return True


def main(argv=None) -> int:
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument("--run-real-data", action="store_true")
    ap.add_argument("--datasets", nargs="+", default=["breast"])
    ap.add_argument("--seeds", type=int, nargs="+", default=[0, 1])
    ap.add_argument("--n-side", type=int, default=20)
    args = ap.parse_args(argv)
    if not H.real_data_enabled(args.run_real_data):
        print(f"Refusing without --run-real-data (or {H.REAL_DATA_ENV}=1).", file=sys.stderr)
        return 2

    import tissueresolve as tr
    from tissueresolve.config import TissueResolveConfig
    from tissueresolve.experimental.spatial_adaptive_smoothing import (
        apply_edge_aware_smoothing, family_adaptive_smoothing)

    OUT.mkdir(parents=True, exist_ok=True)
    rows, runtime_rows, status_rows = [], [], []

    def _fit(D, sc, lam):
        cfg = TissueResolveConfig(); cfg.spatial_solver.lambda_spatial = float(lam)
        t0 = time.perf_counter()
        with warnings.catch_warnings():
            warnings.simplefilter("ignore")
            res = tr.deconv_spatial(sc["Y"], D["ref"], sc["array_row"], sc["array_col"],
                                    sc["lib_sizes"], sc["gene_names"], spot_ids=sc["spot_ids"],
                                    config=cfg, resolution_mode="none", run_neighbourhood=False)
        return res.deconv.proportions, float(time.perf_counter() - t0)

    for ds in args.datasets:
        D = _load_dataset(ds)
        coords_attr = None
        for seed in args.seeds:
            sc = SSP.generate_spatial_scenario(
                D["adata"], D["ref_types"], D["query_donors"], D["mapping"],
                celltype_col=D["celltype_col"], donor_col=D["donor_col"],
                n_side=args.n_side, cells_per_spot=40, seed=seed,
                rare_type=D["rare_type"], domain_families=D["domain_families"])
            truth, coords, domains = sc["truth"], sc["coords"], sc["domain_labels"]
            Xspot = sc["Y"]  # spots x genes

            # base fits
            for lam, name in [(0.1, "default"), (0.02, "weak_smoothing"), (0.0, "no_smoothing")]:
                pred, rt = _fit(D, sc, lam)
                m = _score(truth, pred, coords, domains, D["mapping"], D["rare_type"])
                rows.append({"dataset": ds, "seed": seed, "method": name, "lambda_spatial": lam, **m})
                runtime_rows.append({"dataset": ds, "seed": seed, "method": name,
                                     "runtime_seconds": round(rt, 2)})
                if name == "default":
                    default_pred = pred

            # post-fit experimental smoothers on the DEFAULT fit
            t0 = time.perf_counter()
            ea, _meta = apply_edge_aware_smoothing(default_pred, coords, expression=Xspot,
                                                   base_lambda=0.1)
            rows.append({"dataset": ds, "seed": seed, "method": "edge_aware_smoothing",
                         "lambda_spatial": 0.1,
                         **_score(truth, ea, coords, domains, D["mapping"], D["rare_type"])})
            runtime_rows.append({"dataset": ds, "seed": seed, "method": "edge_aware_smoothing",
                                 "runtime_seconds": round(time.perf_counter() - t0, 2)})

            t0 = time.perf_counter()
            fa, _meta2 = family_adaptive_smoothing(default_pred, coords, D["mapping"],
                                                   base_lambda=0.1)
            rows.append({"dataset": ds, "seed": seed, "method": "family_adaptive_smoothing",
                         "lambda_spatial": 0.1,
                         **_score(truth, fa, coords, domains, D["mapping"], D["rare_type"])})
            runtime_rows.append({"dataset": ds, "seed": seed, "method": "family_adaptive_smoothing",
                                 "runtime_seconds": round(time.perf_counter() - t0, 2)})
            print(f"  {ds} seed={seed} done")

    per_scenario = pd.DataFrame(rows)
    per_scenario.to_csv(OUT / "spatial_per_scenario_metrics.tsv", sep="\t", index=False)
    pd.DataFrame(runtime_rows).to_csv(OUT / "runtime_metrics.tsv", sep="\t", index=False)

    metric_cols = [c for c in per_scenario.columns if c not in
                   ("dataset", "seed", "method", "lambda_spatial")]
    agg = per_scenario.groupby("method")[metric_cols].mean().reset_index()
    agg.to_csv(OUT / "spatial_metrics.tsv", sep="\t", index=False)
    # per-family proxy = conditional RMSE per method
    per_scenario.groupby("method")[["conditional_rmse"]].mean().reset_index().to_csv(
        OUT / "spatial_per_family_metrics.tsv", sep="\t", index=False)

    # promotion gates vs default
    a = agg.set_index("method")
    base = a.loc["default"]
    gate_rows = []
    for method in [m for m in a.index if m != "default"]:
        passed = {}
        for metric in GATES:
            if metric in a.columns:
                passed[metric] = bool(_gate_pass(metric, float(base[metric]), float(a.loc[method, metric])))
        gate_rows.append({"method": method, "n_gates": len(passed),
                          "n_passed": int(sum(passed.values())),
                          "all_passed": all(passed.values()),
                          **{f"gate_{k}": v for k, v in passed.items()}})
    pd.DataFrame(gate_rows).to_csv(OUT / "promotion_gates.tsv", sep="\t", index=False)

    # method status (no silent skips; bulk + deferred recorded explicitly)
    status_rows = [
        {"component": "edge_aware_smoothing", "modality": "spatial",
         "evaluation": "screened (breast, post-fit)", "status": "experimental",
         "end_to_end_wired": False, "notes": "post-fit wrapper; not in solver/default"},
        {"component": "family_adaptive_smoothing", "modality": "spatial",
         "evaluation": "screened (breast, post-fit)", "status": "experimental",
         "end_to_end_wired": False, "notes": "post-fit wrapper; not default"},
        {"component": "weak_smoothing", "modality": "spatial",
         "evaluation": "screened (breast)", "status": "experimental(prior)",
         "end_to_end_wired": True, "notes": "config lambda preset; see weak-smoothing report"},
        {"component": "unresolved_calibration", "modality": "bulk/spatial",
         "evaluation": "unit-tested (mass-conserving)", "status": "experimental",
         "end_to_end_wired": False, "notes": "pipeline wiring deferred"},
        {"component": "diagnostic_auto_mode", "modality": "both",
         "evaluation": "unit-tested (advisory)", "status": "experimental",
         "end_to_end_wired": False, "notes": "advisory recommender; does not override user"},
        {"component": "bulk_donor_stable_weighting", "modality": "bulk",
         "evaluation": "unit-tested (gene ranking)", "status": "experimental",
         "end_to_end_wired": False, "notes": "solver wiring deferred; no default change"},
        {"component": "external_spatial_full_matrix", "modality": "spatial",
         "evaluation": "DEFERRED (runtime)", "status": "deferred",
         "end_to_end_wired": False,
         "notes": "full 2-dataset x 8-scenario x 5-seed x 7-config + RCTD/c2l reruns infeasible in-session"},
    ]
    pd.DataFrame(status_rows).to_csv(OUT / "method_status.tsv", sep="\t", index=False)

    # bulk tables: explicit status rows (not silently empty)
    pd.DataFrame([{"method": "bulk_donor_stable_weighting", "status": "module unit-tested; "
                   "end-to-end bulk accuracy benchmark deferred (needs solver wiring)"}]).to_csv(
        OUT / "bulk_metrics.tsv", sep="\t", index=False)
    pd.DataFrame([{"family": "n/a", "status": "deferred (bulk wiring)"}]).to_csv(
        OUT / "bulk_per_family_metrics.tsv", sep="\t", index=False)

    (OUT / "manifest.json").write_text(json.dumps({
        "datasets": args.datasets, "seeds": args.seeds, "default_lambda": DEFAULT_LAMBDA,
        "default_changed": False, "screen_only": True,
    }, indent=2), encoding="utf-8")
    print(f"\nWrote -> {OUT}/")
    print(agg[["method", "fine_pearson", "broad_pearson", "local_rmse",
               "oversmoothing_score", "boundary_f1", "rare_niche_sensitivity",
               "false_positive_subtype_rate"]].round(3).to_string(index=False))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
