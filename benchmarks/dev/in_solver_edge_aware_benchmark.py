#!/usr/bin/env python
"""Benchmark in-solver edge-aware spatial smoothing (gold-truth synthetic).

Compares, with REAL solver fits on the breast + HLCA/lung synthetic spatial
scenarios (multi-seed): default (lambda=0.1), weak_smoothing (0.02),
no_smoothing (0.0), and in-solver edge_aware at lambda=0.05 and 0.1. Edge-aware
uses the edge-weighted spot graph inside the existing NB-CAR solver (no solver
rewrite). Scores the full spatial metric panel and computes promotion gates vs
default and vs weak_smoothing.

Usage:
  PYTHONPATH=src:. python benchmarks/diagnostics/in_solver_edge_aware_benchmark.py \
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

OUT = REPO / "benchmarks" / "outputs" / "in_solver_edge_aware"

# (label, lambda, edge_aware)
CONFIGS = [
    ("default", 0.1, False),
    ("weak_smoothing", 0.02, False),
    ("no_smoothing", 0.0, False),
    ("edge_aware_0.05", 0.05, True),
    ("edge_aware_0.10", 0.10, True),
]

GATES = {  # vs default; weakness-targeted (>=20% oversmoothing reduction etc.)
    "oversmoothing_score": ("lower20", None),
    "boundary_f1": ("higher_preserve", None),
    "fine_pearson": ("higher_preserve", None),
    "broad_pearson": ("within2", None),
    "local_rmse": ("lower_preserve", None),
    "rare_niche_sensitivity": ("higher_preserve", None),
    "false_positive_subtype_rate": ("lower_preserve", None),
}


def _gate(metric, base, val):
    kind = GATES[metric][0]
    if kind == "lower20":
        return (base - val) / base >= 0.20 if base > 0 else val <= base
    if kind == "lower_preserve":
        return val <= base + 0.002
    if kind == "within2":
        return val >= base - 0.02 * base
    if kind == "higher_preserve":
        return val >= base - 0.01
    return True


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

    OUT.mkdir(parents=True, exist_ok=True)
    rows, runtime_rows, edge_rows = [], [], []

    for ds in args.datasets:
        D = _load_dataset(ds)
        for seed in args.seeds:
            sc = SSP.generate_spatial_scenario(
                D["adata"], D["ref_types"], D["query_donors"], D["mapping"],
                celltype_col=D["celltype_col"], donor_col=D["donor_col"],
                n_side=args.n_side, cells_per_spot=40, seed=seed,
                rare_type=D["rare_type"], domain_families=D["domain_families"])
            truth, coords, domains = sc["truth"], sc["coords"], sc["domain_labels"]
            for label, lam, ea in CONFIGS:
                cfg = TissueResolveConfig()
                cfg.spatial_solver.lambda_spatial = float(lam)
                cfg.spatial_solver.edge_aware = bool(ea)
                t0 = time.perf_counter()
                status = "ok"
                try:
                    with warnings.catch_warnings():
                        warnings.simplefilter("ignore")
                        res = tr.deconv_spatial(
                            sc["Y"], D["ref"], sc["array_row"], sc["array_col"], sc["lib_sizes"],
                            sc["gene_names"], spot_ids=sc["spot_ids"], config=cfg,
                            resolution_mode="none", run_neighbourhood=False)
                    pred = res.deconv.proportions
                except Exception as exc:  # noqa: BLE001
                    status = f"failed: {exc}"
                    runtime_rows.append({"dataset": ds, "seed": seed, "method": label,
                                         "runtime_seconds": float(time.perf_counter() - t0),
                                         "status": status})
                    print(f"  {ds} seed={seed} {label} FAILED: {exc}")
                    continue
                rt = float(time.perf_counter() - t0)
                m = _score(truth, pred, coords, domains, D["mapping"], D["rare_type"])
                rows.append({"dataset": ds, "seed": seed, "method": label,
                             "lambda_spatial": lam, "edge_aware": ea, **m})
                runtime_rows.append({"dataset": ds, "seed": seed, "method": label,
                                     "runtime_seconds": round(rt, 2), "status": status})
                if ea:
                    md = res.deconv.run_metadata
                    edge_rows.append({"dataset": ds, "seed": seed, "method": label,
                                      **{k: md.get(k) for k in
                                         ("edge_weight_min", "edge_weight_max",
                                          "edge_weight_mean", "edge_weight_median")}})
                print(f"  {ds} seed={seed} {label:<16} {rt:5.1f}s fine_r={m['fine_pearson']:.3f} "
                      f"oversmooth={m['oversmoothing_score']:.2f} boundary_f1={m['boundary_f1']:.2f} "
                      f"rare={m['rare_niche_sensitivity']:.2f}")

    per_scenario = pd.DataFrame(rows)
    per_scenario.to_csv(OUT / "per_scenario_metrics.tsv", sep="\t", index=False)
    pd.DataFrame(runtime_rows).to_csv(OUT / "runtime_metrics.tsv", sep="\t", index=False)
    pd.DataFrame(edge_rows).to_csv(OUT / "edge_weight_summary.tsv", sep="\t", index=False)

    mcols = [c for c in per_scenario.columns
             if c not in ("dataset", "seed", "method", "lambda_spatial", "edge_aware")]
    # per dataset
    per_ds = per_scenario.groupby(["dataset", "method"])[mcols].mean().reset_index()
    per_ds.to_csv(OUT / "per_dataset_metrics.tsv", sep="\t", index=False)
    # overall
    agg = per_scenario.groupby("method")[mcols].mean().reset_index()
    agg.to_csv(OUT / "spatial_metrics.tsv", sep="\t", index=False)
    per_scenario.groupby("method")[["conditional_rmse"]].mean().reset_index().to_csv(
        OUT / "per_family_metrics.tsv", sep="\t", index=False)

    # promotion gates PER DATASET vs default (gate requires both datasets)
    gate_rows = []
    for ds in per_ds["dataset"].unique():
        d = per_ds[per_ds.dataset == ds].set_index("method")
        if "default" not in d.index:
            continue
        base = d.loc["default"]
        for method in [m for m in d.index if m != "default"]:
            passed = {met: bool(_gate(met, float(base[met]), float(d.loc[method, met])))
                      for met in GATES if met in d.columns}
            gate_rows.append({"dataset": ds, "method": method,
                              "n_passed": int(sum(passed.values())), "n_gates": len(passed),
                              "all_passed": all(passed.values()),
                              **{f"gate_{k}": v for k, v in passed.items()}})
    pd.DataFrame(gate_rows).to_csv(OUT / "promotion_gates.tsv", sep="\t", index=False)

    # method status (no silent skips)
    statuses = []
    for label, lam, ea in CONFIGS:
        ran = (per_scenario["method"] == label).sum()
        statuses.append({"method": label, "lambda": lam, "edge_aware": ea,
                         "n_runs": int(ran), "status": "executed" if ran else "no_runs"})
    statuses.append({"method": "external (CARD/RCTD/cell2location)", "lambda": "-",
                     "edge_aware": "-", "n_runs": 0,
                     "status": "reused prior reports; no new reruns (runtime) — see spatial reports"})
    pd.DataFrame(statuses).to_csv(OUT / "method_status.tsv", sep="\t", index=False)

    (OUT / "manifest.json").write_text(json.dumps({
        "datasets": args.datasets, "seeds": args.seeds, "configs": [c[0] for c in CONFIGS],
        "default_lambda": 0.1, "default_changed": False,
    }, indent=2), encoding="utf-8")
    print(f"\nWrote -> {OUT}/")
    show = ["dataset", "method", "fine_pearson", "broad_pearson", "local_rmse",
            "oversmoothing_score", "boundary_f1", "rare_niche_sensitivity",
            "false_positive_subtype_rate"]
    print(per_ds[show].round(3).to_string(index=False))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
