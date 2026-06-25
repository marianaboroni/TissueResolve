#!/usr/bin/env python
"""Benchmark combined weak + in-solver edge-aware spatial smoothing (gold-truth).

Real solver fits on breast + HLCA/lung synthetic spatial scenarios (multi-seed),
comparing default / no_smoothing / weak_smoothing / edge_aware_smoothing /
combined_weak_edge_smoothing and combined variants (λ × min_edge_weight). The
combination sets a low global λ AND an edge-weighted spot graph inside the
existing NB-CAR solver (no solver rewrite). Promotion gates vs default,
weak_smoothing, and edge_aware_smoothing.

Usage:
  PYTHONPATH=src:. python benchmarks/diagnostics/combined_spatial_smoothing_benchmark.py \
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

OUT = REPO / "benchmarks" / "outputs" / "combined_spatial_smoothing"

# (label, lambda, edge_aware, min_edge_weight)
CONFIGS = [
    ("default", 0.1, False, 0.05),
    ("no_smoothing", 0.0, False, 0.05),
    ("weak_smoothing", 0.02, False, 0.05),
    ("edge_aware_smoothing", 0.05, True, 0.05),
    ("combined_l0.02_mew0.05", 0.02, True, 0.05),   # the named combined preset
    ("combined_l0.02_mew0.01", 0.02, True, 0.01),
    ("combined_l0.02_mew0.03", 0.02, True, 0.03),
    ("combined_l0.03_mew0.03", 0.03, True, 0.03),
    ("combined_l0.05_mew0.03", 0.05, True, 0.03),
]


def _gate(kind, base, val):
    if kind == "lower20":
        return (base - val) / base >= 0.20 if base > 0 else val <= base
    if kind == "lower_preserve":
        return val <= base + 0.002
    if kind == "within2":
        return val >= base - 0.02 * base
    if kind == "higher_preserve":
        return val >= base - 0.01
    if kind == "lower_strict":
        return val < base
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
            for label, lam, ea, mew in CONFIGS:
                cfg = TissueResolveConfig()
                cfg.spatial_solver.lambda_spatial = float(lam)
                cfg.spatial_solver.edge_aware = bool(ea)
                cfg.spatial_solver.edge_aware_min_weight = float(mew)
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
                rows.append({"dataset": ds, "seed": seed, "method": label, "lambda_spatial": lam,
                             "edge_aware": ea, "min_edge_weight": mew, **m})
                runtime_rows.append({"dataset": ds, "seed": seed, "method": label,
                                     "runtime_seconds": round(rt, 2), "status": status})
                if ea:
                    md = res.deconv.run_metadata
                    edge_rows.append({"dataset": ds, "seed": seed, "method": label,
                                      "min_edge_weight_param": mew,
                                      **{k: md.get(k) for k in
                                         ("edge_weight_min", "edge_weight_max",
                                          "edge_weight_mean", "edge_weight_median")}})
                print(f"  {ds} seed={seed} {label:<24} {rt:5.1f}s fine={m['fine_pearson']:.3f} "
                      f"oversm={m['oversmoothing_score']:.2f} bF1={m['boundary_f1']:.2f} "
                      f"rare={m['rare_niche_sensitivity']:.2f}")

    per_scenario = pd.DataFrame(rows)
    per_scenario.to_csv(OUT / "per_scenario_metrics.tsv", sep="\t", index=False)
    pd.DataFrame(runtime_rows).to_csv(OUT / "runtime_metrics.tsv", sep="\t", index=False)
    pd.DataFrame(edge_rows).to_csv(OUT / "edge_weight_summary.tsv", sep="\t", index=False)

    mcols = [c for c in per_scenario.columns if c not in
             ("dataset", "seed", "method", "lambda_spatial", "edge_aware", "min_edge_weight")]
    per_ds = per_scenario.groupby(["dataset", "method"])[mcols].mean().reset_index()
    per_ds.to_csv(OUT / "per_dataset_metrics.tsv", sep="\t", index=False)
    per_scenario.groupby("method")[mcols].mean().reset_index().to_csv(
        OUT / "spatial_metrics.tsv", sep="\t", index=False)
    per_scenario.groupby("method")[["conditional_rmse"]].mean().reset_index().to_csv(
        OUT / "per_family_metrics.tsv", sep="\t", index=False)

    # promotion gates per dataset vs default, weak, edge_aware
    gate_rows = []
    for ds in per_ds["dataset"].unique():
        d = per_ds[per_ds.dataset == ds].set_index("method")
        if not {"default", "weak_smoothing", "edge_aware_smoothing"} <= set(d.index):
            continue
        base, weak, edge = d.loc["default"], d.loc["weak_smoothing"], d.loc["edge_aware_smoothing"]
        for method in [m for m in d.index if m.startswith("combined")]:
            v = d.loc[method]
            checks = {
                "oversmooth_-20pct_vs_default": _gate("lower20", base.oversmoothing_score, v.oversmoothing_score),
                "oversmooth_better_than_edge_aware": _gate("lower_strict", edge.oversmoothing_score, v.oversmoothing_score),
                "boundary_f1_ge_weak": v.boundary_f1 >= weak.boundary_f1 - 1e-9,
                "boundary_f1_close_to_edge": v.boundary_f1 >= edge.boundary_f1 - 0.01,
                "fine_pearson_preserved_vs_default": _gate("higher_preserve", base.fine_pearson, v.fine_pearson),
                "fine_pearson_within1pct_of_weak": v.fine_pearson >= weak.fine_pearson - 0.01 * max(weak.fine_pearson, 1e-9),
                "broad_pearson_within2": _gate("within2", base.broad_pearson, v.broad_pearson),
                "local_rmse_preserved": _gate("lower_preserve", base.local_rmse, v.local_rmse),
                "rare_sens_preserved": _gate("higher_preserve", base.rare_niche_sensitivity, v.rare_niche_sensitivity),
                "fp_subtype_not_increased": _gate("lower_preserve", base.false_positive_subtype_rate, v.false_positive_subtype_rate),
                "conditional_rmse_le_5pct": v.conditional_rmse <= base.conditional_rmse * 1.05,
            }
            gate_rows.append({"dataset": ds, "method": method,
                              "n_passed": int(sum(checks.values())), "n_gates": len(checks),
                              "all_passed": all(checks.values()),
                              **{f"gate_{k}": bool(x) for k, x in checks.items()}})
    pd.DataFrame(gate_rows).to_csv(OUT / "promotion_gates.tsv", sep="\t", index=False)

    statuses = [{"method": lbl, "lambda": lam, "edge_aware": ea, "min_edge_weight": mew,
                 "n_runs": int((per_scenario["method"] == lbl).sum()),
                 "status": "executed" if (per_scenario["method"] == lbl).any() else "no_runs"}
                for lbl, lam, ea, mew in CONFIGS]
    statuses.append({"method": "external (CARD/RCTD/cell2location)", "lambda": "-", "edge_aware": "-",
                     "min_edge_weight": "-", "n_runs": 0,
                     "status": "reused prior reports; no new reruns (runtime)"})
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
