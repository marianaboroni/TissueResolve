#!/usr/bin/env python
"""Family learnability benchmark: is the within-family fine-state composition
LEARNABLE from bulk on donor-held-out pseudobulks, vs mean-reference baselines?

For each broad family (with >=3 fine states) on breast (+ lung if available), fit the
supervised conditional estimator (donor-held-out) and compare its conditional
within-family RMSE / rare behaviour against:
  - mean-reference conditional split (uniform NNLS)        [baseline]
  - discriminative-weighted NNLS                            [baseline]
  - supervised: ridge / elastic_net / pairwise_ridge        [candidate]
across >=5 seeds. Verdict per family: learnable / partially_learnable / not_learnable /
unstable / diagnostic_only.

Experimental; outputs gitignored; no defaults changed; not wired into the pipeline.

Usage:
  PYTHONPATH=src:. python benchmarks/dev/family_learnability_benchmark.py \
      --run-real-data --datasets breast --seeds 0 1 2 3 4
"""
from __future__ import annotations

import argparse
import json
import sys
import time
from pathlib import Path

import numpy as np
import pandas as pd

REPO = Path(__file__).resolve().parents[2]
sys.path.insert(0, str(REPO / "examples" / "real_breast_cancer" / "scripts"))
sys.path.insert(0, str(REPO))
sys.path.insert(0, str(REPO / "src"))
import _harness as H  # noqa: E402
from benchmarks.spatial.run_weak_smoothing_grid import _load_dataset  # noqa: E402

OUT = REPO / "benchmarks" / "outputs" / "family_learnability"
MODELS = ["ridge", "elastic_net", "pairwise_ridge"]
MIN_STATES = 3
MIN_FAMILY_CELLS = 200


def main(argv=None) -> int:
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument("--run-real-data", action="store_true")
    ap.add_argument("--datasets", nargs="+", default=["breast"])
    ap.add_argument("--seeds", type=int, nargs="+", default=[0, 1, 2, 3, 4])
    ap.add_argument("--feature-mode", default="combined")
    args = ap.parse_args(argv)
    if not H.real_data_enabled(args.run_real_data):
        print(f"Refusing without --run-real-data (or {H.REAL_DATA_ENV}=1).", file=sys.stderr)
        return 2

    from tissueresolve.experimental.conditional_family_estimator import (
        fit_conditional_family_estimator, DEFAULT_GATES)

    OUT.mkdir(parents=True, exist_ok=True)
    per_seed, fail = [], []
    for ds in args.datasets:
        D = _load_dataset(ds)
        adata, mapping, ctc, dc = D["adata"], D["mapping"], D["celltype_col"], D["donor_col"]
        # derive family column from mapping (fine -> broad)
        adata.obs = adata.obs.copy()
        adata.obs["__family"] = adata.obs[ctc].astype(str).map(lambda s: str(mapping.get(s, s)))
        fam_counts = adata.obs.groupby("__family")[ctc].nunique()
        cell_counts = adata.obs["__family"].value_counts()
        families = [f for f in fam_counts.index
                    if fam_counts[f] >= MIN_STATES and cell_counts.get(f, 0) >= MIN_FAMILY_CELLS]
        print(f"[{ds}] families with >={MIN_STATES} states & >={MIN_FAMILY_CELLS} cells: {families}")
        for fam in families:
            for model in MODELS:
                for seed in args.seeds:
                    t0 = time.perf_counter()
                    try:
                        m = fit_conditional_family_estimator(
                            adata, "__family", ctc, dc, fam, feature_mode=args.feature_mode,
                            model=model, n_train_mixtures=800, n_val_mixtures=250,
                            random_state=seed)
                    except Exception as exc:  # noqa: BLE001
                        fail.append({"dataset": ds, "family": fam, "model": model,
                                     "seed": seed, "error": str(exc)})
                        continue
                    rt = time.perf_counter() - t0
                    v = m.validation
                    per_seed.append({
                        "dataset": ds, "family": fam, "n_states": v["n_states"], "model": model,
                        "seed": seed, "rmse_model": v["rmse_model"], "rmse_baseline": v["rmse_baseline"],
                        "delta_rmse_pct": v["delta_rmse_pct"], "pearson": v["pearson"],
                        "rare_fpr_model": v.get("rare_fpr_model"), "rare_fpr_baseline": v.get("rare_fpr_baseline"),
                        "calibration_error": v["calibration_error"], "decision_seed": m.decision,
                        "runtime_s": round(rt, 2)})
                    print(f"  {ds} {fam:<14} {model:<14} s{seed} rmse {v['rmse_baseline']:.3f}->"
                          f"{v['rmse_model']:.3f} (Δ{v['delta_rmse_pct']:+.0f}%) "
                          f"r={v['pearson']:.2f} dec={m.decision}")

    ps = pd.DataFrame(per_seed)
    ps.to_csv(OUT / "per_seed_metrics.tsv", sep="\t", index=False)
    if not fail:
        pd.DataFrame([{"status": "no failures"}]).to_csv(OUT / "model_status.tsv", sep="\t", index=False)
    else:
        pd.DataFrame(fail).to_csv(OUT / "model_status.tsv", sep="\t", index=False)

    # aggregate per (dataset, family, model): mean metrics + reproducibility
    gates = DEFAULT_GATES
    rows = []
    for (ds, fam, model), g in ps.groupby(["dataset", "family", "model"]):
        delta = g["delta_rmse_pct"].to_numpy()
        rmse_cv = float(g["rmse_model"].std() / max(g["rmse_model"].mean(), 1e-9))
        n_improve = int((delta >= gates["min_delta_rmse_pct"]).sum())
        fpr_inc = float((g["rare_fpr_model"].fillna(0) - g["rare_fpr_baseline"].fillna(0)).mean())
        cal = float(g["calibration_error"].mean())
        # family-level decision
        if cal > gates["max_calibration_error"] or fpr_inc > gates["max_rare_fpr_increase"]:
            dec = "diagnostic_only"
        elif rmse_cv > gates["max_seed_instability"]:
            dec = "unstable"
        elif n_improve >= 4 and delta.mean() >= gates["min_delta_rmse_pct"]:
            dec = "learnable"
        elif delta.mean() > gates["partial_delta_rmse_pct"]:
            dec = "partially_learnable"
        else:
            dec = "not_learnable"
        rows.append({"dataset": ds, "family": fam, "n_states": int(g["n_states"].iloc[0]),
                     "model": model, "rmse_baseline": round(g["rmse_baseline"].mean(), 4),
                     "rmse_model": round(g["rmse_model"].mean(), 4),
                     "delta_rmse_pct": round(delta.mean(), 1), "n_seeds_improved": n_improve,
                     "pearson": round(g["pearson"].mean(), 3), "rare_fpr_increase": round(fpr_inc, 3),
                     "calibration_error": round(cal, 4), "rmse_cv_seeds": round(rmse_cv, 3),
                     "runtime_s": round(g["runtime_s"].mean(), 2), "decision": dec})
    fm = pd.DataFrame(rows).sort_values(["dataset", "family", "delta_rmse_pct"], ascending=[True, True, False])
    fm.to_csv(OUT / "per_family_metrics.tsv", sep="\t", index=False)
    # best model per family
    best = fm.sort_values("delta_rmse_pct", ascending=False).groupby(["dataset", "family"]).head(1)
    best.to_csv(OUT / "family_learnability_summary.tsv", sep="\t", index=False)

    verdict = {"families": {}, "any_learnable": False}
    for _, r in best.iterrows():
        verdict["families"][f"{r['dataset']}:{r['family']}"] = {
            "best_model": r["model"], "delta_rmse_pct": r["delta_rmse_pct"],
            "decision": r["decision"], "rare_fpr_increase": r["rare_fpr_increase"]}
        if r["decision"] == "learnable":
            verdict["any_learnable"] = True
    (OUT / "verdict.json").write_text(json.dumps(verdict, indent=2), encoding="utf-8")
    (OUT / "manifest.json").write_text(json.dumps({
        "datasets": args.datasets, "seeds": args.seeds, "models": MODELS,
        "feature_mode": args.feature_mode, "gates": gates,
        "note": "supervised conditional family estimator vs mean-NNLS baseline; "
                "donor-held-out; experimental; defaults unchanged"}, indent=2), encoding="utf-8")
    print(f"\nWrote -> {OUT}/")
    show = ["dataset", "family", "n_states", "model", "rmse_baseline", "rmse_model",
            "delta_rmse_pct", "n_seeds_improved", "rare_fpr_increase", "decision"]
    print("\n=== best model per family ===")
    print(best[show].to_string(index=False))
    print("\nVERDICT:", json.dumps(verdict["families"], indent=2))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
