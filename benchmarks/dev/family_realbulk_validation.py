#!/usr/bin/env python
"""Real-bulk-realism validation gate for the conditional family estimator.

The estimator passed donor-held-out *clean same-platform* pseudobulk. This gate adds
the confounders real bulk has, on the families found learnable:

  Gate A — CROSS-PLATFORM: breast atlas has two chemistries (10x 3' v3 / v2). Train the
           estimator on v3 donors, validate on v2-donor pseudobulk (protocol shift +
           donor shift). Scoreable (pseudobulk truth known).
  Gate B — TECHNICAL CONFOUNDERS: same-platform donor-held-out, but the validation
           bulks get library-size shift + multiplicative noise + gene dropout.

For each (gate, family, model, seed): supervised conditional RMSE vs the mean-NNLS
baseline on the SAME shifted validation bulks (delta %), rare FPR change, decision.

(Real TCGA-TNBC bulk is present but has NO fine-state ground truth and would require
full deconvolution to isolate a family signal → conditional accuracy there is DEFERRED,
not fabricated.)

Experimental; outputs gitignored; no defaults changed.

Usage:
  PYTHONPATH=src:. python benchmarks/dev/family_realbulk_validation.py \
      --run-real-data --seeds 0 1 2 3 4
"""
from __future__ import annotations

import argparse
import json
import sys
from pathlib import Path

import numpy as np
import pandas as pd

REPO = Path(__file__).resolve().parents[2]
sys.path.insert(0, str(REPO / "examples" / "real_breast_cancer" / "scripts"))
sys.path.insert(0, str(REPO))
sys.path.insert(0, str(REPO / "src"))
import _harness as H  # noqa: E402
from benchmarks.spatial.run_weak_smoothing_grid import _load_dataset  # noqa: E402

OUT = REPO / "benchmarks" / "outputs" / "family_realbulk_validation"
MODELS = ["ridge", "elastic_net"]
# families previously found learnable on breast (where present)
FAMILIES = ["T/NK", "Myeloid", "Endothelial", "Epithelial"]
CONFOUND = {"lib_cv": 0.5, "mult_noise_cv": 0.3, "dropout": 0.2}


def main(argv=None) -> int:
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument("--run-real-data", action="store_true")
    ap.add_argument("--seeds", type=int, nargs="+", default=[0, 1, 2, 3, 4])
    args = ap.parse_args(argv)
    if not H.real_data_enabled(args.run_real_data):
        print(f"Refusing without --run-real-data (or {H.REAL_DATA_ENV}=1).", file=sys.stderr)
        return 2

    from tissueresolve.experimental.conditional_family_estimator import (
        fit_conditional_family_estimator, DEFAULT_GATES)

    D = _load_dataset("breast")
    adata, mapping, ctc, dc = D["adata"], D["mapping"], D["celltype_col"], D["donor_col"]
    adata.obs = adata.obs.copy()
    adata.obs["__family"] = adata.obs[ctc].astype(str).map(lambda s: str(mapping.get(s, s)))
    if "assay" not in adata.obs.columns:
        print("No 'assay' column — cross-platform gate not available.", file=sys.stderr)
        return 3
    assay = adata.obs["assay"].astype(str)
    donor = adata.obs[dc].astype(str)
    # platform -> donor sets (donors are platform-specific in this atlas)
    plats = assay.value_counts().index.tolist()[:2]
    pa, pb = plats[0], plats[1]
    don_a = sorted(donor[assay == pa].unique())
    don_b = sorted(donor[assay == pb].unique())
    print(f"platforms: {pa} ({len(don_a)} donors) train -> {pb} ({len(don_b)} donors) test")

    fams = [f for f in FAMILIES if (adata.obs["__family"] == f).sum() >= 200
            and adata.obs.loc[adata.obs["__family"] == f, ctc].nunique() >= 3]
    rows = []
    for fam in fams:
        for model in MODELS:
            for seed in args.seeds:
                # Gate A: cross-platform (train v3 donors, val v2 donors)
                for gate, kw in [("A_cross_platform",
                                  dict(train_donors=don_a, val_donors=don_b)),
                                 ("B_confounders",
                                  dict(val_confound=CONFOUND))]:
                    try:
                        m = fit_conditional_family_estimator(
                            adata, "__family", ctc, dc, fam, feature_mode="combined",
                            model=model, n_train_mixtures=800, n_val_mixtures=250,
                            random_state=seed, **kw)
                        v = m.validation
                        rows.append({"gate": gate, "family": fam, "model": model, "seed": seed,
                                     "n_states": v["n_states"], "rmse_baseline": v["rmse_baseline"],
                                     "rmse_model": v["rmse_model"], "delta_rmse_pct": v["delta_rmse_pct"],
                                     "pearson": v["pearson"], "rare_fpr_model": v.get("rare_fpr_model"),
                                     "rare_fpr_baseline": v.get("rare_fpr_baseline"),
                                     "calibration_error": v["calibration_error"], "decision_seed": m.decision})
                        print(f"  {gate:<18} {fam:<12} {model:<12} s{seed} "
                              f"rmse {v['rmse_baseline']:.3f}->{v['rmse_model']:.3f} "
                              f"(Δ{v['delta_rmse_pct']:+.0f}%) r={v['pearson']:.2f} {m.decision}")
                    except Exception as exc:  # noqa: BLE001
                        rows.append({"gate": gate, "family": fam, "model": model, "seed": seed,
                                     "error": str(exc)})
                        print(f"  {gate} {fam} {model} s{seed} FAILED: {exc}")

    df = pd.DataFrame(rows)
    OUT.mkdir(parents=True, exist_ok=True)
    df.to_csv(OUT / "per_seed_metrics.tsv", sep="\t", index=False)

    g = DEFAULT_GATES
    summ = []
    ok = df[df.get("delta_rmse_pct").notna()] if "delta_rmse_pct" in df.columns else df.iloc[0:0]
    for (gate, fam, model), gg in ok.groupby(["gate", "family", "model"]):
        delta = gg["delta_rmse_pct"].to_numpy()
        n_improve = int((delta >= g["min_delta_rmse_pct"]).sum())
        fpr_inc = float((gg["rare_fpr_model"].fillna(0) - gg["rare_fpr_baseline"].fillna(0)).mean())
        cal = float(gg["calibration_error"].mean())
        rmse_cv = float(gg["rmse_model"].std() / max(gg["rmse_model"].mean(), 1e-9))
        if cal > g["max_calibration_error"] or fpr_inc > g["max_rare_fpr_increase"]:
            dec = "diagnostic_only"
        elif rmse_cv > g["max_seed_instability"]:
            dec = "unstable"
        elif n_improve >= 4 and delta.mean() >= g["min_delta_rmse_pct"]:
            dec = "learnable"
        elif delta.mean() > g["partial_delta_rmse_pct"]:
            dec = "partially_learnable"
        else:
            dec = "not_learnable"
        summ.append({"gate": gate, "family": fam, "model": model,
                     "delta_rmse_pct": round(delta.mean(), 1), "n_seeds_improved": n_improve,
                     "pearson": round(gg["pearson"].mean(), 3), "rare_fpr_increase": round(fpr_inc, 3),
                     "calibration_error": round(cal, 4), "decision": dec})
    sm = pd.DataFrame(summ)
    sm.to_csv(OUT / "summary.tsv", sep="\t", index=False)
    # best model per (gate, family)
    best = sm.sort_values("delta_rmse_pct", ascending=False).groupby(["gate", "family"]).head(1) \
        if not sm.empty else sm
    best.to_csv(OUT / "best_per_gate_family.tsv", sep="\t", index=False)

    verdict = {"gates": {}, "tcga_real_bulk": "deferred: no fine-state ground truth"}
    if not best.empty:
        for _, r in best.iterrows():
            verdict["gates"][f"{r['gate']}:{r['family']}"] = {
                "best_model": r["model"], "delta_rmse_pct": r["delta_rmse_pct"],
                "decision": r["decision"], "rare_fpr_increase": r["rare_fpr_increase"]}
    (OUT / "verdict.json").write_text(json.dumps(verdict, indent=2), encoding="utf-8")
    (OUT / "manifest.json").write_text(json.dumps({
        "dataset": "breast", "train_platform": pa, "test_platform": pb,
        "families": fams, "models": MODELS, "seeds": args.seeds, "confounders": CONFOUND,
        "note": "real-bulk-realism gates (cross-platform + technical confounders); "
                "TCGA conditional accuracy deferred (no ground truth)"}, indent=2), encoding="utf-8")
    print(f"\nWrote -> {OUT}/")
    if not best.empty:
        print("\n=== best model per gate × family ===")
        print(best.to_string(index=False))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
