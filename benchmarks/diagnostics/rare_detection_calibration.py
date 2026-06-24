#!/usr/bin/env python
"""P3 — calibrate the rare-detection layer on breast: reduce the GLM's rare-subtype
false-positive rate without destroying recall.

Uses the exported breast donor-disjoint scenarios. Runs the Poisson GLM, computes
marker-support + detection scores for the rare type across all samples, fits the
calibration on a held-out split of scenarios (never tuned on the evaluation split),
emits the precision–recall curve, picks an operating point, and reports rare
precision / recall / FPR BEFORE vs AFTER gating.

Usage:
  PYTHONPATH=src:. python benchmarks/diagnostics/rare_detection_calibration.py --run-real-data
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

EXT = REPO / "benchmarks" / "outputs" / "holdout_bulk" / "external_inputs"
OUT = REPO / "benchmarks" / "outputs" / "nb_bulk_solver" / "rare_calibration"
HIER = (REPO / "examples" / "real_breast_cancer" / "config" /
        "breast_cancer_cell_type_hierarchy.tsv")
PRESENT_LEVEL = 0.02         # truth mRNA proportion above which the rare type is "present"
CALIB_SCENARIOS = {"imbalanced", "missing_population"}   # fit here; evaluate on the rest


def _rates(pred_rare, truth_rare, present):
    tp = int((pred_rare & present).sum()); fp = int((pred_rare & ~present).sum())
    fn = int((~pred_rare & present).sum()); tn = int((~pred_rare & ~present).sum())
    prec = tp / (tp + fp) if (tp + fp) else float("nan")
    rec = tp / (tp + fn) if (tp + fn) else float("nan")
    fpr = fp / (fp + tn) if (fp + tn) else float("nan")
    return prec, rec, fpr


def main(argv=None) -> int:
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument("--run-real-data", action="store_true")
    ap.add_argument("--rare", default=None,
                    help="Override the rare cell type (default: manifest rare_type).")
    args = ap.parse_args(argv)
    if not H.real_data_enabled(args.run_real_data):
        print(f"Refusing without --run-real-data (or {H.REAL_DATA_ENV}=1).", file=sys.stderr)
        return 2

    import anndata as ad
    from tissueresolve.reference.hierarchy import load_hierarchy_mapping, build_cell_type_hierarchy
    from tissueresolve.config import TissueResolveConfig
    import tissueresolve as tr
    from tissueresolve.experimental.rare_detection import (
        compute_marker_support, rare_detection_probability, calibrate_detection,
        precision_recall_curve, apply_rare_detection_gate)

    manifest = json.loads((EXT / "export_manifest.json").read_text())
    rare = args.rare or manifest["rare_type"]
    OUT.mkdir(parents=True, exist_ok=True)

    adata = ad.read_h5ad(H.REFERENCE_H5AD)
    if "feature_name" in adata.var.columns:
        fn = adata.var["feature_name"].astype(str)
        if fn.nunique() == adata.n_vars:
            adata.var_names = fn.to_numpy()
    rmask = adata.obs["donor_id"].astype(str).isin(set(manifest["reference_donors"])).to_numpy()
    import warnings as _w
    with _w.catch_warnings():
        _w.simplefilter("ignore")
        ref = H.prepare_reference(adata[rmask].copy(), min_cells=30,
                                  estimate_overdispersion=True).reference
    cts = [str(c) for c in ref.cell_types]
    mapping = build_cell_type_hierarchy(cts, load_hierarchy_mapping(HIER))

    # collect rare predicted mass, marker support, truth-present label across scenarios
    rows = []
    cfg = TissueResolveConfig(); cfg.bulk_solver.method = "poisson_glm_experimental"
    for scen_dir in sorted([d for d in EXT.iterdir() if (d / "truth_mrna.tsv").exists()]):
        scen = scen_dir.name
        truth = pd.read_csv(scen_dir / "truth_mrna.tsv", sep="\t", index_col=0)
        bulk = pd.read_csv(scen_dir / "bulk_counts_genes_by_samples.tsv", sep="\t", index_col=0)
        if rare not in truth.columns:
            continue
        with _w.catch_warnings():
            _w.simplefilter("ignore")
            pred = tr.deconv_bulk(bulk, ref, config=cfg, resolution_mode="none").deconv.proportions
        sup = compute_marker_support(bulk, ref, [rare], n_markers=25)
        score = rare_detection_probability(pred[[rare]].rename(columns={rare: rare}),
                                           sup, [rare])
        for s in truth.index:
            rows.append({"scenario": scen, "sample": s,
                         "pred_rare": float(pred.loc[s, rare]),
                         "truth_rare": float(truth.loc[s, rare]),
                         "present": bool(truth.loc[s, rare] > PRESENT_LEVEL),
                         "marker_support": float(sup.loc[s, rare]),
                         "score": float(score.loc[s, rare])})
    df = pd.DataFrame(rows)
    df.to_csv(OUT / "rare_sample_table.tsv", sep="\t", index=False)

    calib = df[df.scenario.isin(CALIB_SCENARIOS)]
    evalset = df[~df.scenario.isin(CALIB_SCENARIOS)]
    cal = calibrate_detection(calib["score"].to_numpy(), calib["present"].to_numpy(),
                              kind="isotonic")
    evalset = evalset.copy()
    evalset["prob"] = cal(evalset["score"].to_numpy())

    # PR curve on the evaluation split
    pr = precision_recall_curve(evalset["prob"].to_numpy(), evalset["present"].to_numpy())
    pr.to_csv(OUT / "rare_precision_recall_curve.tsv", sep="\t", index=False)

    # average precision (PR-AUC proxy) on eval — quantifies identifiability of the rare type
    pr_eval = pr.dropna(subset=["precision", "recall"]).sort_values("recall")
    ap = float(np.trapz(pr_eval["precision"], pr_eval["recall"])) if len(pr_eval) > 1 else float("nan")

    # BEFORE: raw GLM call = predicted mass > PRESENT_LEVEL (the un-gated behaviour)
    present = evalset["present"].to_numpy()
    raw_call = evalset["pred_rare"].to_numpy() > PRESENT_LEVEL
    p0, r0, f0 = _rates(raw_call, None, present)

    # operating point chosen on the CALIBRATION split: max recall s.t. FPR <= target
    pr_cal = precision_recall_curve(cal(calib["score"].to_numpy()), calib["present"].to_numpy())
    feas_cal = pr_cal[(pr_cal["fpr"] <= 0.20) & pr_cal["recall"].notna()]
    thr = float(feas_cal.sort_values("recall", ascending=False)["threshold"].iloc[0]) \
        if not feas_cal.empty else 0.5
    gated = (evalset["prob"].to_numpy() >= thr) & (evalset["marker_support"].to_numpy() >= 0.3)
    p1, r1, f1 = _rates(gated, None, present)

    # best ACHIEVABLE on eval at FPR <= target (upper bound the gate could reach)
    feas_eval = pr[(pr["fpr"] <= 0.20) & pr["recall"].notna()]
    if not feas_eval.empty:
        best = feas_eval.sort_values("recall", ascending=False).iloc[0]
        p2, r2, f2, thr2 = best["precision"], best["recall"], best["fpr"], best["threshold"]
    else:
        p2 = r2 = f2 = thr2 = float("nan")

    summary = pd.DataFrame([
        {"stage": "before (raw GLM mass>0.02)", "precision": p0, "recall": r0, "fpr": f0},
        {"stage": f"after (calib-selected gate prob>={thr:.2f} & support>=0.3)",
         "precision": p1, "recall": r1, "fpr": f1},
        {"stage": f"best achievable on eval (prob>={thr2:.2f}, FPR<=0.20)",
         "precision": p2, "recall": r2, "fpr": f2},
    ])
    print(f"rare type: {rare} | eval average-precision (PR-AUC) = {ap:.3f} "
          f"(base rate {present.mean():.3f})")
    summary.to_csv(OUT / "rare_calibration_summary.tsv", sep="\t", index=False)
    (OUT / "manifest.json").write_text(json.dumps({
        "rare_type": rare, "present_level": PRESENT_LEVEL,
        "calibration_scenarios": sorted(CALIB_SCENARIOS),
        "evaluation_scenarios": sorted(set(df.scenario) - CALIB_SCENARIOS),
        "operating_threshold": thr, "min_marker_support": 0.3,
        "fpr_target": 0.20, "n_eval_samples": int(len(evalset)),
        "eval_average_precision": ap, "eval_base_rate": float(present.mean()),
    }, indent=2), encoding="utf-8")
    print(f"Wrote -> {OUT}/")
    print(summary.round(3).to_string(index=False))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
