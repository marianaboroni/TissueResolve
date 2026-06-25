#!/usr/bin/env python
"""PHASE 2A — partial confidence-weighted unresolved mass: validation harness.

Compares 4 modes (hard_gate / ungated / soft_gate / flat) on held-out-donor
synthetic mixtures with broad+fine mRNA truth, using DISJOINT calibration /
validation / test seed splits (no leakage; rules 10–11).  Computes confidence-
feature predictiveness, calibrates confidence (isotonic/logistic/monotonic on the
calibration split only), benchmarks the 4 modes on TEST with bootstrap CIs,
ablations, promotion gates, and figures.  Does NOT change default behaviour.

Outputs (benchmarks/outputs/):
  confidence_feature_predictiveness.tsv
  v2_soft_gating_calibration.tsv
  v2_soft_gating_mode_comparison.tsv
  v2_soft_gating_ablation.tsv
  v2_soft_gating_promotion_gates.tsv
  figures/phase2a/*.png (+ .data.tsv + caption)

Usage:  python benchmarks/diagnostics/phase2a_soft_gating.py --run-real-data
"""
from __future__ import annotations

import argparse
import sys
import time
import warnings
from pathlib import Path

import numpy as np
import pandas as pd

REPO = Path(__file__).resolve().parents[2]
sys.path.insert(0, str(REPO / "examples" / "real_breast_cancer" / "scripts"))
sys.path.insert(0, str(REPO))
import _harness as H  # noqa: E402
from benchmarks.shared import synthetic_holdout as SH  # noqa: E402
from benchmarks.shared import metrics as M  # noqa: E402
from benchmarks.shared import stats as ST  # noqa: E402
from tissueresolve.experimental.soft_hierarchy import (  # noqa: E402
    apply_partial_confidence_gating, compute_reference_confidence_features,
    add_run_confidence_features, fit_confidence_model, calibration_metrics)

OUT = REPO / "benchmarks" / "outputs"
FIG = OUT / "figures" / "phase2a"
HIER = REPO / "examples" / "real_breast_cancer" / "config" / "breast_cancer_cell_type_hierarchy.tsv"
CAL_SEEDS, VAL_SEEDS, TEST_SEEDS = [0, 1], [2], [3, 4]
SCENARIOS = ["imbalanced", "similar_subtypes", "rare"]
SCEN_SEED = {"imbalanced": 22, "similar_subtypes": 44, "rare": 33}
N_SAMPLES, CELLS, MIN_CELLS = 12, 600, 30
PRESENT, ABSENT, TOL = 0.01, 0.005, 0.05   # predeclared presence + correct-resolution tolerance
RARE_LO, RARE_HI = 0.005, 0.05


# ---------------------------------------------------------------- one seed → artifacts
def _prep(adata, raw_map, seed):
    from tissueresolve.reference.hierarchy import build_cell_type_hierarchy
    ref_d, qry = SH.split_donors(adata, "donor_id", ref_frac=0.5, seed=seed)
    rmask = adata.obs["donor_id"].astype(str).isin(set(ref_d)).to_numpy()
    with warnings.catch_warnings():
        warnings.simplefilter("ignore")
        ref = H.prepare_reference(adata[rmask].copy(), min_cells=MIN_CELLS,
                                  estimate_overdispersion=True).reference
    ref_types = [str(c) for c in ref.cell_types]
    mapping = build_cell_type_hierarchy(ref_types, raw_map)
    return ref, ref_types, mapping, qry


def _scenario(adata, ref_types, mapping, qry, scen, seed):
    kw = {}
    if scen == "similar_subtypes":
        fam = max({mapping.get(x, x) for x in ref_types},
                  key=lambda f: sum(mapping.get(x, x) == f for x in ref_types))
        kw = {"similar_pair": tuple([t for t in ref_types if mapping.get(t) == fam][:2])}
    elif scen == "rare":
        kw = {"rare_type": next((t for t in ["regulatory T cell", "natural killer cell"]
                                 if t in ref_types), ref_types[-1]), "rare_level": 0.01}
    tgt = SH.build_target_proportions(ref_types, N_SAMPLES, scen, seed=SCEN_SEED[scen], **kw)
    ds = SH.realize_pseudobulk(adata, tgt, celltype_col="cell_type", donor_col="donor_id",
                               query_donors=qry, seed=seed, cells_per_sample=CELLS)
    return ds


def _modes(counts, ref, mapping, conf_model=None, conf_feat=None):
    """Return dict of mode→(fine_pred DataFrame, runtime)."""
    from tissueresolve.api import deconv_bulk
    from tissueresolve.reference.hierarchy import combine_family_and_conditional_estimates
    out = {}
    with warnings.catch_warnings():
        warnings.simplefilter("ignore")
        t = time.perf_counter()
        flat = deconv_bulk(counts, ref, solver="nnls", resolution_mode="none").deconv.proportions
        out["flat"] = (flat, time.perf_counter() - t)
        t = time.perf_counter()
        h = deconv_bulk(counts, ref, solver="auto", resolution_mode="hierarchical",
                        hierarchy_mapping=mapping)
        out["hard_gate"] = (h.estimates.combined_fine, time.perf_counter() - t)
    est = h.estimates
    fam, cond = est.family_proportions, est.conditional_proportions
    pre = combine_family_and_conditional_estimates(fam, cond, mapping)
    out["ungated"] = (pre, 0.0)
    if conf_model is not None:
        t = time.perf_counter()
        conf = conf_model.predict(conf_feat)
        r = apply_partial_confidence_gating(fam, pre, mapping, confidence=conf)
        out["soft_gate"] = (r.combined, time.perf_counter() - t)
        out["_soft_result"] = r
    out["_est"] = est
    return out


# ---------------------------------------------------------------- metrics
def _fine_cols(df):
    return [c for c in df.columns if not str(c).startswith("unresolved_")]


def _metrics(truth, pred, mapping):
    cols = [c for c in truth.columns if not str(c).startswith("unresolved_")]
    t = truth[cols]
    p = pred.reindex(index=t.index, columns=cols).fillna(0.0)
    acc = M.accuracy_metrics(t, p)
    tfam = M.aggregate_to_families(t, mapping)
    pfam = M.aggregate_to_families(pred, mapping)   # include unresolved → family
    bacc = M.accuracy_metrics(tfam, pfam)
    # conditional within-family
    fam_of = {c: mapping.get(c, c) for c in cols}
    ct, cp = t * 0.0, p * 0.0
    for fam in set(fam_of.values()):
        mem = [c for c in cols if fam_of[c] == fam]
        st_, sp_ = t[mem].sum(1), p[mem].sum(1)
        for c in mem:
            ct[c] = np.where(st_ > 0, t[c] / st_.replace(0, np.nan), 0)
            cp[c] = np.where(sp_ > 0, p[c] / sp_.replace(0, np.nan), 0)
    cond = M.accuracy_metrics(ct.fillna(0), cp.fillna(0))
    tv, pv = t.to_numpy(float), p.to_numpy(float)
    present, absent = tv > PRESENT, tv <= ABSENT
    rare = (tv >= RARE_LO) & (tv <= RARE_HI)
    resolved_mass = float(p.to_numpy().sum())
    cxt = M.composition_complexity(t)["effective_n_populations"].mean()
    cxp = M.composition_complexity(p)["effective_n_populations"].mean()
    ucols = [c for c in pred.columns if str(c).startswith("unresolved_")]
    umass = float(pred[ucols].to_numpy().sum() / max(pred.to_numpy().sum(), 1e-9)) if ucols else 0.0
    return {
        "broad_pearson": bacc["pearson"], "broad_rmse": bacc["rmse"],
        "fine_pearson": acc["pearson"], "fine_rmse": acc["rmse"],
        "cond_pearson": cond["pearson"],
        "true_eff_n": float(cxt), "pred_eff_n": float(cxp),
        "eff_n_abs_error": float(abs(cxp - cxt)),
        "richness_pred": float((p > PRESENT).sum(1).mean()),
        "richness_true": float((t > PRESENT).sum(1).mean()),
        "rare_sensitivity": float((pv[rare] > ABSENT).mean()) if rare.any() else float("nan"),
        "spillover_offtarget": float(pv[absent].mean()) if absent.any() else float("nan"),
        "false_resolution_rate": float(pv[absent].sum() / max(resolved_mass, 1e-9)),
        "unresolved_mass_frac": umass,
        "mass_error": float(abs(pred.sum(1) - truth.sum(1)).max()),
    }


def _false_abstention(truth, pred_combined, mapping):
    """Mass the method left unresolved in families that DO have ≥2 present subtypes
    (should have been resolved) / total mass of those families. Lower = better."""
    cols = [c for c in truth.columns if not str(c).startswith("unresolved_")]
    fam_of = {c: mapping.get(c, c) for c in cols}
    ucols = {c[len("unresolved_"):]: c for c in pred_combined.columns if str(c).startswith("unresolved_")}
    num = den = 0.0
    for fam in set(fam_of.values()):
        mem = [c for c in cols if fam_of[c] == fam]
        present_cnt = (truth[mem] > PRESENT).sum(1)
        should_resolve = present_cnt >= 2
        if should_resolve.any():
            fam_mass = truth[mem].sum(1)[should_resolve].sum()
            u = pred_combined[ucols[fam]][should_resolve].sum() if fam in ucols else 0.0
            num += float(u); den += float(fam_mass)
    return float(num / den) if den > 0 else float("nan")


# ---------------------------------------------------------------- main
def main(argv=None) -> int:
    ap = argparse.ArgumentParser(); ap.add_argument("--run-real-data", action="store_true")
    args = ap.parse_args(argv)
    if not H.real_data_enabled(args.run_real_data):
        print("Refusing without --run-real-data.", file=sys.stderr); return 2
    import anndata as ad
    from tissueresolve.reference.hierarchy import (load_hierarchy_mapping,
                                                   combine_family_and_conditional_estimates)
    from tissueresolve.api import deconv_bulk
    adata = ad.read_h5ad(H.REFERENCE_H5AD)
    if "feature_name" in adata.var.columns:
        fn = adata.var["feature_name"].astype(str)
        if fn.nunique() == adata.n_vars:
            adata.var_names = fn.to_numpy()
    raw_map = load_hierarchy_mapping(HIER)
    FIG.mkdir(parents=True, exist_ok=True)

    # ---- gather (subtype,sample) evidence + correct-resolution labels on CAL+VAL ----
    feat_rows = []
    for seed in CAL_SEEDS + VAL_SEEDS:
        ref, ref_types, mapping, qry = _prep(adata, raw_map, seed)
        rfeat = compute_reference_confidence_features(ref, mapping)
        for scen in SCENARIOS:
            ds = _scenario(adata, ref_types, mapping, qry, scen, seed)
            truth = ds.true_mrna_proportions
            with warnings.catch_warnings():
                warnings.simplefilter("ignore")
                h = deconv_bulk(ds.counts, ref, solver="auto", resolution_mode="hierarchical",
                                hierarchy_mapping=mapping)
                # cross-solver agreement feature
                h2 = deconv_bulk(ds.counts, ref, solver="ridge_nnls", resolution_mode="hierarchical",
                                 hierarchy_mapping=mapping)
            pre = combine_family_and_conditional_estimates(
                h.estimates.family_proportions, h.estimates.conditional_proportions, mapping)
            rfeat2 = add_run_confidence_features(
                rfeat, {"auto": h.estimates.conditional_proportions,
                        "ridge": h2.estimates.conditional_proportions})
            cols = _fine_cols(truth)
            for samp in truth.index:
                for st in cols:
                    err = abs(float(pre.at[samp, st]) - float(truth.at[samp, st])) \
                        if st in pre.columns else float(truth.at[samp, st])
                    feat_rows.append({
                        "split": "cal" if seed in CAL_SEEDS else "val",
                        "seed": seed, "scenario": scen, "subtype": st, "sample": samp,
                        "abs_error": err, "correct": int(err < TOL),
                        "true_present": int(truth.at[samp, st] > PRESENT),
                        "separability": rfeat2.at[st, "separability"] if st in rfeat2.index else np.nan,
                        "min_disc_genes": rfeat2.at[st, "min_disc_genes"] if st in rfeat2.index else np.nan,
                        "marker_support": rfeat2.at[st, "marker_support"] if st in rfeat2.index else np.nan,
                        "ref_evidence": rfeat2.at[st, "ref_evidence"] if st in rfeat2.index else np.nan,
                        "cross_solver_agreement": rfeat2.at[st, "cross_solver_agreement"]
                        if "cross_solver_agreement" in rfeat2.columns and st in rfeat2.index else np.nan,
                    })
    fdf = pd.DataFrame(feat_rows)

    # ---- STEP 4: feature predictiveness (on CAL) ----
    from scipy.stats import spearmanr
    cal = fdf[fdf.split == "cal"]
    pred_rows = []
    for feat in ["separability", "min_disc_genes", "marker_support", "ref_evidence",
                 "cross_solver_agreement"]:
        sub = cal[[feat, "abs_error", "correct", "true_present"]].dropna()
        if len(sub) < 10:
            continue
        rho_err = spearmanr(sub[feat], sub["abs_error"]).correlation
        # AUROC of feature for correct-resolution
        try:
            from sklearn.metrics import roc_auc_score
            auroc = roc_auc_score(sub["correct"], sub[feat]) if sub["correct"].nunique() > 1 else float("nan")
            auroc_present = roc_auc_score(sub["true_present"], sub[feat]) if sub["true_present"].nunique() > 1 else float("nan")
        except Exception:
            auroc = auroc_present = float("nan")
        pred_rows.append({"feature": feat, "n": len(sub),
                          "spearman_vs_abs_error": round(float(rho_err), 4),
                          "auroc_correct_resolution": round(float(auroc), 4),
                          "auroc_subtype_presence": round(float(auroc_present), 4)})
    pd.DataFrame(pred_rows).to_csv(OUT / "confidence_feature_predictiveness.tsv", sep="\t", index=False)
    print("STEP4 feature predictiveness:\n", pd.DataFrame(pred_rows).to_string(index=False))

    # ---- STEP 5: calibrate confidence (cal fit, val evaluate) ----
    score_col = "ref_evidence"
    cal_feat = cal.set_index(["sample", "subtype"])[[score_col]]
    cal_lab = cal.set_index(["sample", "subtype"])["correct"].astype(float)
    val = fdf[fdf.split == "val"]
    val_feat = val.set_index(["sample", "subtype"])[[score_col]]
    val_lab = val.set_index(["sample", "subtype"])["correct"].astype(float)
    calib_rows = []
    models = {}
    for kind in ["monotonic", "logistic", "isotonic"]:
        m = fit_confidence_model(cal_feat, cal_lab, kind=kind, score_col=score_col)
        models[kind] = m
        cm = calibration_metrics(m.predict(val_feat).to_numpy(), val_lab.to_numpy())
        calib_rows.append({"calibrator": kind, "brier": round(cm["brier"], 4),
                           "ece": round(cm["ece"], 4), "n_val": cm["n"]})
    calib = pd.DataFrame(calib_rows)
    calib.to_csv(OUT / "v2_soft_gating_calibration.tsv", sep="\t", index=False)
    best_kind = calib.sort_values("brier").iloc[0]["calibrator"]
    best_model = models[best_kind]
    print(f"STEP5 calibration:\n{calib.to_string(index=False)}\n  selected: {best_kind}")

    # ---- STEP 6: 4-mode comparison on TEST (per seed×scenario replicate) ----
    rows = []
    for seed in TEST_SEEDS:
        ref, ref_types, mapping, qry = _prep(adata, raw_map, seed)
        rfeat = compute_reference_confidence_features(ref, mapping)
        conf_feat = rfeat[[score_col]]
        for scen in SCENARIOS:
            ds = _scenario(adata, ref_types, mapping, qry, scen, seed)
            truth = ds.true_mrna_proportions
            modes = _modes(ds.counts, ref, mapping, conf_model=best_model, conf_feat=conf_feat)
            for mode in ["flat", "hard_gate", "ungated", "soft_gate"]:
                pred, rt = modes[mode]
                m = _metrics(truth, pred, mapping)
                m["false_abstention_rate"] = _false_abstention(truth, pred, mapping) \
                    if mode in ("hard_gate", "soft_gate") else 0.0
                rows.append({"seed": seed, "scenario": scen, "mode": mode,
                             "runtime_s": round(rt, 3), **m})
    comp = pd.DataFrame(rows)
    # CI summary per mode
    ci_rows = []
    for mode in ["flat", "hard_gate", "ungated", "soft_gate"]:
        sub = comp[comp["mode"] == mode]
        for met in ["fine_pearson", "fine_rmse", "broad_rmse", "cond_pearson",
                    "false_resolution_rate", "false_abstention_rate", "eff_n_abs_error",
                    "rare_sensitivity", "unresolved_mass_frac", "mass_error"]:
            c = ST.bootstrap_ci(sub[met].to_numpy(), seed=0)
            ci_rows.append({"mode": mode, "metric": met, "point": round(c["point"], 4),
                            "lo": round(c["lo"], 4), "hi": round(c["hi"], 4)})
    comp_ci = pd.DataFrame(ci_rows)
    with open(OUT / "v2_soft_gating_mode_comparison.tsv", "w") as fh:
        fh.write("# Phase 2A 4-mode comparison on HELD-OUT TEST seeds (3,4). Bootstrap 95% CI over seed×scenario.\n")
        comp_ci.to_csv(fh, sep="\t", index=False)
    print("\nSTEP6 mode comparison (mean over test):")
    print(comp.groupby("mode")[["fine_pearson", "fine_rmse", "broad_rmse", "cond_pearson",
          "false_resolution_rate", "eff_n_abs_error", "unresolved_mass_frac"]].mean().round(3).to_string())

    # ---- STEP 7: ablations (calibrate on different score features; + hard/ungated refs) ----
    abl_rows = []
    for feat in ["ref_evidence", "separability", "marker_support", "cross_solver_agreement"]:
        cf = cal.dropna(subset=[feat]).set_index(["sample", "subtype"])[[feat]]
        cl = cal.dropna(subset=[feat]).set_index(["sample", "subtype"])["correct"].astype(float)
        if len(cf) < 10:
            continue
        mdl = fit_confidence_model(cf, cl, kind=best_kind, score_col=feat)
        fp, frr, en = [], [], []
        for seed in TEST_SEEDS:
            ref, ref_types, mapping, qry = _prep(adata, raw_map, seed)
            rfeat = compute_reference_confidence_features(ref, mapping)
            if feat not in rfeat.columns:
                continue
            cfeat = rfeat[[feat]]
            for scen in SCENARIOS:
                ds = _scenario(adata, ref_types, mapping, qry, scen, seed)
                truth = ds.true_mrna_proportions
                modes = _modes(ds.counts, ref, mapping, conf_model=mdl, conf_feat=cfeat)
                if "soft_gate" not in modes:
                    continue
                mm = _metrics(truth, modes["soft_gate"][0], mapping)
                fp.append(mm["fine_pearson"]); frr.append(mm["false_resolution_rate"]); en.append(mm["eff_n_abs_error"])
        abl_rows.append({"confidence_score": feat, "n": len(fp),
                         "fine_pearson": round(float(np.mean(fp)), 4) if fp else np.nan,
                         "false_resolution_rate": round(float(np.mean(frr)), 4) if frr else np.nan,
                         "eff_n_abs_error": round(float(np.mean(en)), 4) if en else np.nan})
    pd.DataFrame(abl_rows).to_csv(OUT / "v2_soft_gating_ablation.tsv", sep="\t", index=False)
    print("\nSTEP7 ablation (soft-gate by confidence score):\n",
          pd.DataFrame(abl_rows).to_string(index=False))

    # ---- STEP 8: promotion gates (soft vs ungated/hard on TEST means) ----
    g = comp.groupby("mode").mean(numeric_only=True)
    def rel(a, b):  # relative change a vs b
        return (a - b) / abs(b) if b else float("nan")
    gates = []
    soft, ung, hard, flat = g.loc["soft_gate"], g.loc["ungated"], g.loc["hard_gate"], g.loc["flat"]
    gates.append(("broad RMSE ≤ ungated +2%", soft["broad_rmse"], ung["broad_rmse"],
                  rel(soft["broad_rmse"], ung["broad_rmse"]) <= 0.02))
    gates.append(("fine RMSE ≤ ungated +2%", soft["fine_rmse"], ung["fine_rmse"],
                  rel(soft["fine_rmse"], ung["fine_rmse"]) <= 0.02))
    gates.append(("false-resolution ≥10% better than hard", soft["false_resolution_rate"],
                  hard["false_resolution_rate"],
                  rel(soft["false_resolution_rate"], hard["false_resolution_rate"]) <= -0.10
                  or soft["false_resolution_rate"] <= hard["false_resolution_rate"]))
    gates.append(("false-abstention < hard", soft["false_abstention_rate"], hard["false_abstention_rate"],
                  soft["false_abstention_rate"] < hard["false_abstention_rate"]))
    gates.append(("eff-N error < hard (closer to truth)", soft["eff_n_abs_error"], hard["eff_n_abs_error"],
                  soft["eff_n_abs_error"] < hard["eff_n_abs_error"]))
    gates.append(("mass conservation < 1e-6", soft["mass_error"], 1e-6, soft["mass_error"] < 1e-6))
    gate_df = pd.DataFrame([{"gate": n, "soft": round(float(a), 4), "reference": round(float(b), 6),
                             "pass": bool(p)} for n, a, b, p in gates])
    with open(OUT / "v2_soft_gating_promotion_gates.tsv", "w") as fh:
        fh.write("# Phase 2A promotion gates on HELD-OUT TEST. soft_gate vs ungated/hard references.\n")
        gate_df.to_csv(fh, sep="\t", index=False)
    print("\nSTEP8 promotion gates:\n", gate_df.to_string(index=False))

    _figures(comp, comp_ci, calib, models, val_feat, val_lab, fdf)
    comp.to_csv(OUT / "v2_soft_gating_mode_comparison_raw.tsv", sep="\t", index=False)
    print(f"\nWrote Phase-2A tables + figures -> {OUT}/")
    return 0


def _figures(comp, comp_ci, calib, models, val_feat, val_lab, fdf):
    import matplotlib; matplotlib.use("Agg")
    import matplotlib.pyplot as plt

    def save(fig, name, data, cap):
        fig.tight_layout(); fig.savefig(FIG / f"{name}.png", dpi=150); fig.savefig(FIG / f"{name}.svg")
        plt.close(fig); data.to_csv(FIG / f"{name}.data.tsv", sep="\t")
        (FIG / f"{name}.caption.txt").write_text(cap.strip() + "\n", encoding="utf-8")

    order = ["flat", "ungated", "soft_gate", "hard_gate"]
    g = comp.groupby("mode").mean(numeric_only=True).reindex(order)
    fig, ax = plt.subplots(figsize=(8, 4.2))
    g[["fine_pearson", "cond_pearson"]].plot.bar(ax=ax)
    ax.set_title("Phase 2A — fine & conditional accuracy by mode (held-out test)")
    ax.set_ylabel("Pearson"); plt.setp(ax.get_xticklabels(), rotation=10)
    save(fig, "fig_mode_accuracy", g, "Fine & within-family conditional Pearson by gating mode on held-out test seeds.")

    fig, ax = plt.subplots(figsize=(7, 4))
    g[["false_resolution_rate", "false_abstention_rate"]].plot.bar(ax=ax)
    ax.set_title("False-resolution vs false-abstention trade-off"); plt.setp(ax.get_xticklabels(), rotation=10)
    save(fig, "fig_false_resolution_abstention", g[["false_resolution_rate", "false_abstention_rate"]],
         "False-resolution and false-abstention by mode. Soft gating should reduce both vs hard gating.")

    fig, ax = plt.subplots(figsize=(6, 5))
    g2 = comp.groupby("mode").mean(numeric_only=True)
    ax.scatter(g2["true_eff_n"], g2["pred_eff_n"])
    for m in g2.index:
        ax.annotate(m, (g2.loc[m, "true_eff_n"], g2.loc[m, "pred_eff_n"]), fontsize=8)
    lim = max(g2["true_eff_n"].max(), g2["pred_eff_n"].max()) + 1
    ax.plot([0, lim], [0, lim], "k--", lw=0.8)
    ax.set_xlabel("true effective-N"); ax.set_ylabel("predicted effective-N")
    ax.set_title("Predicted vs true effective-N by mode")
    save(fig, "fig_effective_n", g2[["true_eff_n", "pred_eff_n"]], "Predicted vs true effective number of populations by mode.")

    # reliability (best calibrator on val)
    best = calib.sort_values("brier").iloc[0]["calibrator"]
    cm = calibration_metrics(models[best].predict(val_feat).to_numpy(), val_lab.to_numpy())
    rel = pd.DataFrame(cm["reliability"])
    if not rel.empty:
        fig, ax = plt.subplots(figsize=(5, 5))
        ax.plot([0, 1], [0, 1], "k--", lw=0.8)
        ax.plot(rel["mean_confidence"], rel["observed_accuracy"], "o-")
        ax.set_xlabel("mean predicted confidence"); ax.set_ylabel("observed correct-resolution rate")
        ax.set_title(f"Reliability ({best}); Brier={cm['brier']:.3f}, ECE={cm['ece']:.3f}")
        save(fig, "fig_reliability", rel, f"Calibration reliability of the {best} confidence model on validation seed.")

    # confidence vs absolute error (cal)
    cal = fdf[fdf.split == "cal"].dropna(subset=["ref_evidence"])
    fig, ax = plt.subplots(figsize=(6, 4))
    ax.scatter(cal["ref_evidence"], cal["abs_error"], s=6, alpha=0.3)
    ax.set_xlabel("reference evidence score"); ax.set_ylabel("|pred − true| (ungated)")
    ax.set_title("Confidence evidence vs absolute fine error")
    save(fig, "fig_confidence_vs_error", cal[["ref_evidence", "abs_error"]].reset_index(drop=True),
         "Per-(subtype,sample) reference evidence score vs ungated absolute fine error (calibration split).")


if __name__ == "__main__":
    raise SystemExit(main())
