#!/usr/bin/env python
"""Evaluate the high-granularity strategy (ReferenceCalibration + CandidateRestriction
+ ConsensusStability) on two tissues (breast + HLCA lung), applied BEFORE the
validated partial confidence-weighted soft gate.

Order: π_f → q_raw → calibrate → restrict candidates → multi-panel deconv →
consensus → stability/candidate-aware confidence → soft gate (once) → unresolved.

Design (donor- AND seed-disjoint; rules 10–11): per seed a reference/query donor
split; donor-stability panels use REFERENCE donors only; the soft-gating
confidence model is fit on CALIBRATION seeds only; reference calibration uses only
the query expression (never truth/test donors); all reported metrics from held-out
TEST seeds only.

Ablated methods (all share the SAME soft gate unless noted):
  baseline_softgate     current hierarchical fine + soft gate
  calib_softgate        ReferenceCalibration only
  candidate_softgate    CandidateRestriction only
  consensus_softgate    Multi-panel consensus only
  calib_candidate        ReferenceCalibration + CandidateRestriction
  candidate_consensus    CandidateRestriction + ConsensusStability
  full_softgate          Calibration + CandidateRestriction + Consensus
  flat_nnls / hard_gate / ungated / refiner_negbaseline   (reference baselines)

Outputs (git-ignored benchmarks/outputs/): high_granularity_<tissue>_modes.tsv,
high_granularity_<tissue>_perfamily.tsv, high_granularity_gates.tsv.
Nothing is committed.

Usage:  python high_granularity_benchmark.py --run-real-data --tissue both
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
from tissueresolve.experimental.soft_hierarchy import (  # noqa: E402
    apply_partial_confidence_gating, compute_reference_confidence_features, fit_confidence_model)
from tissueresolve.experimental.soft_hierarchy.fine_refiner import build_donor_residuals  # noqa: E402
from tissueresolve.experimental.soft_hierarchy.high_granularity import (  # noqa: E402
    HighGranularityStrategy, HighGranularityConfig)

OUT = REPO / "benchmarks" / "outputs"
CAL, TEST = [0, 1], [2, 3]
SCENARIOS = ["balanced", "imbalanced", "rare", "similar_subtypes"]
N_SAMPLES, CELLS, MIN_CELLS = 10, 400, 20
PRESENT, ABSENT, TOL = 0.01, 0.005, 0.05

_BREAST_SUPPLEMENT = {
    "IgG plasma cell": "B/Plasma", "naive B cell": "B/Plasma",
    "unswitched memory B cell": "B/Plasma",
    "activated CD8-positive, alpha-beta T cell": "T/NK", "gamma-delta T cell": "T/NK",
    "mast cell": "Myeloid", "myeloid dendritic cell": "Myeloid",
    "neutrophil": "Myeloid", "plasmacytoid dendritic cell": "Myeloid",
}
TISSUES = {
    "breast": {"h5ad": REPO / "examples/real_breast_cancer/data/reference/breast_cancer_sc_reference.h5ad",
               "hmap": REPO / "examples/real_breast_cancer/config/breast_cancer_cell_type_hierarchy.tsv",
               "fine_col": "cell_type", "donor_col": "donor_id", "feature_name": True,
               "supplement": _BREAST_SUPPLEMENT},
    "lung": {"h5ad": REPO / "examples/second_tissue_lung/data/derived/hlca_subset.h5ad",
             "hmap": REPO / "examples/second_tissue_lung/data/derived/hlca_hierarchy.tsv",
             "fine_col": "cell_type_fine", "donor_col": "donor_id", "feature_name": True},
}

# refiner modes for the high-granularity strategy
MODES = {
    "calib_softgate": dict(mode="reference_calibrated_consensus", calibration="scale",
                           candidate_rule="none", panels=("current",), consensus="mean"),
    "candidate_softgate": dict(mode="candidate_consensus", calibration="none",
                               candidate_rule="hybrid", panels=("current",), consensus="mean"),
    "consensus_softgate": dict(mode="candidate_consensus", calibration="none",
                               candidate_rule="none",
                               panels=("current", "marker", "query_detectable", "donor_stable"),
                               consensus="median"),
    "calib_candidate": dict(mode="reference_calibrated_consensus", calibration="scale",
                            candidate_rule="hybrid", panels=("current",), consensus="mean"),
    "candidate_consensus": dict(mode="candidate_consensus", calibration="none",
                                candidate_rule="hybrid",
                                panels=("current", "marker", "query_detectable", "donor_stable"),
                                consensus="conservative"),
    "full_softgate": dict(mode="reference_calibrated_consensus", calibration="scale",
                          candidate_rule="hybrid",
                          panels=("current", "marker", "query_detectable", "donor_stable"),
                          consensus="conservative"),
}


def _load(cfg):
    import anndata as ad
    adata = ad.read_h5ad(cfg["h5ad"])
    if cfg["feature_name"] and "feature_name" in adata.var.columns:
        adata.var_names = adata.var["feature_name"].astype(str).to_numpy()
        adata.var_names_make_unique()
    hmap = pd.read_csv(cfg["hmap"], sep="\t")
    mapping = dict(zip(hmap["fine_cell_type"].astype(str), hmap["broad_cell_type"].astype(str)))
    mapping.update(cfg.get("supplement", {}))
    return adata, mapping


def _prep(adata, mapping, fine_col, donor_col, seed):
    from tissueresolve.reference.hierarchy import build_cell_type_hierarchy
    ref_d, qry = SH.split_donors(adata, donor_col, ref_frac=0.5, seed=seed)
    rmask = adata.obs[donor_col].astype(str).isin(set(ref_d)).to_numpy()
    with warnings.catch_warnings():
        warnings.simplefilter("ignore")
        ref = H.prepare_reference(adata[rmask].copy(), cell_type_col=fine_col,
                                  min_cells=MIN_CELLS, estimate_overdispersion=False).reference
    ref_types = [str(c) for c in ref.cell_types]
    mp = dict(build_cell_type_hierarchy(ref_types, mapping))
    return ref, ref_types, mp, qry, ref_d


def _scenario(adata, ref_types, mp, qry, fine_col, donor_col, scen, seed):
    kw = {}
    if scen == "similar_subtypes":
        fam = max({mp.get(x, x) for x in ref_types},
                  key=lambda f: sum(mp.get(x, x) == f for x in ref_types))
        kw = {"similar_pair": tuple([t for t in ref_types if mp.get(t) == fam][:2])}
    elif scen == "rare":
        kw = {"rare_type": ref_types[-1], "rare_level": 0.01}
    tgt = SH.build_target_proportions(ref_types, N_SAMPLES, scen, seed=1000 + seed, **kw)
    return SH.realize_pseudobulk(adata, tgt, celltype_col=fine_col, donor_col=donor_col,
                                 query_donors=qry, seed=seed, cells_per_sample=CELLS)


def _conditional(df, mp, cols):
    fam_of = {c: mp.get(c, c) for c in cols}
    out = df[cols] * 0.0
    for fam in set(fam_of.values()):
        mem = [c for c in cols if fam_of[c] == fam]
        s = df[mem].sum(axis=1)
        for c in mem:
            out[c] = np.where(s > 0, df[c] / s.replace(0, np.nan), 0.0)
    return out.fillna(0.0)


def _pairwise_spillover(tcond, pcond, mp, cols):
    fam_of = {c: mp.get(c, c) for c in cols}
    leak, n = 0.0, 0
    for fam in set(fam_of.values()):
        mem = [c for c in cols if fam_of[c] == fam]
        if len(mem) < 2:
            continue
        absent = ~(tcond[mem] > PRESENT)
        leaked = (pcond[mem].to_numpy() * absent.to_numpy()).sum(axis=1)
        active = (tcond[mem] > PRESENT).any(axis=1).to_numpy()
        if active.any():
            leak += float(leaked[active].sum()); n += int(active.sum())
    return leak / n if n else 0.0


def _metrics(truth, pred, mp):
    cols = [c for c in truth.columns if not str(c).startswith("unresolved_")]
    t = truth[cols]; p = pred.reindex(index=t.index, columns=cols).fillna(0.0)
    acc = M.accuracy_metrics(t, p)
    tfam, pfam = M.aggregate_to_families(t, mp), M.aggregate_to_families(pred, mp)
    bacc = M.accuracy_metrics(tfam, pfam)
    tcond, pcond = _conditional(t, mp, cols), _conditional(p, mp, cols)
    cond = M.accuracy_metrics(tcond, pcond)
    tv, pv = t.to_numpy(float), p.to_numpy(float)
    absent = tv <= ABSENT; rare = (tv >= 0.005) & (tv <= 0.05)
    rare_detected = (pv[rare] > ABSENT)
    cx, cxt = M.composition_complexity(p), M.composition_complexity(t)
    ucols = [c for c in pred.columns if str(c).startswith("unresolved_")]
    umass = float(pred[ucols].to_numpy().sum() / max(pred.to_numpy().sum(), 1e-9)) if ucols else 0.0
    # rare precision: of predicted-present-but-truly-rare/absent, how many are truly rare
    pred_present_absent = pv[absent] > PRESENT
    return {
        "cond_rmse": cond["rmse"], "cond_pearson": cond["pearson"], "cond_spearman": cond["spearman"],
        "fine_rmse": acc["rmse"], "fine_pearson": acc["pearson"], "broad_rmse": bacc["rmse"],
        "pairwise_spillover": _pairwise_spillover(tcond, pcond, mp, cols),
        "false_resolution_rate": float(pv[absent].sum() / max(pv.sum(), 1e-9)),
        "false_positive_rate": float(pred_present_absent.mean()) if absent.any() else 0.0,
        "rare_sensitivity": float(rare_detected.mean()) if rare.any() else float("nan"),
        "rare_precision": float((rare.sum()) / max((pv > ABSENT).sum() - (tv > ABSENT).sum() + rare.sum(), 1))
        if rare.any() else float("nan"),
        "eff_n_pred": float(cx["effective_n_populations"].mean()),
        "eff_n_abs_error": float(abs(cx["effective_n_populations"].mean() - cxt["effective_n_populations"].mean())),
        "entropy_pred": float(cx["entropy"].mean()) if "entropy" in cx.columns else float("nan"),
        "unresolved_mass_frac": umass,
        "mass_error": float(abs(pred.sum(1) - truth.sum(1)).max()),
    }


def _perfamily_cond_rmse(truth, pred, mp):
    cols = [c for c in truth.columns if not str(c).startswith("unresolved_")]
    t = truth[cols]; p = pred.reindex(index=t.index, columns=cols).fillna(0.0)
    tcond, pcond = _conditional(t, mp, cols), _conditional(p, mp, cols)
    fam_of = {c: mp.get(c, c) for c in cols}
    out = {}
    for fam in set(fam_of.values()):
        mem = [c for c in cols if fam_of[c] == fam]
        if len(mem) < 2:
            continue
        present = (t[mem].sum(axis=1) > PRESENT).to_numpy()
        if present.any():
            err = tcond[mem].to_numpy()[present] - pcond[mem].to_numpy()[present]
            out[fam] = float(np.sqrt((err ** 2).mean()))
    return out


def _fit_confidence(adata, mapping, cfg):
    from tissueresolve.api import deconv_bulk
    from tissueresolve.reference.hierarchy import combine_family_and_conditional_estimates
    fine_col, donor_col = cfg["fine_col"], cfg["donor_col"]
    rows = []
    for seed in CAL:
        ref, ref_types, mp, qry, _ = _prep(adata, mapping, fine_col, donor_col, seed)
        rfeat = compute_reference_confidence_features(ref, mp)
        for scen in SCENARIOS:
            ds = _scenario(adata, ref_types, mp, qry, fine_col, donor_col, scen, seed)
            truth = ds.true_mrna_proportions
            with warnings.catch_warnings():
                warnings.simplefilter("ignore")
                h = deconv_bulk(ds.counts, ref, solver="nnls",
                                resolution_mode="hierarchical", hierarchy_mapping=mp)
            est = h.estimates
            pre = combine_family_and_conditional_estimates(
                est.family_proportions, est.conditional_proportions, mp)
            for samp in truth.index:
                for st in truth.columns:
                    if st not in rfeat.index:
                        continue
                    err = abs(float(pre.at[samp, st]) - float(truth.at[samp, st])) \
                        if st in pre.columns else float(truth.at[samp, st])
                    rows.append({"subtype": st, "sample": f"{seed}_{samp}",
                                 "ref_evidence": rfeat.at[st, "ref_evidence"],
                                 "correct": int(err < TOL)})
    fdf = pd.DataFrame(rows).dropna(subset=["ref_evidence"])
    return fit_confidence_model(
        fdf.set_index(["sample", "subtype"])[["ref_evidence"]],
        fdf.set_index(["sample", "subtype"])["correct"].astype(float),
        kind="logistic", score_col="ref_evidence")


def run_tissue(tissue, adata, mapping):
    from tissueresolve.api import deconv_bulk
    from tissueresolve.reference.hierarchy import combine_family_and_conditional_estimates
    cfg = TISSUES[tissue]
    fine_col, donor_col = cfg["fine_col"], cfg["donor_col"]
    print(f"\n[{tissue}] fitting confidence model on seeds {CAL} ...", flush=True)
    cmodel = _fit_confidence(adata, mapping, cfg)
    rows, perfam = [], []
    for seed in TEST:
        ref, ref_types, mp, qry, ref_d = _prep(adata, mapping, fine_col, donor_col, seed)
        conf_feat = compute_reference_confidence_features(ref, mp)[["ref_evidence"]]
        base_conf = cmodel.predict(conf_feat)
        rmask = adata.obs[donor_col].astype(str).isin(set(ref_d)).to_numpy()
        with warnings.catch_warnings():
            warnings.simplefilter("ignore")
            donor_res = build_donor_residuals(adata[rmask].copy(), ref, mp,
                                              celltype_col=fine_col, donor_col=donor_col)
        fdc = (adata[rmask].obs.groupby(
            adata[rmask].obs[fine_col].astype(str).map(mp))[donor_col].nunique().to_dict())
        for scen in SCENARIOS:
            ds = _scenario(adata, ref_types, mp, qry, fine_col, donor_col, scen, seed)
            truth = ds.true_mrna_proportions
            t0 = time.time()
            with warnings.catch_warnings():
                warnings.simplefilter("ignore")
                flat = deconv_bulk(ds.counts, ref, solver="nnls", resolution_mode="none").deconv.proportions
                h = deconv_bulk(ds.counts, ref, solver="nnls",
                                resolution_mode="hierarchical", hierarchy_mapping=mp)
            est = h.estimates; base_t = time.time() - t0
            raw_cond, fam_props = est.conditional_proportions, est.family_proportions
            pre = combine_family_and_conditional_estimates(fam_props, raw_cond, mp)
            soft_base = apply_partial_confidence_gating(fam_props, pre, mp, confidence=base_conf).combined

            for mode, pred in [("flat_nnls", flat), ("hard_gate", est.combined_fine),
                               ("ungated", pre), ("baseline_softgate", soft_base)]:
                rows.append({"tissue": tissue, "seed": seed, "scenario": scen, "mode": mode,
                             "runtime_s": base_t, **_metrics(truth, pred, mp)})

            for mode, kw in MODES.items():
                tr = time.time()
                strat = HighGranularityStrategy(HighGranularityConfig(**kw))
                rr = strat.refine(raw_cond, ref, mp, query=ds.counts,
                                  family_donor_counts=fdc, donor_residuals=donor_res)
                abs_ref = combine_family_and_conditional_estimates(fam_props, rr.refined_conditional, mp)
                # Component 5: confidence modulated by stability (subtype-level), gate ONCE
                conf = (base_conf * rr.stability_confidence.reindex(base_conf.index).fillna(1.0)).clip(0, 1)
                gated = apply_partial_confidence_gating(fam_props, abs_ref, mp, confidence=conf).combined
                rows.append({"tissue": tissue, "seed": seed, "scenario": scen, "mode": mode,
                             "runtime_s": base_t + (time.time() - tr), **_metrics(truth, gated, mp)})
                for fam, v in _perfamily_cond_rmse(truth, gated, mp).items():
                    perfam.append({"tissue": tissue, "seed": seed, "scenario": scen,
                                   "mode": mode, "family": fam, "cond_rmse": v})
            # baseline per-family for comparison
            for fam, v in _perfamily_cond_rmse(truth, soft_base, mp).items():
                perfam.append({"tissue": tissue, "seed": seed, "scenario": scen,
                               "mode": "baseline_softgate", "family": fam, "cond_rmse": v})
        print(f"[{tissue}] seed {seed} done", flush=True)
    return pd.DataFrame(rows), pd.DataFrame(perfam)


def evaluate_gates(modes_df, perfam_df):
    g = modes_df.groupby(["tissue", "mode"]).mean(numeric_only=True)
    rel = lambda a, b: (a - b) / abs(b) if b else float("inf")
    refiner_modes = list(MODES)
    rows, best_modes = [], {}
    for tissue in modes_df["tissue"].unique():
        base = g.loc[(tissue, "baseline_softgate")]
        best_mode = min(refiner_modes,
                        key=lambda m: g.loc[(tissue, m), "cond_rmse"] if (tissue, m) in g.index else 1e9)
        best = g.loc[(tissue, best_mode)]
        best_modes[tissue] = best_mode
        pf = perfam_df[perfam_df.tissue == tissue]
        base_fam = pf[pf["mode"] == "baseline_softgate"].groupby("family")["cond_rmse"].mean()
        best_fam = pf[pf["mode"] == best_mode].groupby("family")["cond_rmse"].mean()
        common = best_fam.index.intersection(base_fam.index)
        n_improved = int((best_fam[common] < base_fam[common] - 1e-6).sum()) if len(common) else 0
        checks = {
            "cond_rmse_rel_improve>=5%": rel(best["cond_rmse"], base["cond_rmse"]) <= -0.05,
            "pairwise_spillover_noninferior": best["pairwise_spillover"] <= base["pairwise_spillover"] + 0.01,
            "false_positive_noninferior": best["false_positive_rate"] <= base["false_positive_rate"] + 0.02,
            "rare_sensitivity_preserved": (best["rare_sensitivity"] >= base["rare_sensitivity"] - 0.02)
            if np.isfinite(best["rare_sensitivity"]) and np.isfinite(base["rare_sensitivity"]) else False,
            "broad_rmse_not_worse>2%": rel(best["broad_rmse"], base["broad_rmse"]) <= 0.02,
            "eff_n_not_inflated": best["eff_n_pred"] <= 1.5 * base["eff_n_pred"],
            "mass_conserved": best["mass_error"] < 1e-6,
            "runtime_acceptable": best["runtime_s"] < 5 * base["runtime_s"] + 10,
            "multi_family_improvement>1": n_improved > 1,
        }
        for name, ok in checks.items():
            rows.append({"tissue": tissue, "best_mode": best_mode, "gate": name, "pass": bool(ok)})
    gate_df = pd.DataFrame(rows)
    both = gate_df.pivot_table(index="gate", columns="tissue", values="pass", aggfunc="all")
    gate_df.attrs["both"] = both
    gate_df.attrs["best_modes"] = best_modes
    return gate_df, g


def main(argv=None) -> int:
    ap = argparse.ArgumentParser()
    ap.add_argument("--run-real-data", action="store_true")
    ap.add_argument("--tissue", choices=["breast", "lung", "both"], default="both")
    args = ap.parse_args(argv)
    if not args.run_real_data:
        print("Refusing without --run-real-data.", file=sys.stderr); return 2
    OUT.mkdir(parents=True, exist_ok=True)
    tissues = ["breast", "lung"] if args.tissue == "both" else [args.tissue]
    all_modes, all_perfam = [], []
    for tissue in tissues:
        cfg = TISSUES[tissue]
        if not Path(cfg["h5ad"]).exists():
            print(f"Missing {cfg['h5ad']}; skipping {tissue}.", file=sys.stderr); continue
        adata, mapping = _load(cfg)
        print(f"[{tissue}] {adata.shape}; donors={adata.obs[cfg['donor_col']].nunique()} "
              f"fine={adata.obs[cfg['fine_col']].nunique()}", flush=True)
        mdf, pf = run_tissue(tissue, adata, mapping)
        mdf.to_csv(OUT / f"high_granularity_{tissue}_modes.tsv", sep="\t", index=False)
        pf.to_csv(OUT / f"high_granularity_{tissue}_perfamily.tsv", sep="\t", index=False)
        all_modes.append(mdf); all_perfam.append(pf)
        print(f"\n=== [{tissue}] mode means ===")
        print(mdf.groupby("mode").mean(numeric_only=True)[
            ["cond_rmse", "cond_pearson", "fine_rmse", "broad_rmse", "pairwise_spillover",
             "rare_sensitivity", "false_positive_rate", "eff_n_pred", "runtime_s"]].round(4).to_string())
    if not all_modes:
        print("No tissue ran.", file=sys.stderr); return 2
    modes_df = pd.concat(all_modes, ignore_index=True)
    perfam_df = pd.concat(all_perfam, ignore_index=True)
    gate_df, g = evaluate_gates(modes_df, perfam_df)
    gate_df.to_csv(OUT / "high_granularity_gates.tsv", sep="\t", index=False)
    print("\n=== PROMOTION GATES (best variant vs baseline_softgate, per tissue) ===")
    print(gate_df.to_string(index=False))
    print(f"\nbest modes per tissue: {gate_df.attrs['best_modes']}")
    both = gate_df.attrs["both"]
    print("\n=== gate passes in BOTH tissues ===")
    print(both.to_string())
    all_pass = bool(both.all(axis=1).all()) and len(both.columns) >= 2
    print(f"\nALL GATES PASS IN BOTH TISSUES: {all_pass}")
    print("DECISION: " + ("INTEGRATE" if all_pass else
                          "KEEP EXPERIMENTAL — gates not met on both tissues; do not integrate"))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
