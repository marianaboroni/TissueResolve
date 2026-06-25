#!/usr/bin/env python
"""Evaluate the experimental FineGranularityRefiner on two tissues (breast + lung).

Tests whether a compact fine-level refinement layer improves *conditional
within-family* subtype prediction in collinear families, applied BEFORE the
validated partial confidence-weighted soft gate (refine → θ=π·q_refined →
confidence → soft gate, applied exactly once).

Design (donor- AND seed-disjoint; rules 12–13):
  * reference/query donor split per seed (synthetic_holdout.split_donors);
  * contrast weights + donor-stability residuals use REFERENCE donors only;
  * confidence + spillover calibration fit on CALIBRATION seeds only;
  * all reported metrics from held-out TEST seeds only.

Methods compared (all but flat/hard/ungated share the SAME soft gate):
  baseline_softgate  — current hierarchical fine + soft gating
  contrast_softgate  — contrast-weighted WNNLS + soft gating
  residual_softgate  — residual-contrast refinement + soft gating
  spillover_softgate — contrast-weighted + family spillover calibration + soft gating
  flat_nnls          — flat fine NNLS (reference)
  hard_gate          — legacy binary gate (reference)
  ungated            — diagnostic, no abstention (reference)

Outputs (benchmarks/outputs/): fine_refiner_<tissue>_modes.tsv,
fine_refiner_<tissue>_perfamily.tsv, fine_refiner_gates.tsv.
NOTHING is committed; inputs are git-ignored real data.

Usage:
  python fine_refiner_benchmark.py --run-real-data --tissue both
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
    apply_partial_confidence_gating, compute_reference_confidence_features, fit_confidence_model)
from tissueresolve.experimental.soft_hierarchy.fine_refiner import (  # noqa: E402
    FineGranularityRefiner, FineRefinerConfig, SpilloverCalibrator, build_donor_residuals)

OUT = REPO / "benchmarks" / "outputs"
CAL, TEST = [0, 1], [2, 3]
SCENARIOS = ["balanced", "imbalanced", "rare", "similar_subtypes"]
N_SAMPLES, CELLS, MIN_CELLS = 10, 400, 20
PRESENT, ABSENT, TOL = 0.01, 0.005, 0.05

# Reference cell types absent from the shipped breast hierarchy tsv, mapped into
# its existing broad families here (the config file is left untouched until/unless
# the method passes the promotion gates).  All assignments use the same families
# already present in the tsv.
_BREAST_SUPPLEMENT = {
    "IgG plasma cell": "B/Plasma",
    "naive B cell": "B/Plasma",
    "unswitched memory B cell": "B/Plasma",
    "activated CD8-positive, alpha-beta T cell": "T/NK",
    "gamma-delta T cell": "T/NK",
    "mast cell": "Myeloid",
    "myeloid dendritic cell": "Myeloid",
    "neutrophil": "Myeloid",
    "plasmacytoid dendritic cell": "Myeloid",
}

TISSUES = {
    "breast": {
        "h5ad": REPO / "examples/real_breast_cancer/data/reference/breast_cancer_sc_reference.h5ad",
        "hmap": REPO / "examples/real_breast_cancer/config/breast_cancer_cell_type_hierarchy.tsv",
        "fine_col": "cell_type", "donor_col": "donor_id", "feature_name": True,
        "supplement": _BREAST_SUPPLEMENT,
    },
    "lung": {
        "h5ad": REPO / "examples/second_tissue_lung/data/derived/hlca_subset.h5ad",
        "hmap": REPO / "examples/second_tissue_lung/data/derived/hlca_hierarchy.tsv",
        "fine_col": "cell_type_fine", "donor_col": "donor_id", "feature_name": True,
    },
}


# --------------------------------------------------------------------- helpers
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
        prep = H.prepare_reference(adata[rmask].copy(), cell_type_col=fine_col,
                                   min_cells=MIN_CELLS, estimate_overdispersion=False)
    ref = prep.reference
    ref_types = [str(c) for c in ref.cell_types]
    mp = {k: v for k, v in build_cell_type_hierarchy(ref_types, mapping).items()}
    return ref, ref_types, mp, qry, ref_d


def _scenario(adata, ref_types, mp, qry, fine_col, donor_col, scen, seed):
    kw = {}
    if scen == "similar_subtypes":
        fam = max({mp.get(x, x) for x in ref_types},
                  key=lambda f: sum(mp.get(x, x) == f for x in ref_types))
        pair = [t for t in ref_types if mp.get(t) == fam][:2]
        kw = {"similar_pair": tuple(pair)}
    elif scen == "rare":
        kw = {"rare_type": ref_types[-1], "rare_level": 0.01}
    tgt = SH.build_target_proportions(ref_types, N_SAMPLES, scen, seed=1000 + seed, **kw)
    return SH.realize_pseudobulk(adata, tgt, celltype_col=fine_col, donor_col=donor_col,
                                 query_donors=qry, seed=seed, cells_per_sample=CELLS)


def _conditional(df, mp, cols):
    """Within-family conditional proportions (members renormalised to sum 1)."""
    fam_of = {c: mp.get(c, c) for c in cols}
    out = df[cols] * 0.0
    for fam in set(fam_of.values()):
        mem = [c for c in cols if fam_of[c] == fam]
        s = df[mem].sum(axis=1)
        for c in mem:
            out[c] = np.where(s > 0, df[c] / s.replace(0, np.nan), 0.0)
    return out.fillna(0.0)


def _pairwise_spillover(tcond, pcond, mp, cols):
    """Mean conditional mass leaked to truly-absent subtypes within families."""
    fam_of = {c: mp.get(c, c) for c in cols}
    leak, n = 0.0, 0
    for fam in set(fam_of.values()):
        mem = [c for c in cols if fam_of[c] == fam]
        if len(mem) < 2:
            continue
        present = (tcond[mem] > PRESENT)
        absent = ~present
        leaked = (pcond[mem].to_numpy() * absent.to_numpy()).sum(axis=1)
        active = present.any(axis=1).to_numpy()
        if active.any():
            leak += float(leaked[active].sum()); n += int(active.sum())
    return leak / n if n else 0.0


def _metrics(truth, pred, mp, mode):
    cols = [c for c in truth.columns if not str(c).startswith("unresolved_")]
    t = truth[cols]; p = pred.reindex(index=t.index, columns=cols).fillna(0.0)
    acc = M.accuracy_metrics(t, p)
    tfam, pfam = M.aggregate_to_families(t, mp), M.aggregate_to_families(pred, mp)
    bacc = M.accuracy_metrics(tfam, pfam)
    tcond, pcond = _conditional(t, mp, cols), _conditional(p, mp, cols)
    cond = M.accuracy_metrics(tcond, pcond)
    tv, pv = t.to_numpy(float), p.to_numpy(float)
    absent = tv <= ABSENT; rare = (tv >= 0.005) & (tv <= 0.05)
    cx, cxt = M.composition_complexity(p), M.composition_complexity(t)
    ucols = [c for c in pred.columns if str(c).startswith("unresolved_")]
    umass = float(pred[ucols].to_numpy().sum() / max(pred.to_numpy().sum(), 1e-9)) if ucols else 0.0
    row = {
        "cond_rmse": cond["rmse"], "cond_pearson": cond["pearson"], "cond_spearman": cond["spearman"],
        "fine_rmse": acc["rmse"], "fine_pearson": acc["pearson"], "broad_rmse": bacc["rmse"],
        "pairwise_spillover": _pairwise_spillover(tcond, pcond, mp, cols),
        "false_resolution_rate": float(pv[absent].sum() / max(pv.sum(), 1e-9)),
        "false_positive_rate": float((pv[absent] > PRESENT).mean()) if absent.any() else 0.0,
        "rare_sensitivity": float((pv[rare] > ABSENT).mean()) if rare.any() else float("nan"),
        "eff_n_pred": float(cx["effective_n_populations"].mean()),
        "eff_n_abs_error": float(abs(cx["effective_n_populations"].mean() - cxt["effective_n_populations"].mean())),
        "entropy_pred": float(cx["entropy"].mean()) if "entropy" in cx.columns else float("nan"),
        "unresolved_mass_frac": umass,
        "mass_error": float(abs(pred.sum(1) - truth.sum(1)).max()),
    }
    return row


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
        if not present.any():
            continue
        err = (tcond[mem].to_numpy()[present] - pcond[mem].to_numpy()[present])
        out[fam] = float(np.sqrt((err ** 2).mean()))
    return out


# --------------------------------------------------------------------- calibration
def _calibrate(adata, mapping, cfg, tissue):
    """Fit confidence model + spillover calibrator on CALIBRATION seeds only."""
    from tissueresolve.api import deconv_bulk
    from tissueresolve.reference.hierarchy import combine_family_and_conditional_estimates
    fine_col, donor_col = cfg["fine_col"], cfg["donor_col"]
    feat_rows = []
    pred_by_fam: dict = {}
    true_by_fam: dict = {}
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
            # confidence-calibration rows
            for samp in truth.index:
                for st in truth.columns:
                    if st not in rfeat.index:
                        continue
                    err = abs(float(pre.at[samp, st]) - float(truth.at[samp, st])) \
                        if st in pre.columns else float(truth.at[samp, st])
                    feat_rows.append({"subtype": st, "sample": f"{seed}_{samp}",
                                      "ref_evidence": rfeat.at[st, "ref_evidence"],
                                      "correct": int(err < TOL)})
            # spillover-calibration rows (conditional pred vs true, per family)
            cols = [c for c in truth.columns if c in ref_types]
            tcond = _conditional(truth, mp, cols)
            pcond = _conditional(pre.reindex(columns=cols).fillna(0.0), mp, cols)
            fam_of = {c: mp.get(c, c) for c in cols}
            for fam in set(fam_of.values()):
                mem = [c for c in cols if fam_of[c] == fam]
                if len(mem) < 2:
                    continue
                present = (truth[mem].sum(axis=1) > PRESENT).to_numpy()
                if not present.any():
                    continue
                pred_by_fam.setdefault(fam, []).append(pcond[mem].iloc[present].copy())
                true_by_fam.setdefault(fam, []).append(tcond[mem].iloc[present].copy())
    fdf = pd.DataFrame(feat_rows).dropna(subset=["ref_evidence"])
    cmodel = fit_confidence_model(
        fdf.set_index(["sample", "subtype"])[["ref_evidence"]],
        fdf.set_index(["sample", "subtype"])["correct"].astype(float),
        kind="logistic", score_col="ref_evidence")
    pred_f = {f: pd.concat(v, axis=0) for f, v in pred_by_fam.items()}
    true_f = {f: pd.concat(v, axis=0) for f, v in true_by_fam.items()}
    calibrator = SpilloverCalibrator(kind="ridge", ridge_alpha=1.0).fit(pred_f, true_f)
    return cmodel, calibrator


# --------------------------------------------------------------------- per-tissue run
def run_tissue(tissue, adata, mapping):
    from tissueresolve.api import deconv_bulk
    from tissueresolve.reference.hierarchy import combine_family_and_conditional_estimates
    cfg = TISSUES[tissue]
    fine_col, donor_col = cfg["fine_col"], cfg["donor_col"]
    print(f"\n[{tissue}] calibrating confidence + spillover on seeds {CAL} ...", flush=True)
    cmodel, calibrator = _calibrate(adata, mapping, cfg, tissue)
    print(f"[{tissue}] spillover-calibrated families: {sorted(calibrator.maps_.keys())}", flush=True)

    modes = {
        "baseline_softgate": FineRefinerConfig(mode="none"),
        "contrast_softgate": FineRefinerConfig(mode="contrast_weighted"),
        "residual_softgate": FineRefinerConfig(mode="residual_contrast"),
        "spillover_softgate": FineRefinerConfig(mode="spillover_calibrated", calibration="ridge"),
    }
    rows, perfam_rows = [], []
    for seed in TEST:
        ref, ref_types, mp, qry, ref_d = _prep(adata, mapping, fine_col, donor_col, seed)
        conf_feat = compute_reference_confidence_features(ref, mp)[["ref_evidence"]]
        conf = cmodel.predict(conf_feat)
        # donor-stability residuals from REFERENCE donors only (no test-donor leak)
        rmask = adata.obs[donor_col].astype(str).isin(set(ref_d)).to_numpy()
        with warnings.catch_warnings():
            warnings.simplefilter("ignore")
            donor_res = build_donor_residuals(adata[rmask].copy(), ref, mp,
                                              celltype_col=fine_col, donor_col=donor_col)
        fam_donor_counts = (adata[rmask].obs.groupby(
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
            est = h.estimates
            base_t = time.time() - t0
            raw_cond = est.conditional_proportions
            fam_props = est.family_proportions
            pre = combine_family_and_conditional_estimates(fam_props, raw_cond, mp)

            # reference modes (no refiner)
            for mode, pred in [("flat_nnls", flat), ("hard_gate", est.combined_fine),
                               ("ungated", pre)]:
                rows.append({"tissue": tissue, "seed": seed, "scenario": scen,
                             "mode": mode, "runtime_s": base_t, **_metrics(truth, pred, mp, mode)})

            # refiner modes (refine → combine → soft gate)
            for mode, rcfg in modes.items():
                tr = time.time()
                refiner = FineGranularityRefiner(
                    rcfg, calibrator=calibrator if rcfg.mode == "spillover_calibrated" else None)
                rr = refiner.refine(raw_cond, ref, mp, query=ds.counts,
                                    family_donor_counts=fam_donor_counts,
                                    donor_residuals=donor_res)
                abs_refined = combine_family_and_conditional_estimates(
                    fam_props, rr.refined_conditional, mp)
                gated = apply_partial_confidence_gating(
                    fam_props, abs_refined, mp, confidence=conf).combined
                rt = base_t + (time.time() - tr)
                rows.append({"tissue": tissue, "seed": seed, "scenario": scen,
                             "mode": mode, "runtime_s": rt, **_metrics(truth, gated, mp, mode)})
                # per-family conditional RMSE
                for fam, v in _perfamily_cond_rmse(truth, gated, mp).items():
                    perfam_rows.append({"tissue": tissue, "seed": seed, "scenario": scen,
                                        "mode": mode, "family": fam, "cond_rmse": v})
        print(f"[{tissue}] seed {seed} done", flush=True)
    return pd.DataFrame(rows), pd.DataFrame(perfam_rows)


# --------------------------------------------------------------------- gates
def evaluate_gates(modes_df, perfam_df):
    """Best-refiner vs baseline_softgate, per tissue + combined gate table."""
    g = modes_df.groupby(["tissue", "mode"]).mean(numeric_only=True)
    rel = lambda a, b: (a - b) / abs(b) if b else float("inf")
    refiner_modes = ["contrast_softgate", "residual_softgate", "spillover_softgate"]
    gate_rows, per_tissue_best = [], {}
    for tissue in modes_df["tissue"].unique():
        base = g.loc[(tissue, "baseline_softgate")]
        # best refiner = lowest conditional RMSE
        best_mode = min(refiner_modes,
                        key=lambda m: g.loc[(tissue, m), "cond_rmse"]
                        if (tissue, m) in g.index else float("inf"))
        best = g.loc[(tissue, best_mode)]
        per_tissue_best[tissue] = best_mode
        # per-family improvement count
        pf = perfam_df[(perfam_df.tissue == tissue)]
        base_fam = pf[pf["mode"] == "baseline_softgate"].groupby("family")["cond_rmse"].mean() \
            if "baseline_softgate" in pf["mode"].unique() else pd.Series(dtype=float)
        best_fam = pf[pf["mode"] == best_mode].groupby("family")["cond_rmse"].mean()
        # baseline per-family from modes (baseline has no perfam rows → recompute via best vs itself)
        n_fam_improved = 0
        if not base_fam.empty:
            common = best_fam.index.intersection(base_fam.index)
            n_fam_improved = int((best_fam[common] < base_fam[common] - 1e-6).sum())
        checks = {
            "cond_rmse_rel_improve>=5%": rel(best["cond_rmse"], base["cond_rmse"]) <= -0.05,
            "pairwise_spillover_decreases": best["pairwise_spillover"] < base["pairwise_spillover"],
            "rare_sensitivity_noninferior": (best["rare_sensitivity"] >= base["rare_sensitivity"] - 0.02)
                if np.isfinite(best["rare_sensitivity"]) and np.isfinite(base["rare_sensitivity"]) else False,
            "false_positive_not_worse": best["false_positive_rate"] <= base["false_positive_rate"] + 0.02,
            "broad_rmse_not_worse>2%": rel(best["broad_rmse"], base["broad_rmse"]) <= 0.02,
            "eff_n_not_collapsed": best["eff_n_pred"] >= 0.6 * base["eff_n_pred"],
            "mass_conserved": best["mass_error"] < 1e-6,
            "runtime_acceptable": best["runtime_s"] < 5 * base["runtime_s"] + 5,
            "multi_family_improvement>1": n_fam_improved > 1,
        }
        for name, ok in checks.items():
            gate_rows.append({"tissue": tissue, "best_mode": best_mode, "gate": name,
                              "pass": bool(ok)})
    gate_df = pd.DataFrame(gate_rows)
    # tissue-agnostic gate: improvement must occur in BOTH tissues
    both = (gate_df.pivot_table(index="gate", columns="tissue", values="pass", aggfunc="all"))
    gate_df.attrs["both_tissues"] = both
    gate_df.attrs["best_modes"] = per_tissue_best
    return gate_df, g


def main(argv=None) -> int:
    ap = argparse.ArgumentParser()
    ap.add_argument("--run-real-data", action="store_true")
    ap.add_argument("--tissue", choices=["breast", "lung", "both"], default="both")
    args = ap.parse_args(argv)
    if not args.run_real_data:
        print("Refusing without --run-real-data.", file=sys.stderr)
        return 2
    OUT.mkdir(parents=True, exist_ok=True)
    tissues = ["breast", "lung"] if args.tissue == "both" else [args.tissue]
    all_modes, all_perfam = [], []
    for tissue in tissues:
        cfg = TISSUES[tissue]
        if not Path(cfg["h5ad"]).exists():
            print(f"Missing {cfg['h5ad']}; skipping {tissue}.", file=sys.stderr)
            continue
        adata, mapping = _load(cfg)
        print(f"[{tissue}] {adata.shape}; donors={adata.obs[cfg['donor_col']].nunique()} "
              f"fine={adata.obs[cfg['fine_col']].nunique()}", flush=True)
        mdf, pf = run_tissue(tissue, adata, mapping)
        mdf.to_csv(OUT / f"fine_refiner_{tissue}_modes.tsv", sep="\t", index=False)
        pf.to_csv(OUT / f"fine_refiner_{tissue}_perfamily.tsv", sep="\t", index=False)
        all_modes.append(mdf); all_perfam.append(pf)
        print(f"\n=== [{tissue}] mode means ===")
        print(mdf.groupby("mode").mean(numeric_only=True)[
            ["cond_rmse", "cond_pearson", "fine_rmse", "broad_rmse", "pairwise_spillover",
             "rare_sensitivity", "false_positive_rate", "eff_n_pred", "runtime_s"]].round(4).to_string())
    if not all_modes:
        print("No tissue ran.", file=sys.stderr)
        return 2
    modes_df = pd.concat(all_modes, ignore_index=True)
    perfam_df = pd.concat(all_perfam, ignore_index=True)
    gate_df, g = evaluate_gates(modes_df, perfam_df)
    gate_df.to_csv(OUT / "fine_refiner_gates.tsv", sep="\t", index=False)
    print("\n=== PROMOTION GATES (best refiner vs baseline_softgate, per tissue) ===")
    print(gate_df.to_string(index=False))
    print(f"\nbest modes per tissue: {gate_df.attrs['best_modes']}")
    both = gate_df.attrs["both_tissues"]
    print("\n=== gate passes in BOTH tissues ===")
    print(both.to_string())
    all_pass_both = bool(both.all(axis=1).all()) and len(both.columns) >= 2
    print(f"\nALL GATES PASS IN BOTH TISSUES: {all_pass_both}")
    print("DECISION: " + ("INTEGRATE (passes all gates on both tissues)"
                          if all_pass_both else
                          "KEEP EXPERIMENTAL — gates not met on both tissues; do not integrate"))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
