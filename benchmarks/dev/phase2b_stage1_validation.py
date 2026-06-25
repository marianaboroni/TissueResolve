#!/usr/bin/env python
"""PHASE 2B Stage 1 — independent (expanded) validation of soft gating.

Breast atlas only (no adequate second tissue available — recorded as a limitation;
2-tissue replication therefore UNMET → promotion blocked). Expands Phase 2A:
calibration seeds {0,1}, validation {2}, TEST seeds {3..12} (10), scenarios
{balanced, imbalanced, similar_subtypes, rare, missing_population,
reduced_gene_overlap}, modes {flat, hard_gate, ungated, soft_gate}. Donor- AND
seed-disjoint (rules 3,4). Evaluates the PROSPECTIVE gates declared in
PHASE2B_IMPLEMENTATION_AND_VALIDATION_PLAN.md (revised false-resolution criterion =
non-inferior to ungated).

Outputs: benchmarks/outputs/phase2a_independent_validation.tsv (raw per-replicate),
phase2a_tissue_summary.tsv (per-mode CI), phase2a_resolution_metrics.tsv.

Usage:  python benchmarks/diagnostics/phase2b_stage1_validation.py --run-real-data
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
    fit_confidence_model)

OUT = REPO / "benchmarks" / "outputs"
HIER = REPO / "examples" / "real_breast_cancer" / "config" / "breast_cancer_cell_type_hierarchy.tsv"
CAL, VAL, TEST = [0, 1], [2], list(range(3, 13))
SCENARIOS = ["balanced", "imbalanced", "similar_subtypes", "rare", "missing_population",
             "reduced_gene_overlap"]
SCEN_SEED = {"balanced": 11, "imbalanced": 22, "similar_subtypes": 44, "rare": 33,
             "missing_population": 55, "reduced_gene_overlap": 66}
N_SAMPLES, CELLS, MIN_CELLS = 12, 600, 30
PRESENT, ABSENT, TOL = 0.01, 0.005, 0.05
SCORE = "ref_evidence"


def _prep(adata, raw_map, seed):
    from tissueresolve.reference.hierarchy import build_cell_type_hierarchy
    ref_d, qry = SH.split_donors(adata, "donor_id", ref_frac=0.5, seed=seed)
    rmask = adata.obs["donor_id"].astype(str).isin(set(ref_d)).to_numpy()
    # restrict reference cells to mapped cell types (some donor subsets surface
    # atlas types like 'mast cell' not in the hierarchy mapping)
    mapped = adata.obs["cell_type"].astype(str).isin(set(raw_map)).to_numpy()
    rmask = rmask & mapped
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
    elif scen == "missing_population":
        kw = {"missing_type": next((t for t in ref_types if mapping.get(t) == "Epithelial"),
                                   ref_types[0])}
    base = "missing_population" if scen == "reduced_gene_overlap" else scen
    sseed = SCEN_SEED[scen]
    tgt = SH.build_target_proportions(ref_types, N_SAMPLES,
                                      "imbalanced" if scen == "reduced_gene_overlap" else base,
                                      seed=sseed, **({} if scen == "reduced_gene_overlap" else kw))
    ds = SH.realize_pseudobulk(adata, tgt, celltype_col="cell_type", donor_col="donor_id",
                               query_donors=qry, seed=seed, cells_per_sample=CELLS)
    counts = ds.counts
    if scen == "reduced_gene_overlap":
        rng = np.random.default_rng(seed * 100 + 7)
        keep = rng.choice(counts.index, size=int(0.4 * len(counts.index)), replace=False)
        counts = counts.loc[sorted(keep)]
    return ds.true_mrna_proportions, counts


def _resolution_metrics(truth, combined, mapping):
    """Truth-based unresolved precision/recall + false-resolution/abstention.
    should-abstain family = <2 present subtypes (truth>PRESENT)."""
    cols = [c for c in truth.columns if not str(c).startswith("unresolved_")]
    fam_of = {c: mapping.get(c, c) for c in cols}
    ucols = {c[len("unresolved_"):]: c for c in combined.columns if str(c).startswith("unresolved_")}
    u_should = u_total = should_abstain_mass = res_correct = n_fam = 0.0
    fa_num = fa_den = 0.0
    for fam in set(fam_of.values()):
        mem = [c for c in cols if fam_of[c] == fam]
        present_cnt = (truth[mem] > PRESENT).sum(1)
        fam_mass = truth[mem].sum(1)
        u = combined[ucols[fam]] if fam in ucols else pd.Series(0.0, index=truth.index)
        res = combined[mem].sum(1) if all(m in combined.columns for m in mem) else pd.Series(0.0, index=truth.index)
        should_abstain = present_cnt < 2
        should_resolve = present_cnt >= 2
        u_total += float(u.sum())
        u_should += float(u[should_abstain].sum())
        should_abstain_mass += float(fam_mass[should_abstain].sum())
        # false abstention: unresolved mass in should-resolve families
        fa_num += float(u[should_resolve].sum()); fa_den += float(fam_mass[should_resolve].sum())
        # correct-resolution-level: per (family,sample) decision matches
        decided_resolve = res > u
        ok = ((should_resolve & decided_resolve) | (should_abstain & ~decided_resolve))
        res_correct += float(ok.sum()); n_fam += float(len(truth.index))
    return {
        "unresolved_precision": float(u_should / u_total) if u_total > 1e-9 else float("nan"),
        "unresolved_recall": float(u_should / should_abstain_mass) if should_abstain_mass > 1e-9 else float("nan"),
        "false_abstention_rate": float(fa_num / fa_den) if fa_den > 1e-9 else float("nan"),
        "correct_resolution_level": float(res_correct / n_fam) if n_fam else float("nan"),
    }


def _metrics(truth, pred, mapping, mode):
    cols = [c for c in truth.columns if not str(c).startswith("unresolved_")]
    t = truth[cols]; p = pred.reindex(index=t.index, columns=cols).fillna(0.0)
    acc = M.accuracy_metrics(t, p); comp = M.compositional_metrics(t, p)
    tfam = M.aggregate_to_families(t, mapping); pfam = M.aggregate_to_families(pred, mapping)
    bacc = M.accuracy_metrics(tfam, pfam); bcomp = M.compositional_metrics(tfam, pfam)
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
    absent = tv <= ABSENT
    rare = (tv >= 0.005) & (tv <= 0.05)
    cx = M.composition_complexity(p)
    cxt = M.composition_complexity(t)
    ucols = [c for c in pred.columns if str(c).startswith("unresolved_")]
    umass = float(pred[ucols].to_numpy().sum() / max(pred.to_numpy().sum(), 1e-9)) if ucols else 0.0
    row = {
        "broad_pearson": bacc["pearson"], "broad_spearman": bacc["spearman"],
        "broad_rmse": bacc["rmse"], "broad_mae": bacc["mae"],
        "broad_jsd": bcomp["jsd_mean"], "broad_aitchison": bcomp["aitchison_mean"],
        "fine_pearson": acc["pearson"], "fine_rmse": acc["rmse"],
        "cond_pearson": cond["pearson"], "cond_rmse": cond["rmse"],
        "rare_sensitivity": float((pv[rare] > ABSENT).mean()) if rare.any() else float("nan"),
        "false_positive_detection": float((pv[absent] > PRESENT).mean()) if absent.any() else float("nan"),
        "false_resolution_rate": float(pv[absent].sum() / max(pv.sum(), 1e-9)),
        "true_eff_n": float(cxt["effective_n_populations"].mean()),
        "pred_eff_n": float(cx["effective_n_populations"].mean()),
        "eff_n_abs_error": float(abs(cx["effective_n_populations"].mean() - cxt["effective_n_populations"].mean())),
        "entropy_pred": float(cx["shannon_entropy_nats"].mean()),
        "entropy_true": float(cxt["shannon_entropy_nats"].mean()),
        "richness_pred": float((p > PRESENT).sum(1).mean()),
        "gini_pred": float(cx["gini"].mean()),
        "dominant_pred": float(cx["dominant_fraction"].mean()),
        "unresolved_mass_frac": umass,
        "mass_error": float(abs(pred.sum(1) - truth.sum(1)).max()),
        "n_negative": int((pred.to_numpy() < -1e-9).sum()),
    }
    if mode in ("hard_gate", "soft_gate"):
        row.update(_resolution_metrics(truth, pred, mapping))
    else:
        row.update({"unresolved_precision": float("nan"), "unresolved_recall": float("nan"),
                    "false_abstention_rate": 0.0, "correct_resolution_level": float("nan")})
    return row


def main(argv=None) -> int:
    ap = argparse.ArgumentParser(); ap.add_argument("--run-real-data", action="store_true")
    args = ap.parse_args(argv)
    if not H.real_data_enabled(args.run_real_data):
        print("Refusing without --run-real-data.", file=sys.stderr); return 2
    import anndata as ad
    from tissueresolve.api import deconv_bulk
    from tissueresolve.reference.hierarchy import (load_hierarchy_mapping,
                                                   combine_family_and_conditional_estimates)
    adata = ad.read_h5ad(H.REFERENCE_H5AD)
    if "feature_name" in adata.var.columns:
        fn = adata.var["feature_name"].astype(str)
        if fn.nunique() == adata.n_vars:
            adata.var_names = fn.to_numpy()
    raw_map = load_hierarchy_mapping(HIER)

    # ---- calibrate confidence on CAL seeds (ref_evidence → correct-resolution) ----
    feat_rows = []
    for seed in CAL:
        ref, ref_types, mapping, qry = _prep(adata, raw_map, seed)
        rfeat = compute_reference_confidence_features(ref, mapping)
        for scen in SCENARIOS:
            truth, counts = _scenario(adata, ref_types, mapping, qry, scen, seed)
            with warnings.catch_warnings():
                warnings.simplefilter("ignore")
                h = deconv_bulk(counts, ref, solver="auto", resolution_mode="hierarchical",
                                hierarchy_mapping=mapping)
            pre = combine_family_and_conditional_estimates(
                h.estimates.family_proportions, h.estimates.conditional_proportions, mapping)
            for samp in truth.index:
                for st in [c for c in truth.columns]:
                    if st not in rfeat.index:
                        continue
                    err = abs(float(pre.at[samp, st]) - float(truth.at[samp, st])) if st in pre.columns else float(truth.at[samp, st])
                    feat_rows.append({"subtype": st, "sample": samp, SCORE: rfeat.at[st, SCORE],
                                      "correct": int(err < TOL)})
    fdf = pd.DataFrame(feat_rows).dropna(subset=[SCORE])
    conf_model = fit_confidence_model(
        fdf.set_index(["sample", "subtype"])[[SCORE]],
        fdf.set_index(["sample", "subtype"])["correct"].astype(float),
        kind="logistic", score_col=SCORE)
    print(f"Calibrated confidence (logistic) on {len(fdf)} cal pairs.", flush=True)

    # ---- TEST ----
    rows = []
    for seed in TEST:
        ref, ref_types, mapping, qry = _prep(adata, raw_map, seed)
        conf_feat = compute_reference_confidence_features(ref, mapping)[[SCORE]]
        for scen in SCENARIOS:
            truth, counts = _scenario(adata, ref_types, mapping, qry, scen, seed)
            with warnings.catch_warnings():
                warnings.simplefilter("ignore")
                t0 = time.perf_counter()
                flat = deconv_bulk(counts, ref, solver="nnls", resolution_mode="none").deconv.proportions
                rt_flat = time.perf_counter() - t0
                t0 = time.perf_counter()
                h = deconv_bulk(counts, ref, solver="auto", resolution_mode="hierarchical",
                                hierarchy_mapping=mapping)
                rt_h = time.perf_counter() - t0
            est = h.estimates
            pre = combine_family_and_conditional_estimates(est.family_proportions,
                                                           est.conditional_proportions, mapping)
            conf = conf_model.predict(conf_feat)
            soft = apply_partial_confidence_gating(est.family_proportions, pre, mapping, confidence=conf).combined
            preds = {"flat": (flat, rt_flat), "hard_gate": (est.combined_fine, rt_h),
                     "ungated": (pre, 0.0), "soft_gate": (soft, 0.0)}
            for mode, (pred, rt) in preds.items():
                rows.append({"tissue": "breast", "seed": seed, "scenario": scen, "mode": mode,
                             "runtime_s": round(rt, 3), **_metrics(truth, pred, mapping, mode)})
        print(f"  seed{seed} done", flush=True)

    df = pd.DataFrame(rows)
    OUT.mkdir(parents=True, exist_ok=True)
    with open(OUT / "phase2a_independent_validation.tsv", "w") as fh:
        fh.write("# Phase 2A independent (expanded) validation. Breast only — 2-tissue replication UNMET (no 2nd atlas).\n")
        df.to_csv(fh, sep="\t", index=False)

    # per-mode CI summary
    gate_metrics = ["broad_rmse", "fine_rmse", "fine_pearson", "cond_pearson", "cond_rmse",
                    "false_resolution_rate", "false_abstention_rate", "eff_n_abs_error",
                    "rare_sensitivity", "unresolved_mass_frac", "unresolved_precision",
                    "mass_error", "broad_jsd", "broad_aitchison"]
    ci = []
    for mode in ["flat", "hard_gate", "ungated", "soft_gate"]:
        sub = df[df["mode"] == mode]
        for met in gate_metrics:
            c = ST.bootstrap_ci(sub[met].dropna().to_numpy(), seed=0)
            ci.append({"mode": mode, "metric": met, "point": round(c["point"], 4),
                       "lo": round(c["lo"], 4), "hi": round(c["hi"], 4), "n": c["n"]})
    cidf = pd.DataFrame(ci)
    with open(OUT / "phase2a_tissue_summary.tsv", "w") as fh:
        fh.write("# Per-mode bootstrap 95% CI over 10 seeds × 6 scenarios (breast). Single tissue.\n")
        cidf.to_csv(fh, sep="\t", index=False)
    res = df[df["mode"].isin(["hard_gate", "soft_gate"])].groupby("mode")[
        ["unresolved_precision", "unresolved_recall", "false_resolution_rate",
         "false_abstention_rate", "correct_resolution_level", "unresolved_mass_frac"]].mean().round(4)
    res.to_csv(OUT / "phase2a_resolution_metrics.tsv", sep="\t")

    # ---- prospective gates (soft vs ungated + vs hard) ----
    g = df.groupby("mode").mean(numeric_only=True)
    soft, ung, hard = g.loc["soft_gate"], g.loc["ungated"], g.loc["hard_gate"]
    rel = lambda a, b: (a - b) / abs(b) if b else float("inf")
    gates = [
        ("vs ungated: broad RMSE ≤ +2%", soft["broad_rmse"], ung["broad_rmse"], rel(soft["broad_rmse"], ung["broad_rmse"]) <= 0.02),
        ("vs ungated: fine RMSE ≤ +2%", soft["fine_rmse"], ung["fine_rmse"], rel(soft["fine_rmse"], ung["fine_rmse"]) <= 0.02),
        ("vs ungated: false-resolution ≤ +10% (revised)", soft["false_resolution_rate"], ung["false_resolution_rate"], rel(soft["false_resolution_rate"], ung["false_resolution_rate"]) <= 0.10),
        ("vs ungated: eff-N error not worse", soft["eff_n_abs_error"], ung["eff_n_abs_error"], soft["eff_n_abs_error"] <= ung["eff_n_abs_error"] * 1.05),
        ("vs ungated: mass conserved <1e-6", soft["mass_error"], 1e-6, soft["mass_error"] < 1e-6),
        ("vs hard: false-abstention materially lower", soft["false_abstention_rate"], hard["false_abstention_rate"], soft["false_abstention_rate"] < 0.5 * hard["false_abstention_rate"]),
        ("vs hard: rare sensitivity materially higher", soft["rare_sensitivity"], hard["rare_sensitivity"], soft["rare_sensitivity"] > 1.5 * hard["rare_sensitivity"]),
        ("vs hard: eff-N error substantially lower", soft["eff_n_abs_error"], hard["eff_n_abs_error"], soft["eff_n_abs_error"] < 0.6 * hard["eff_n_abs_error"]),
        ("vs hard: fine RMSE lower", soft["fine_rmse"], hard["fine_rmse"], soft["fine_rmse"] < hard["fine_rmse"]),
    ]
    gate_df = pd.DataFrame([{"gate": n, "soft": round(float(a), 4), "ref": round(float(b), 6),
                             "pass": bool(p)} for n, a, b, p in gates])
    gate_df.to_csv(OUT / "phase2a_stage1_gates.tsv", sep="\t", index=False)
    n_pass = int(gate_df["pass"].sum())
    print("\n=== Stage-1 mode means (breast) ===")
    print(g.loc[["flat", "ungated", "soft_gate", "hard_gate"],
          ["fine_pearson", "fine_rmse", "broad_rmse", "false_resolution_rate",
           "false_abstention_rate", "eff_n_abs_error", "rare_sensitivity"]].round(3).to_string())
    print("\n=== PROSPECTIVE GATES ===")
    print(gate_df.to_string(index=False))
    print(f"\nGATES PASSED: {n_pass}/{len(gate_df)}  | 2-tissue replication: UNMET (no 2nd atlas)")
    print(f"Wrote Stage-1 tables -> {OUT}/")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
