#!/usr/bin/env python
"""PHASE 2B Stage 2 — soft joint broad–fine hierarchy benchmark (controlled).

Tests the hypothesis that the hierarchical-vs-flat gap is the *sequential*
factorisation.  All methods run on IDENTICAL preprocessing (shared-gene CPM-log1p
S and y) so the comparison isolates the factorisation:

  flat            : NNLS(S, y) renormalised
  sequential      : broad NNLS(B,y) → π ; within-family conditional NNLS → combine
  joint_alt       : fit_alternating_soft_hierarchy (soft consistency π≈Aθ)
  joint_soft      : joint_alt + partial confidence gating (orthogonal layer)
  oracle_broad    : combine TRUE family mass × sequential conditional

Outputs: benchmarks/outputs/phase2b_joint_hierarchy_metrics.tsv,
phase2b_error_decomposition.tsv, phase2b_ablation_metrics.tsv,
phase2_promotion_gates.tsv.  Experimental; promotion blocked (1 tissue).

Usage:  python benchmarks/diagnostics/phase2b_stage2_joint.py --run-real-data
"""
from __future__ import annotations

import argparse
import sys
import warnings
from pathlib import Path

import numpy as np
import pandas as pd
from scipy.optimize import nnls

REPO = Path(__file__).resolve().parents[2]
sys.path.insert(0, str(REPO / "examples" / "real_breast_cancer" / "scripts"))
sys.path.insert(0, str(REPO))
import _harness as H  # noqa: E402
from benchmarks.shared import synthetic_holdout as SH  # noqa: E402
from benchmarks.shared import metrics as M  # noqa: E402
from benchmarks.shared import stats as ST  # noqa: E402
from tissueresolve.experimental.soft_hierarchy import (  # noqa: E402
    fit_alternating_soft_hierarchy, apply_partial_confidence_gating,
    compute_reference_confidence_features, fit_confidence_model,
    compute_error_decomposition)

OUT = REPO / "benchmarks" / "outputs"
HIER = REPO / "examples" / "real_breast_cancer" / "config" / "breast_cancer_cell_type_hierarchy.tsv"
TEST_SEEDS = [3, 4, 5, 6, 7, 8]
SCENARIOS = ["imbalanced", "similar_subtypes", "rare"]
SCEN_SEED = {"imbalanced": 22, "similar_subtypes": 44, "rare": 33}
N_SAMPLES, CELLS, MIN_CELLS = 12, 600, 30


def _simplex(M, b):
    t, _ = nnls(M, b); s = t.sum()
    return t / s if s > 0 else np.full(M.shape[1], 1.0 / M.shape[1])


def _prep(adata, raw_map, seed):
    from tissueresolve.reference.hierarchy import build_cell_type_hierarchy
    ref_d, qry = SH.split_donors(adata, "donor_id", ref_frac=0.5, seed=seed)
    mapped = adata.obs["cell_type"].astype(str).isin(set(raw_map)).to_numpy()
    rmask = adata.obs["donor_id"].astype(str).isin(set(ref_d)).to_numpy() & mapped
    with warnings.catch_warnings():
        warnings.simplefilter("ignore")
        ref = H.prepare_reference(adata[rmask].copy(), min_cells=MIN_CELLS,
                                  estimate_overdispersion=True).reference
    ref_types = [str(c) for c in ref.cell_types]
    return ref, ref_types, build_cell_type_hierarchy(ref_types, raw_map), qry


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
    return SH.realize_pseudobulk(adata, tgt, celltype_col="cell_type", donor_col="donor_id",
                                 query_donors=qry, seed=seed, cells_per_sample=CELLS)


def _build_Sy(ref, counts):
    """Shared-gene CPM-log1p fine signatures S (G×K) and samples y (n×G)."""
    rg = [str(g) for g in ref.gene_names]; gi = {g: i for i, g in enumerate(rg)}
    shared = [g for g in map(str, counts.index) if g in gi]
    R = np.asarray(ref.as_R_cpm())[:, [gi[g] for g in shared]]   # K×G
    S = np.log1p(R).T                                            # G×K
    Y = counts.loc[shared].to_numpy(float).T                    # n×G(counts)
    Y = np.log1p(Y / np.clip(Y.sum(1, keepdims=True), 1, None) * 1e6)
    return S, Y, shared


def _sequential(S, Y, ref_types, mapping):
    """Sequential factorisation on identical inputs (broad → conditional → combine)."""
    fams = sorted({mapping.get(t, t) for t in ref_types})
    fidx = {f: i for i, f in enumerate(fams)}
    A = np.zeros((len(fams), len(ref_types)))
    for j, k in enumerate(ref_types):
        A[fidx[mapping.get(k, k)], j] = 1.0
    B = S @ (A.T / np.clip(A.sum(1, keepdims=True), 1, None).T)
    fine = np.zeros((Y.shape[0], len(ref_types)))
    broad = np.zeros((Y.shape[0], len(fams)))
    for s in range(Y.shape[0]):
        pi = _simplex(B, Y[s]); broad[s] = pi
        # conditional within each family from a fine NNLS, then × family mass
        th_flat = _simplex(S, Y[s])
        for f, fi in fidx.items():
            mem = [j for j, k in enumerate(ref_types) if mapping.get(k, k) == f]
            sub = th_flat[mem]; ssum = sub.sum()
            cond = sub / ssum if ssum > 0 else np.full(len(mem), 1 / len(mem))
            for jj, j in enumerate(mem):
                fine[s, j] = pi[fi] * cond[jj]
    return (pd.DataFrame(broad, columns=fams), pd.DataFrame(fine, columns=ref_types))


def _metrics(truth, pred_df, mapping):
    cols = list(truth.columns)
    p = pred_df.reindex(index=truth.index, columns=cols).fillna(0.0)
    acc = M.accuracy_metrics(truth, p); comp = M.compositional_metrics(truth, p)
    tfam = M.aggregate_to_families(truth, mapping); pfam = M.aggregate_to_families(p, mapping)
    bacc = M.accuracy_metrics(tfam, pfam)
    fam_of = {c: mapping.get(c, c) for c in cols}
    ct, cp = truth * 0.0, p * 0.0
    for fam in set(fam_of.values()):
        mem = [c for c in cols if fam_of[c] == fam]
        st_, sp_ = truth[mem].sum(1), p[mem].sum(1)
        for c in mem:
            ct[c] = np.where(st_ > 0, truth[c] / st_.replace(0, np.nan), 0)
            cp[c] = np.where(sp_ > 0, p[c] / sp_.replace(0, np.nan), 0)
    cond = M.accuracy_metrics(ct.fillna(0), cp.fillna(0))
    tv, pv = truth.to_numpy(float), p.to_numpy(float)
    absent = tv <= 0.005; rare = (tv >= 0.005) & (tv <= 0.05)
    cx = M.composition_complexity(p); cxt = M.composition_complexity(truth)
    return {"broad_pearson": bacc["pearson"], "broad_rmse": bacc["rmse"],
            "fine_pearson": acc["pearson"], "fine_rmse": acc["rmse"],
            "cond_pearson": cond["pearson"], "cond_rmse": cond["rmse"],
            "jsd": comp["jsd_mean"], "aitchison": comp["aitchison_mean"],
            "rare_sensitivity": float((pv[rare] > 0.005).mean()) if rare.any() else float("nan"),
            "false_resolution_rate": float(pv[absent].sum() / max(pv.sum(), 1e-9)),
            "eff_n_abs_error": float(abs(cx["effective_n_populations"].mean() - cxt["effective_n_populations"].mean())),
            "mass_error": float(abs(p.sum(1) - truth.sum(1)).max())}


def main(argv=None) -> int:
    ap = argparse.ArgumentParser(); ap.add_argument("--run-real-data", action="store_true")
    args = ap.parse_args(argv)
    if not H.real_data_enabled(args.run_real_data):
        print("Refusing without --run-real-data.", file=sys.stderr); return 2
    import anndata as ad
    from tissueresolve.reference.hierarchy import load_hierarchy_mapping
    adata = ad.read_h5ad(H.REFERENCE_H5AD)
    if "feature_name" in adata.var.columns:
        fn = adata.var["feature_name"].astype(str)
        if fn.nunique() == adata.n_vars:
            adata.var_names = fn.to_numpy()
    raw_map = load_hierarchy_mapping(HIER)

    rows, dec_rows, abl_rows = [], [], []
    LAMBDAS = [0.0, 0.1, 1.0, 5.0]
    for seed in TEST_SEEDS:
        ref, ref_types, mapping, qry = _prep(adata, raw_map, seed)
        conf_feat = compute_reference_confidence_features(ref, mapping)[["ref_evidence"]]
        for scen in SCENARIOS:
            ds = _scenario(adata, ref_types, mapping, qry, scen, seed)
            truth = ds.true_mrna_proportions.reindex(columns=ref_types).fillna(0.0)
            S, Y, shared = _build_Sy(ref, ds.counts)
            ids = list(truth.index)
            # flat
            flat = pd.DataFrame([_simplex(S, Y[s]) for s in range(Y.shape[0])],
                                index=ids, columns=ref_types)
            # sequential (same inputs)
            seq_b, seq_f = _sequential(S, Y, ref_types, mapping)
            seq_b.index = ids; seq_f.index = ids
            # joint (alternating, lambda_h=1)
            jr = fit_alternating_soft_hierarchy(Y, S, mapping, ref_types, lambda_h=1.0,
                                                lambda_2=0.0, sample_ids=ids)
            joint = jr.reconciled_fine_proportions
            # oracle broad × sequential conditional
            true_fam = M.aggregate_to_families(truth, mapping)
            from tissueresolve.reference.hierarchy import combine_family_and_conditional_estimates
            # conditional from sequential fine
            cond = seq_f.copy()
            for fam in seq_b.columns:
                mem = [c for c in ref_types if mapping.get(c, c) == fam]
                s_ = seq_f[mem].sum(1)
                for c in mem:
                    cond[c] = np.where(s_ > 0, seq_f[c] / s_.replace(0, np.nan), 0)
            oracle = combine_family_and_conditional_estimates(
                true_fam.reindex(columns=seq_b.columns).fillna(0.0), cond.fillna(0), mapping)
            for mode, pred in [("flat", flat), ("sequential", seq_f), ("joint_alt", joint),
                               ("oracle_broad", oracle)]:
                rows.append({"seed": seed, "scenario": scen, "mode": mode,
                             **_metrics(truth, pred, mapping)})
            # error decomposition: sequential vs joint
            dec_seq = compute_error_decomposition(truth, seq_b, seq_f, mapping)
            dec_joint = compute_error_decomposition(truth, jr.broad_proportions, joint, mapping)
            dec_rows.append({"seed": seed, "scenario": scen, "mode": "sequential", **dec_seq})
            dec_rows.append({"seed": seed, "scenario": scen, "mode": "joint_alt", **dec_joint})
        # lambda_h ablation (seed-level, scenario=imbalanced)
        ds = _scenario(adata, ref_types, mapping, qry, "imbalanced", seed)
        truth = ds.true_mrna_proportions.reindex(columns=ref_types).fillna(0.0)
        S, Y, _ = _build_Sy(ref, ds.counts); ids = list(truth.index)
        for lh in LAMBDAS:
            jr = fit_alternating_soft_hierarchy(Y, S, mapping, ref_types, lambda_h=lh, sample_ids=ids)
            m = _metrics(truth, jr.reconciled_fine_proportions, mapping)
            abl_rows.append({"seed": seed, "lambda_h": lh, "fine_rmse": m["fine_rmse"],
                             "fine_pearson": m["fine_pearson"], "cond_rmse": m["cond_rmse"],
                             "consistency": float(jr.hierarchy_consistency_error.mean())})
        print(f"  seed{seed} done", flush=True)

    df = pd.DataFrame(rows)
    OUT.mkdir(parents=True, exist_ok=True)
    with open(OUT / "phase2b_joint_hierarchy_metrics.tsv", "w") as fh:
        fh.write("# Phase 2B Stage-2 controlled comparison (identical preprocessing). Experimental; 1 tissue.\n")
        df.to_csv(fh, sep="\t", index=False)
    pd.DataFrame(dec_rows).to_csv(OUT / "phase2b_error_decomposition.tsv", sep="\t", index=False)
    pd.DataFrame(abl_rows).to_csv(OUT / "phase2b_ablation_metrics.tsv", sep="\t", index=False)

    g = df.groupby("mode").mean(numeric_only=True)
    # Stage-2 promotion gates: joint vs sequential + vs flat
    seq, joint, flat = g.loc["sequential"], g.loc["joint_alt"], g.loc["flat"]
    rel = lambda a, b: (a - b) / abs(b) if b else float("inf")
    gates = [
        ("vs sequential: fine RMSE ≥5% better", joint["fine_rmse"], seq["fine_rmse"],
         rel(joint["fine_rmse"], seq["fine_rmse"]) <= -0.05),
        ("vs sequential: conditional RMSE improves", joint["cond_rmse"], seq["cond_rmse"],
         joint["cond_rmse"] < seq["cond_rmse"]),
        ("vs sequential: broad RMSE ≤ +2%", joint["broad_rmse"], seq["broad_rmse"],
         rel(joint["broad_rmse"], seq["broad_rmse"]) <= 0.02),
        ("vs sequential: rare sensitivity non-inferior", joint["rare_sensitivity"], seq["rare_sensitivity"],
         joint["rare_sensitivity"] >= seq["rare_sensitivity"] * 0.95),
        ("vs flat: broad RMSE non-inferior (≤+2%)", joint["broad_rmse"], flat["broad_rmse"],
         rel(joint["broad_rmse"], flat["broad_rmse"]) <= 0.02),
        ("vs flat: fine RMSE gap reduced vs sequential", joint["fine_rmse"], flat["fine_rmse"],
         (joint["fine_rmse"] - flat["fine_rmse"]) < (seq["fine_rmse"] - flat["fine_rmse"])),
        ("mass conserved <1e-6", joint["mass_error"], 1e-6, joint["mass_error"] < 1e-6),
    ]
    gate_df = pd.DataFrame([{"gate": n, "joint": round(float(a), 4), "ref": round(float(b), 6),
                             "pass": bool(p)} for n, a, b, p in gates])
    with open(OUT / "phase2_promotion_gates.tsv", "w") as fh:
        fh.write("# Phase 2B Stage-2 promotion gates (joint vs sequential/flat). Experimental; promotion blocked (1 tissue).\n")
        gate_df.to_csv(fh, sep="\t", index=False)
    print("\n=== Stage-2 mode means ===")
    print(g.loc[["flat", "sequential", "joint_alt", "oracle_broad"],
          ["fine_pearson", "fine_rmse", "broad_rmse", "cond_rmse", "rare_sensitivity",
           "eff_n_abs_error"]].round(3).to_string())
    print("\n=== error decomposition (mean) ===")
    print(pd.DataFrame(dec_rows).groupby("mode")[["broad_error", "conditional_fine_error",
          "consistency_error", "total_fine_error"]].mean().round(4).to_string())
    print("\n=== lambda_h ablation (mean) ===")
    print(pd.DataFrame(abl_rows).groupby("lambda_h")[["fine_rmse", "fine_pearson",
          "cond_rmse", "consistency"]].mean().round(4).to_string())
    print("\n=== STAGE-2 PROMOTION GATES ===")
    print(gate_df.to_string(index=False))
    print(f"\nGATES PASSED: {int(gate_df['pass'].sum())}/{len(gate_df)} | promotion blocked (1 tissue)")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
