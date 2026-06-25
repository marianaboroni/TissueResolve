#!/usr/bin/env python
"""STAGE 4 — second-tissue (HLCA lung) generalization validation.

Tests whether the breast findings generalize to a distinct tissue, using the
donor-disjoint HLCA-core subset (raw counts, broad+fine labels, ~40 donors):

  Part A — identifiability: per within-family pair, signature collinearity +
           shared-lineage fraction + 2-signature mixture-recovery MAE (the
           deconvolution-relevant ceiling). Does shared-lineage dominate in lung?
  Part B — soft-gating prospective gates (held-out-donor pseudobulk, 5 donor
           splits × 3 scenarios): flat / hard / ungated / soft. Does soft gating
           replicate (vs ungated non-inferior; vs hard materially better)?

Donor- AND seed-disjoint calibration {0,1} / test {2,3,4} (rules 4,10,11).
Outputs: benchmarks/outputs/lung_identifiability_family_summary.tsv,
lung_softgating_modes.tsv, lung_softgating_gates.tsv.
Usage:  python 02_validate_generalization.py --run-real-data
"""
from __future__ import annotations

import argparse
import sys
import warnings
from itertools import combinations
from pathlib import Path

import numpy as np
import pandas as pd
from scipy.optimize import nnls

REPO = Path(__file__).resolve().parents[3]
sys.path.insert(0, str(REPO / "examples" / "real_breast_cancer" / "scripts"))
sys.path.insert(0, str(REPO))
import _harness as H  # noqa: E402
from benchmarks.shared import synthetic_holdout as SH  # noqa: E402
from benchmarks.shared import metrics as M  # noqa: E402
from benchmarks.shared import stats as ST  # noqa: E402
from tissueresolve.experimental.soft_hierarchy import (  # noqa: E402
    apply_partial_confidence_gating, compute_reference_confidence_features, fit_confidence_model)

SUBSET = Path(__file__).resolve().parents[1] / "data" / "derived" / "hlca_subset.h5ad"
HMAP = Path(__file__).resolve().parents[1] / "data" / "derived" / "hlca_hierarchy.tsv"
OUT = REPO / "benchmarks" / "outputs"
CAL, TEST = [0, 1], [2, 3, 4]
SCEN_SEED = {"imbalanced": 22, "similar_subtypes": 44, "rare": 33}
N_SAMPLES, CELLS, MIN_CELLS = 12, 600, 20
PRESENT, ABSENT, TOL = 0.01, 0.005, 0.05
RATIOS, N_NOISE, DEPTH = [0.25, 0.5, 0.75], 10, 5e5


def _logcpm(X):
    import scipy.sparse as sp
    X = X.toarray() if sp.issparse(X) else np.asarray(X, float)
    return np.log1p(X / np.clip(X.sum(1, keepdims=True), 1, None) * 1e4)


def part_A_identifiability(adata, mapping):
    """Per within-family pair: signature collinearity + shared-lineage + mixture recovery."""
    import scipy.sparse as sp
    ct = adata.obs["cell_type_fine"].astype(str).to_numpy()
    Xcpm_genes = None
    fam_members = {}
    for c in sorted(set(ct)):
        fam_members.setdefault(mapping.get(c, c), []).append(c)
    fam_members = {f: m for f, m in fam_members.items() if len(m) >= 2}
    # linear CPM signatures per subtype
    def sig(stype):
        ii = np.where(ct == stype)[0]
        Xi = adata.X[ii]; Xi = Xi.toarray() if sp.issparse(Xi) else np.asarray(Xi, float)
        cpm = Xi / np.clip(Xi.sum(1, keepdims=True), 1, None) * 1e6
        return cpm.mean(0)
    rng = np.random.default_rng(0); rows = []
    for fam, members in fam_members.items():
        sigs = {m: sig(m) for m in members}
        for a, b in combinations(members, 2):
            Sa, Sb = sigs[a], sigs[b]; Smat = np.vstack([Sa, Sb]).T
            Bsh = 0.5 * (Sa + Sb)
            shared = float(np.linalg.norm(Bsh) / (0.5 * (np.linalg.norm(Sa) + np.linalg.norm(Sb)) + 1e-12))
            errs = []
            for r in RATIOS:
                mu = Smat @ np.array([r, 1 - r])
                for _ in range(N_NOISE):
                    y = rng.poisson(np.clip(mu / mu.sum() * DEPTH, 0, None)).astype(float)
                    th, _ = nnls(Smat, y); s = th.sum(); th = th / s if s > 0 else np.array([r, 1 - r])
                    errs.append(abs(th[0] - r))
            rows.append({"family": fam, "subtype_a": a, "subtype_b": b,
                         "signature_pearson": round(float(np.corrcoef(Sa, Sb)[0, 1]), 4),
                         "shared_lineage_fraction": round(shared, 4),
                         "ratio_recovery_mae": round(float(np.mean(errs)), 4),
                         "condition_number": round(float(np.linalg.cond(Smat)), 2)})
    pairs = pd.DataFrame(rows)
    fam = pairs.groupby("family").agg(
        n_pairs=("signature_pearson", "size"),
        median_signature_pearson=("signature_pearson", "median"),
        median_shared_lineage_fraction=("shared_lineage_fraction", "median"),
        median_ratio_recovery_mae=("ratio_recovery_mae", "median")).reset_index()
    return pairs, fam


def _prep(adata, mapping, seed):
    from tissueresolve.reference.hierarchy import build_cell_type_hierarchy
    ref_d, qry = SH.split_donors(adata, "donor_id", ref_frac=0.5, seed=seed)
    rmask = adata.obs["donor_id"].astype(str).isin(set(ref_d)).to_numpy()
    with warnings.catch_warnings():
        warnings.simplefilter("ignore")
        ref = H.prepare_reference(adata[rmask].copy(), cell_type_col="cell_type_fine",
                                  min_cells=MIN_CELLS, estimate_overdispersion=True).reference
    ref_types = [str(c) for c in ref.cell_types]
    mp = build_cell_type_hierarchy(ref_types, mapping)
    return ref, ref_types, mp, qry


def _scenario(adata, ref_types, mp, qry, scen, seed):
    kw = {}
    if scen == "similar_subtypes":
        fam = max({mp.get(x, x) for x in ref_types},
                  key=lambda f: sum(mp.get(x, x) == f for x in ref_types))
        kw = {"similar_pair": tuple([t for t in ref_types if mp.get(t) == fam][:2])}
    elif scen == "rare":
        kw = {"rare_type": ref_types[np.argmin([1])] if False else ref_types[-1], "rare_level": 0.01}
    tgt = SH.build_target_proportions(ref_types, N_SAMPLES, scen, seed=SCEN_SEED[scen], **kw)
    return SH.realize_pseudobulk(adata, tgt, celltype_col="cell_type_fine", donor_col="donor_id",
                                 query_donors=qry, seed=seed, cells_per_sample=CELLS)


def _metrics(truth, pred, mp, mode):
    cols = [c for c in truth.columns if not str(c).startswith("unresolved_")]
    t = truth[cols]; p = pred.reindex(index=t.index, columns=cols).fillna(0.0)
    acc = M.accuracy_metrics(t, p)
    tfam, pfam = M.aggregate_to_families(t, mp), M.aggregate_to_families(pred, mp)
    bacc = M.accuracy_metrics(tfam, pfam)
    fam_of = {c: mp.get(c, c) for c in cols}
    ct_, cp_ = t * 0.0, p * 0.0
    for fam in set(fam_of.values()):
        mem = [c for c in cols if fam_of[c] == fam]
        st_, sp_ = t[mem].sum(1), p[mem].sum(1)
        for c in mem:
            ct_[c] = np.where(st_ > 0, t[c] / st_.replace(0, np.nan), 0)
            cp_[c] = np.where(sp_ > 0, p[c] / sp_.replace(0, np.nan), 0)
    cond = M.accuracy_metrics(ct_.fillna(0), cp_.fillna(0))
    tv, pv = t.to_numpy(float), p.to_numpy(float)
    absent = tv <= ABSENT; rare = (tv >= 0.005) & (tv <= 0.05)
    cx, cxt = M.composition_complexity(p), M.composition_complexity(t)
    ucols = [c for c in pred.columns if str(c).startswith("unresolved_")]
    umass = float(pred[ucols].to_numpy().sum() / max(pred.to_numpy().sum(), 1e-9)) if ucols else 0.0
    row = {"broad_rmse": bacc["rmse"], "fine_rmse": acc["rmse"], "fine_pearson": acc["pearson"],
           "cond_pearson": cond["pearson"],
           "false_resolution_rate": float(pv[absent].sum() / max(pv.sum(), 1e-9)),
           "rare_sensitivity": float((pv[rare] > ABSENT).mean()) if rare.any() else float("nan"),
           "eff_n_abs_error": float(abs(cx["effective_n_populations"].mean() - cxt["effective_n_populations"].mean())),
           "unresolved_mass_frac": umass, "mass_error": float(abs(pred.sum(1) - truth.sum(1)).max())}
    if mode in ("hard_gate", "soft_gate"):
        fam_of2 = {c: mp.get(c, c) for c in cols}
        num = den = 0.0
        ucol = {c[len("unresolved_"):]: c for c in pred.columns if str(c).startswith("unresolved_")}
        for fam in set(fam_of2.values()):
            mem = [c for c in cols if fam_of2[c] == fam]
            should_resolve = (truth[mem] > PRESENT).sum(1) >= 2
            if should_resolve.any() and fam in ucol:
                num += float(pred[ucol[fam]][should_resolve].sum()); den += float(truth[mem].sum(1)[should_resolve].sum())
        row["false_abstention_rate"] = float(num / den) if den > 1e-9 else 0.0
    else:
        row["false_abstention_rate"] = 0.0
    return row


def part_B_softgating(adata, mapping):
    from tissueresolve.api import deconv_bulk
    from tissueresolve.reference.hierarchy import combine_family_and_conditional_estimates
    # calibrate confidence on CAL seeds
    feat = []
    for seed in CAL:
        ref, ref_types, mp, qry = _prep(adata, mapping, seed)
        rfeat = compute_reference_confidence_features(ref, mp)
        for scen in SCEN_SEED:
            ds = _scenario(adata, ref_types, mp, qry, scen, seed); truth = ds.true_mrna_proportions
            with warnings.catch_warnings():
                warnings.simplefilter("ignore")
                h = deconv_bulk(ds.counts, ref, solver="auto", resolution_mode="hierarchical", hierarchy_mapping=mp)
            pre = combine_family_and_conditional_estimates(h.estimates.family_proportions,
                                                           h.estimates.conditional_proportions, mp)
            for samp in truth.index:
                for st in truth.columns:
                    if st not in rfeat.index: continue
                    err = abs(float(pre.at[samp, st]) - float(truth.at[samp, st])) if st in pre.columns else float(truth.at[samp, st])
                    feat.append({"subtype": st, "sample": samp, "ref_evidence": rfeat.at[st, "ref_evidence"], "correct": int(err < TOL)})
    fdf = pd.DataFrame(feat).dropna(subset=["ref_evidence"])
    cmodel = fit_confidence_model(fdf.set_index(["sample", "subtype"])[["ref_evidence"]],
                                  fdf.set_index(["sample", "subtype"])["correct"].astype(float),
                                  kind="logistic", score_col="ref_evidence")
    print(f"  calibrated confidence on {len(fdf)} cal pairs", flush=True)
    rows = []
    for seed in TEST:
        ref, ref_types, mp, qry = _prep(adata, mapping, seed)
        conf_feat = compute_reference_confidence_features(ref, mp)[["ref_evidence"]]
        for scen in SCEN_SEED:
            ds = _scenario(adata, ref_types, mp, qry, scen, seed); truth = ds.true_mrna_proportions
            with warnings.catch_warnings():
                warnings.simplefilter("ignore")
                flat = deconv_bulk(ds.counts, ref, solver="nnls", resolution_mode="none").deconv.proportions
                h = deconv_bulk(ds.counts, ref, solver="auto", resolution_mode="hierarchical", hierarchy_mapping=mp)
            est = h.estimates
            pre = combine_family_and_conditional_estimates(est.family_proportions, est.conditional_proportions, mp)
            soft = apply_partial_confidence_gating(est.family_proportions, pre, mp,
                                                   confidence=cmodel.predict(conf_feat)).combined
            for mode, pred in [("flat", flat), ("hard_gate", est.combined_fine), ("ungated", pre), ("soft_gate", soft)]:
                rows.append({"seed": seed, "scenario": scen, "mode": mode, **_metrics(truth, pred, mp, mode)})
        print(f"  seed{seed} done", flush=True)
    return pd.DataFrame(rows)


def main(argv=None) -> int:
    ap = argparse.ArgumentParser(); ap.add_argument("--run-real-data", action="store_true")
    args = ap.parse_args(argv)
    if not args.run_real_data:
        print("Refusing without --run-real-data.", file=sys.stderr); return 2
    if not SUBSET.exists():
        print(f"Missing {SUBSET}. Run 01_audit_and_subsample.py first.", file=sys.stderr); return 2
    import anndata as ad
    adata = ad.read_h5ad(SUBSET)
    if "feature_name" in adata.var.columns:
        adata.var_names = adata.var["feature_name"].astype(str).to_numpy(); adata.var_names_make_unique()
    hmap = pd.read_csv(HMAP, sep="\t")
    mapping = dict(zip(hmap["fine_cell_type"].astype(str), hmap["broad_cell_type"].astype(str)))
    print(f"lung subset: {adata.shape}; donors={adata.obs.donor_id.nunique()} "
          f"broad={adata.obs.cell_type_broad.nunique()} fine={adata.obs.cell_type_fine.nunique()}", flush=True)

    print("Part A — identifiability ...", flush=True)
    pairs, famA = part_A_identifiability(adata, mapping)
    pairs.to_csv(OUT / "lung_identifiability_pairs.tsv", sep="\t", index=False)
    famA.to_csv(OUT / "lung_identifiability_family_summary.tsv", sep="\t", index=False)
    print(famA.to_string(index=False), flush=True)

    print("\nPart B — soft-gating gates ...", flush=True)
    comp = part_B_softgating(adata, mapping)
    comp.to_csv(OUT / "lung_softgating_modes.tsv", sep="\t", index=False)
    g = comp.groupby("mode").mean(numeric_only=True)
    soft, ung, hard = g.loc["soft_gate"], g.loc["ungated"], g.loc["hard_gate"]
    rel = lambda a, b: (a - b) / abs(b) if b else float("inf")
    gates = [
        ("vs ungated: broad RMSE ≤+2%", soft["broad_rmse"], ung["broad_rmse"], rel(soft["broad_rmse"], ung["broad_rmse"]) <= 0.02),
        ("vs ungated: fine RMSE ≤+2%", soft["fine_rmse"], ung["fine_rmse"], rel(soft["fine_rmse"], ung["fine_rmse"]) <= 0.02),
        ("vs ungated: false-resolution ≤+10%", soft["false_resolution_rate"], ung["false_resolution_rate"], rel(soft["false_resolution_rate"], ung["false_resolution_rate"]) <= 0.10),
        ("vs ungated: mass conserved", soft["mass_error"], 1e-6, soft["mass_error"] < 1e-6),
        ("vs hard: false-abstention < 0.5×", soft["false_abstention_rate"], hard["false_abstention_rate"], soft["false_abstention_rate"] < 0.5 * hard["false_abstention_rate"]),
        ("vs hard: rare sensitivity > 1.5×", soft["rare_sensitivity"], hard["rare_sensitivity"], soft["rare_sensitivity"] > 1.5 * hard["rare_sensitivity"]),
        ("vs hard: eff-N error < 0.6×", soft["eff_n_abs_error"], hard["eff_n_abs_error"], soft["eff_n_abs_error"] < 0.6 * hard["eff_n_abs_error"]),
        ("vs hard: fine RMSE lower", soft["fine_rmse"], hard["fine_rmse"], soft["fine_rmse"] < hard["fine_rmse"]),
    ]
    gate_df = pd.DataFrame([{"gate": n, "soft": round(float(a), 4), "ref": round(float(b), 6), "pass": bool(p)} for n, a, b, p in gates])
    gate_df.to_csv(OUT / "lung_softgating_gates.tsv", sep="\t", index=False)
    print("\n=== lung soft-gating mode means ===")
    print(g.loc[["flat", "ungated", "soft_gate", "hard_gate"],
          ["fine_pearson", "fine_rmse", "broad_rmse", "false_resolution_rate",
           "false_abstention_rate", "eff_n_abs_error", "rare_sensitivity"]].round(3).to_string())
    print("\n=== lung soft-gating PROSPECTIVE GATES ===")
    print(gate_df.to_string(index=False))
    print(f"\nGATES PASSED: {int(gate_df['pass'].sum())}/{len(gate_df)} | "
          f"median shared-lineage fraction (lung): {pairs['shared_lineage_fraction'].median():.3f}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
