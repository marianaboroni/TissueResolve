#!/usr/bin/env python
"""STAGE 1 — subpopulation identifiability ceiling (diagnosis only).

Estimates the UPPER BOUND of within-family subtype discriminability using simple
classifiers under **donor-held-out** cross-validation and different gene panels,
BEFORE attempting any new solver/contrast.  If subtypes are not separable even by
an oracle classifier on donor-stable, query-detectable genes, no deconvolution
method can resolve them — preserve unresolved mass and report broad-only.

Panels: all reference genes; query-detectable (genes seen in TCGA-TNBC bulk);
donor-stable (low cross-donor DE variance).  Classifiers: logistic regression,
LDA, nearest-centroid (correlation).  Donor-held-out via GroupKFold by donor_id.

Outputs:
  benchmarks/outputs/subtype_identifiability_pairs.tsv
  benchmarks/outputs/subtype_identifiability_family_summary.tsv
  benchmarks/outputs/subtype_query_detectability.tsv
  benchmarks/outputs/subtype_donor_stability.tsv
  figures/identifiability/*.png (+ .data.tsv + caption)

Usage:  python benchmarks/diagnostics/identifiability_ceiling.py --run-real-data
"""
from __future__ import annotations

import argparse
import sys
import warnings
from itertools import combinations
from pathlib import Path

import numpy as np
import pandas as pd

REPO = Path(__file__).resolve().parents[2]
sys.path.insert(0, str(REPO / "examples" / "real_breast_cancer" / "scripts"))
sys.path.insert(0, str(REPO))
import _harness as H  # noqa: E402

OUT = REPO / "benchmarks" / "outputs"
FIG = OUT / "figures" / "identifiability"
HIER = REPO / "examples" / "real_breast_cancer" / "config" / "breast_cancer_cell_type_hierarchy.tsv"
TCGA = H.DATA_DIR / "bulk_tcga_tnbc" / "tcga_tnbc_counts.tsv"
MAX_CELLS_PER_TYPE = 400          # subsample for tractable CV
MIN_CELLS_PER_DONOR = 10
MIN_DONORS = 3
TOP_GENES = 300
SEED = 0
THRESH = {"identifiable": 0.90, "conditional": 0.75, "ambiguous": 0.60}


def _logcpm(X):
    X = np.asarray(X.todense()) if hasattr(X, "todense") else np.asarray(X, float)
    lib = np.clip(X.sum(1, keepdims=True), 1, None)
    return np.log1p(X / lib * 1e4)


def _auroc(y, score):
    from sklearn.metrics import roc_auc_score, average_precision_score, balanced_accuracy_score
    try:
        return (roc_auc_score(y, score), average_precision_score(y, score))
    except Exception:
        return (float("nan"), float("nan"))


def _donor_holdout_auroc(Xl, y, donors, genes_idx, clf_kind):
    """Mean donor-held-out AUROC/AUPRC/balanced-acc + donor-to-donor variance."""
    from sklearn.model_selection import GroupKFold
    from sklearn.linear_model import LogisticRegression
    from sklearn.discriminant_analysis import LinearDiscriminantAnalysis
    from sklearn.neighbors import NearestCentroid
    from sklearn.metrics import balanced_accuracy_score
    uniq = np.unique(donors)
    if len(uniq) < MIN_DONORS:
        return None
    k = min(5, len(uniq))
    gkf = GroupKFold(n_splits=k)
    X = Xl[:, genes_idx]
    aurocs, auprcs, baccs = [], [], []
    for tr, te in gkf.split(X, y, groups=donors):
        if len(np.unique(y[tr])) < 2 or len(np.unique(y[te])) < 2:
            continue
        if clf_kind == "logreg":
            m = LogisticRegression(max_iter=500, C=1.0)
            m.fit(X[tr], y[tr]); sc = m.predict_proba(X[te])[:, 1]; pred = (sc >= 0.5).astype(int)
        elif clf_kind == "lda":
            m = LinearDiscriminantAnalysis()
            m.fit(X[tr], y[tr])
            sc = m.predict_proba(X[te])[:, 1] if hasattr(m, "predict_proba") else m.decision_function(X[te])
            pred = m.predict(X[te])
        else:  # nearest-centroid (correlation-like)
            m = NearestCentroid(metric="euclidean")
            m.fit(X[tr], y[tr]); pred = m.predict(X[te])
            # score = signed distance proxy
            c0, c1 = m.centroids_
            sc = -(np.linalg.norm(X[te] - c1, axis=1) - np.linalg.norm(X[te] - c0, axis=1))
        a, p = _auroc(y[te], sc)
        aurocs.append(a); auprcs.append(p)
        baccs.append(balanced_accuracy_score(y[te], pred))
    if not aurocs:
        return None
    return {"auroc": float(np.nanmean(aurocs)), "auprc": float(np.nanmean(auprcs)),
            "balanced_acc": float(np.nanmean(baccs)),
            "donor_auroc_sd": float(np.nanstd(aurocs)), "n_folds": len(aurocs)}


def main(argv=None) -> int:
    ap = argparse.ArgumentParser(); ap.add_argument("--run-real-data", action="store_true")
    args = ap.parse_args(argv)
    if not H.real_data_enabled(args.run_real_data):
        print("Refusing without --run-real-data.", file=sys.stderr); return 2
    import anndata as ad
    from tissueresolve.reference.hierarchy import load_hierarchy_mapping, build_cell_type_hierarchy

    adata = ad.read_h5ad(H.REFERENCE_H5AD)
    if "feature_name" in adata.var.columns:
        fn = adata.var["feature_name"].astype(str)
        if fn.nunique() == adata.n_vars:
            adata.var_names = fn.to_numpy()
    genes = [str(g) for g in adata.var_names]
    raw_map = load_hierarchy_mapping(HIER)
    ct = adata.obs["cell_type"].astype(str).to_numpy()
    don = adata.obs["donor_id"].astype(str).to_numpy()
    mapped = np.isin(ct, list(raw_map))
    mapping = build_cell_type_hierarchy(sorted(set(ct[mapped])), raw_map)

    # query-detectable genes from TCGA bulk
    query_genes = set(genes)
    if TCGA.exists():
        bulk = pd.read_csv(TCGA, sep="\t", index_col=0, comment="#")
        det = (bulk > 0).mean(axis=1)
        bulk_detected = set(map(str, det.index[det >= 0.5]))
        query_genes = {g for g in genes if g in bulk_detected}
    print(f"query-detectable genes (TCGA ≥50% samples): {len(query_genes)}/{len(genes)}", flush=True)

    Xl = _logcpm(adata.X)                       # log-CPM, cells × genes
    gidx = {g: i for i, g in enumerate(genes)}
    qcols = np.array([gidx[g] for g in genes if g in query_genes])

    # families with ≥2 subtypes
    fam_members = {}
    for c in sorted(set(ct[mapped])):
        fam_members.setdefault(mapping.get(c, c), []).append(c)
    fam_members = {f: m for f, m in fam_members.items() if len(m) >= 2}

    rng = np.random.default_rng(SEED)
    pair_rows, qdet_rows = [], []
    for fam, members in fam_members.items():
        for a, b in combinations(members, 2):
            ia = np.where((ct == a) & mapped)[0]; ib = np.where((ct == b) & mapped)[0]
            if ia.size > MAX_CELLS_PER_TYPE: ia = rng.choice(ia, MAX_CELLS_PER_TYPE, replace=False)
            if ib.size > MAX_CELLS_PER_TYPE: ib = rng.choice(ib, MAX_CELLS_PER_TYPE, replace=False)
            idx = np.concatenate([ia, ib])
            y = np.concatenate([np.zeros(len(ia)), np.ones(len(ib))]).astype(int)
            donors = don[idx]
            # require ≥MIN_DONORS donors with both classes? at least MIN_DONORS donors total
            if len(np.unique(donors)) < MIN_DONORS:
                continue
            Xpair = Xl[idx]
            # top discriminative genes (by abs mean diff) within ALL genes and within query genes
            mu0, mu1 = Xpair[y == 0].mean(0), Xpair[y == 1].mean(0)
            diff = np.abs(mu1 - mu0)
            top_all = np.argsort(diff)[::-1][:TOP_GENES]
            qmask = np.zeros(len(genes), bool); qmask[qcols] = True
            top_q = np.array([g for g in np.argsort(diff)[::-1] if qmask[g]][:TOP_GENES])
            # signature similarity on the pair means
            cos = float(np.dot(mu0, mu1) / (np.linalg.norm(mu0) * np.linalg.norm(mu1) + 1e-12))
            pear = float(np.corrcoef(mu0, mu1)[0, 1])
            cond = float(np.linalg.cond(np.vstack([mu0, mu1]))) if min(len(mu0), 2) >= 2 else float("nan")
            row = {"family": fam, "subtype_a": a, "subtype_b": b,
                   "n_a": int(len(ia)), "n_b": int(len(ib)),
                   "n_donors": int(len(np.unique(donors))),
                   "cosine_similarity": cos, "signature_pearson": pear, "condition_number": cond,
                   "n_query_top_genes": int(len(top_q))}
            for panel, gi in [("all_genes", top_all), ("query_detectable", top_q)]:
                if len(gi) < 5:
                    continue
                for clf in ["logreg", "nearest_centroid"]:
                    r = _donor_holdout_auroc(Xl[idx], y, donors, gi, clf)
                    if r:
                        row[f"{panel}__{clf}__auroc"] = round(r["auroc"], 4)
                        row[f"{panel}__{clf}__auprc"] = round(r["auprc"], 4)
                        row[f"{panel}__{clf}__bacc"] = round(r["balanced_acc"], 4)
                        row[f"{panel}__{clf}__donor_sd"] = round(r["donor_auroc_sd"], 4)
            pair_rows.append(row)
        print(f"  family {fam}: {len(members)} subtypes, {len(list(combinations(members,2)))} pairs", flush=True)

    pairs = pd.DataFrame(pair_rows)
    OUT.mkdir(parents=True, exist_ok=True)
    pairs.to_csv(OUT / "subtype_identifiability_pairs.tsv", sep="\t", index=False)

    # primary metric: query-detectable, logreg, donor-held-out AUROC
    pcol = "query_detectable__logreg__auroc"
    acol = "all_genes__logreg__auroc"

    def _bucket(a):
        if not np.isfinite(a): return "unknown"
        if a >= THRESH["identifiable"]: return "identifiable"
        if a >= THRESH["conditional"]: return "conditional"
        if a >= THRESH["ambiguous"]: return "ambiguous"
        return "not_identifiable"

    pairs["bucket_query_logreg"] = pairs[pcol].map(_bucket) if pcol in pairs else "unknown"
    # family summary
    fam_rows = []
    for fam, g in pairs.groupby("family"):
        med_q = float(g[pcol].median()) if pcol in g else float("nan")
        med_all = float(g[acol].median()) if acol in g else float("nan")
        frac_id = float((g[pcol] >= THRESH["conditional"]).mean()) if pcol in g else float("nan")
        frac_unid = float((g[pcol] < THRESH["ambiguous"]).mean()) if pcol in g else float("nan")
        rec = ("fine" if med_q >= THRESH["identifiable"]
               else "selected_fine" if med_q >= THRESH["conditional"]
               else "broad_only")
        fam_rows.append({"family": fam, "n_pairs": len(g),
                         "median_auroc_query_logreg": round(med_q, 4),
                         "median_auroc_all_logreg": round(med_all, 4),
                         "frac_pairs_conditional_or_better": round(frac_id, 3),
                         "frac_pairs_not_identifiable(<0.60)": round(frac_unid, 3),
                         "recommended_resolution": rec})
    famsum = pd.DataFrame(fam_rows).sort_values("median_auroc_query_logreg")
    famsum.to_csv(OUT / "subtype_identifiability_family_summary.tsv", sep="\t", index=False)

    # query detectability per subtype (fraction of its top markers query-detectable)
    for fam, members in fam_members.items():
        for st in members:
            ii = np.where((ct == st) & mapped)[0]
            if ii.size == 0: continue
            mu = Xl[ii].mean(0); top = np.argsort(mu)[::-1][:50]
            qfrac = float(np.mean([genes[g] in query_genes for g in top]))
            qdet_rows.append({"family": fam, "subtype": st, "n_cells": int(ii.size),
                              "n_donors": int(len(np.unique(don[ii]))),
                              "top50_marker_query_detectable_frac": round(qfrac, 3)})
    pd.DataFrame(qdet_rows).to_csv(OUT / "subtype_query_detectability.tsv", sep="\t", index=False)
    # donor stability: per pair, donor_sd of AUROC (already captured); export compact
    ds_cols = [c for c in pairs.columns if c.endswith("donor_sd")]
    pairs[["family", "subtype_a", "subtype_b"] + ds_cols].to_csv(
        OUT / "subtype_donor_stability.tsv", sep="\t", index=False)

    # ---- STOP/GO GATE ----
    valid = pairs[pcol].dropna() if pcol in pairs else pd.Series(dtype=float)
    frac_unident = float((valid < THRESH["ambiguous"]).mean()) if len(valid) else float("nan")
    gate_force_solver = frac_unident <= 0.5  # if >50% pairs <0.60 → do NOT force with solver

    _figures(pairs, famsum, pcol, acol)
    print("\n=== FAMILY identifiability (donor-held-out, query-detectable, logreg) ===")
    print(famsum.to_string(index=False))
    print(f"\nPairs total: {len(pairs)} | with query-logreg AUROC: {len(valid)}")
    print(f"Fraction of pairs NOT identifiable (AUROC<0.60): {frac_unident:.2f}")
    print(f"\nSTOP/GO: {'most pairs unidentifiable → DO NOT force with a complex solver/contrast; report broad-only' if not gate_force_solver else 'some headroom exists → contrastive signatures MAY help'}")
    print(f"Wrote Stage-1 tables + figures -> {OUT}/")
    return 0


def _figures(pairs, famsum, pcol, acol):
    import matplotlib; matplotlib.use("Agg")
    import matplotlib.pyplot as plt
    FIG.mkdir(parents=True, exist_ok=True)

    def save(fig, name, data, cap):
        fig.tight_layout(); fig.savefig(FIG / f"{name}.png", dpi=150); fig.savefig(FIG / f"{name}.svg")
        plt.close(fig); data.to_csv(FIG / f"{name}.data.tsv", sep="\t")
        (FIG / f"{name}.caption.txt").write_text(cap.strip() + "\n", encoding="utf-8")

    if pcol in pairs:
        v = pairs[pcol].dropna()
        fig, ax = plt.subplots(figsize=(6, 4))
        ax.hist(v, bins=20, color="#2c6fbb", alpha=0.85)
        for t, c in [(0.6, "red"), (0.75, "orange"), (0.9, "green")]:
            ax.axvline(t, color=c, ls="--", lw=1)
        ax.set_xlabel("donor-held-out AUROC (query-detectable, logreg)"); ax.set_ylabel("# subtype pairs")
        ax.set_title("Within-family subtype identifiability ceiling")
        save(fig, "fig_auroc_distribution", v.to_frame("auroc"),
             "Donor-held-out pairwise AUROC on query-detectable genes. Lines: 0.60 "
             "(ambiguous), 0.75 (conditional), 0.90 (identifiable).")
        # reference-only vs query-detectable
        if acol in pairs:
            d = pairs[[acol, pcol]].dropna()
            fig, ax = plt.subplots(figsize=(5, 5))
            ax.scatter(d[acol], d[pcol], alpha=0.7)
            ax.plot([0.4, 1], [0.4, 1], "k--", lw=0.8)
            ax.set_xlabel("AUROC (all genes)"); ax.set_ylabel("AUROC (query-detectable)")
            ax.set_title("Reference-only vs query-detectable identifiability")
            save(fig, "fig_ref_vs_query", d, "Each point = a within-family subtype pair.")
    # recommended resolution by family
    fig, ax = plt.subplots(figsize=(7, 4))
    colors = {"fine": "#2ca02c", "selected_fine": "#ff7f0e", "broad_only": "#888"}
    ax.barh(range(len(famsum)), famsum["median_auroc_query_logreg"],
            color=[colors[r] for r in famsum["recommended_resolution"]])
    ax.set_yticks(range(len(famsum))); ax.set_yticklabels(famsum["family"])
    ax.axvline(0.75, color="orange", ls="--"); ax.axvline(0.90, color="green", ls="--")
    ax.set_xlabel("median donor-held-out AUROC (query-detectable)")
    ax.set_title("Recommended resolution by family")
    save(fig, "fig_recommended_resolution", famsum.set_index("family"),
         "Median within-family pair AUROC → recommended resolution (green=fine, "
         "orange=selected-fine, grey=broad-only).")
    # signature similarity vs held-out accuracy
    if pcol in pairs:
        d = pairs[["signature_pearson", pcol]].dropna()
        fig, ax = plt.subplots(figsize=(5.5, 4))
        ax.scatter(d["signature_pearson"], d[pcol], alpha=0.7)
        ax.set_xlabel("subtype signature Pearson (similarity)"); ax.set_ylabel("donor-held-out AUROC")
        ax.set_title("Signature similarity vs held-out discrimination")
        save(fig, "fig_similarity_vs_auroc", d, "High signature similarity → low discriminability.")


if __name__ == "__main__":
    raise SystemExit(main())
