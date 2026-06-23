#!/usr/bin/env python
"""STAGE 1b — DECONVOLUTION identifiability ceiling via 2-signature mixture recovery.

The cell-classification AUROC (identifiability_ceiling.py) is ~1.0 for every
family — but that is the wrong ceiling for *deconvolution*.  Deconvolution must
recover subtype PROPORTIONS from a mixed signal, which is governed by the
collinearity of the subtype *signatures*, not cell separability.  Here, for each
within-family pair we build the two reference signatures, synthesise mixtures at
known ratios with NB-like noise, recover the ratio by NNLS, and measure the
recovery error — the true deconvolution ceiling.  Also reports the shared-lineage
fraction ‖B‖/‖S‖ (Stage-3 motivation: does the shared program dominate the contrast?).

Outputs: benchmarks/outputs/mixture_recovery_ceiling_pairs.tsv,
mixture_recovery_family_summary.tsv; figure.
Usage:  python benchmarks/diagnostics/mixture_recovery_ceiling.py --run-real-data
"""
from __future__ import annotations

import argparse
import sys
from itertools import combinations
from pathlib import Path

import numpy as np
import pandas as pd
from scipy.optimize import nnls

REPO = Path(__file__).resolve().parents[2]
sys.path.insert(0, str(REPO / "examples" / "real_breast_cancer" / "scripts"))
sys.path.insert(0, str(REPO))
import _harness as H  # noqa: E402

OUT = REPO / "benchmarks" / "outputs"
FIG = OUT / "figures" / "identifiability"
HIER = REPO / "examples" / "real_breast_cancer" / "config" / "breast_cancer_cell_type_hierarchy.tsv"
TCGA = H.DATA_DIR / "bulk_tcga_tnbc" / "tcga_tnbc_counts.tsv"
RATIOS = [0.1, 0.25, 0.5, 0.75, 0.9]
N_NOISE = 20
DEPTH = 5e5
SEED = 0
THRESH = {"identifiable": 0.10, "conditional": 0.20}   # ratio-recovery MAE thresholds


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
    mapped = np.isin(ct, list(raw_map))
    mapping = build_cell_type_hierarchy(sorted(set(ct[mapped])), raw_map)
    X = adata.X
    X = np.asarray(X.todense()) if hasattr(X, "todense") else np.asarray(X, float)
    # CPM signatures per subtype (linear, the deconvolution space)
    qset = set(genes)
    if TCGA.exists():
        bulk = pd.read_csv(TCGA, sep="\t", index_col=0, comment="#")
        det = (bulk > 0).mean(axis=1); qset = {g for g in genes if g in set(map(str, det.index[det >= 0.5]))}
    qidx = np.array([i for i, g in enumerate(genes) if g in qset])

    def sig(st):
        ii = np.where((ct == st) & mapped)[0]
        cpm = X[ii] / np.clip(X[ii].sum(1, keepdims=True), 1, None) * 1e6
        return cpm.mean(0)[qidx]

    fam_members = {}
    for c in sorted(set(ct[mapped])):
        fam_members.setdefault(mapping.get(c, c), []).append(c)
    fam_members = {f: m for f, m in fam_members.items() if len(m) >= 2}

    rng = np.random.default_rng(SEED)
    rows = []
    for fam, members in fam_members.items():
        sigs = {m: sig(m) for m in members}
        for a, b in combinations(members, 2):
            Sa, Sb = sigs[a], sigs[b]
            Smat = np.vstack([Sa, Sb]).T               # G×2
            # shared-lineage fraction: B = mean(Sa,Sb); D = S - B
            Bshared = 0.5 * (Sa + Sb)
            shared_frac = float(np.linalg.norm(Bshared) /
                                (0.5 * (np.linalg.norm(Sa) + np.linalg.norm(Sb)) + 1e-12))
            errs = []
            for r in RATIOS:
                theta = np.array([r, 1 - r])
                mu = Smat @ theta
                for _ in range(N_NOISE):
                    y = rng.poisson(np.clip(mu / mu.sum() * DEPTH, 0, None)).astype(float)
                    th, _ = nnls(Smat, y)
                    s = th.sum(); th = th / s if s > 0 else theta
                    errs.append(abs(th[0] - r))
            mae = float(np.mean(errs))
            rows.append({"family": fam, "subtype_a": a, "subtype_b": b,
                         "ratio_recovery_mae": round(mae, 4),
                         "signature_pearson": round(float(np.corrcoef(Sa, Sb)[0, 1]), 4),
                         "condition_number": round(float(np.linalg.cond(Smat)), 2),
                         "shared_lineage_fraction": round(shared_frac, 4)})
        print(f"  {fam}: done", flush=True)
    pairs = pd.DataFrame(rows)
    OUT.mkdir(parents=True, exist_ok=True)
    pairs.to_csv(OUT / "mixture_recovery_ceiling_pairs.tsv", sep="\t", index=False)

    def _bucket(mae):
        return ("fine" if mae <= THRESH["identifiable"]
                else "selected_fine" if mae <= THRESH["conditional"] else "broad_only")
    fam_rows = []
    for fam, g in pairs.groupby("family"):
        med = float(g["ratio_recovery_mae"].median())
        fam_rows.append({"family": fam, "n_pairs": len(g),
                         "median_ratio_recovery_mae": round(med, 4),
                         "median_signature_pearson": round(float(g["signature_pearson"].median()), 4),
                         "median_condition_number": round(float(g["condition_number"].median()), 2),
                         "median_shared_lineage_fraction": round(float(g["shared_lineage_fraction"].median()), 4),
                         "deconv_recommended_resolution": _bucket(med),
                         "frac_pairs_broad_only": round(float((g["ratio_recovery_mae"] > THRESH["conditional"]).mean()), 3)})
    famsum = pd.DataFrame(fam_rows).sort_values("median_ratio_recovery_mae", ascending=False)
    famsum.to_csv(OUT / "mixture_recovery_family_summary.tsv", sep="\t", index=False)

    # figure: ratio-recovery MAE vs signature correlation
    import matplotlib; matplotlib.use("Agg")
    import matplotlib.pyplot as plt
    FIG.mkdir(parents=True, exist_ok=True)
    fig, ax = plt.subplots(figsize=(6, 4.5))
    sc = ax.scatter(pairs["signature_pearson"], pairs["ratio_recovery_mae"],
                    c=pairs["shared_lineage_fraction"], cmap="viridis", alpha=0.8)
    ax.axhline(0.10, color="green", ls="--"); ax.axhline(0.20, color="orange", ls="--")
    ax.set_xlabel("subtype signature Pearson (collinearity)")
    ax.set_ylabel("2-signature ratio-recovery MAE (deconvolution ceiling)")
    ax.set_title("Deconvolution identifiability ceiling (mixture recovery)")
    fig.colorbar(sc, ax=ax, label="shared-lineage fraction")
    fig.tight_layout(); fig.savefig(FIG / "fig_mixture_recovery_ceiling.png", dpi=150)
    fig.savefig(FIG / "fig_mixture_recovery_ceiling.svg"); plt.close(fig)
    pairs.to_csv(FIG / "fig_mixture_recovery_ceiling.data.tsv", sep="\t", index=False)
    (FIG / "fig_mixture_recovery_ceiling.caption.txt").write_text(
        "2-signature NNLS ratio-recovery error vs subtype signature collinearity. "
        "High collinearity (shared-lineage dominance) → high deconvolution error even "
        "when cells are perfectly classifiable. Green/orange = 0.10/0.20 MAE thresholds.\n")

    print("\n=== DECONVOLUTION ceiling by family (mixture recovery) ===")
    print(famsum.to_string(index=False))
    frac_broad = float((pairs["ratio_recovery_mae"] > THRESH["conditional"]).mean())
    print(f"\nPairs requiring broad-only (recovery MAE>0.20): {frac_broad:.2f}")
    print(f"Median shared-lineage fraction: {pairs['shared_lineage_fraction'].median():.3f} "
          f"(high → contrast is a small residual → Stage-3 motivation)")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
