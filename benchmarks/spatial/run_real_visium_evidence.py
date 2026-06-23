#!/usr/bin/env python
"""Real Visium breast-cancer evidence — NO GROUND TRUTH (PART 13).

Runs deconvolution on the real 10x Visium breast-cancer section and reports
*evidence* metrics only — NEVER accuracy.  Per the brief: concordance is not
accuracy; marker recovery is a proxy; high Moran's I / smoothness is not
automatically better.

Metrics per method: marker recovery (proxy), Moran's I + Geary's C (structure),
proportion entropy, dominant fraction, near-zero fraction, reconstruction Pearson
(observed vs reference-reconstructed CPM, a self-consistency check), runtime,
output completeness.  Cross-method: pairwise fine + family-level concordance.

Also exports sc reference + spot counts so RCTD / cell2location can be run on the
SAME spots (cross-method concordance with established tools).

Usage:  python benchmarks/spatial/run_real_visium_evidence.py --run-real-data
"""
from __future__ import annotations

import argparse
import json
import sys
import time
import warnings
from pathlib import Path

import numpy as np
import pandas as pd

REPO = Path(__file__).resolve().parents[2]
RBC = REPO / "examples" / "real_breast_cancer"
sys.path.insert(0, str(RBC / "scripts"))
sys.path.insert(0, str(REPO))
import _harness as H  # noqa: E402
from benchmarks.shared import metrics as M  # noqa: E402
from benchmarks.shared import spatial_metrics as SM  # noqa: E402
from benchmarks.shared import spatial_multimetric_ranking as RANK  # noqa: E402

OUT = REPO / "benchmarks" / "outputs" / "spatial_real_visium"
RAW = OUT / "raw"
EI = OUT / "external_inputs"
MAX_SPOTS = 500
SEED = 0
N_MARKERS = 15


def _load_scenario():
    import anndata as ad
    import scipy.sparse as sp
    from tissueresolve.results import ReferenceSignature
    from tissueresolve.reference.hierarchy import load_hierarchy_mapping, build_cell_type_hierarchy
    ref = ReferenceSignature.load(RBC / "outputs" / "reference" / "breast_cancer_reference")
    adata = ad.read_h5ad(RBC / "data" / "spatial" / "human_breast_cancer_1.h5ad")
    n_total = adata.n_obs
    rng = np.random.default_rng(SEED)
    if n_total > MAX_SPOTS:
        keep = np.sort(rng.choice(n_total, size=MAX_SPOTS, replace=False))
        adata = adata[keep].copy()
    Y = adata.X
    Y = np.asarray(Y.todense()) if sp.issparse(Y) else np.asarray(Y)
    gene_names = [str(g) for g in adata.var_names]
    lib = Y.sum(axis=1).astype("float32")
    with warnings.catch_warnings():
        warnings.simplefilter("ignore")
        mapping = build_cell_type_hierarchy(list(ref.cell_types),
                                            load_hierarchy_mapping(RBC / "config" /
                                            "breast_cancer_cell_type_hierarchy.tsv"))
    coords = np.c_[adata.obs["array_row"].to_numpy(), adata.obs["array_col"].to_numpy()].astype(float)
    return {"ref": ref, "Y": Y, "gene_names": gene_names, "lib": lib,
            "array_row": adata.obs["array_row"].to_numpy(),
            "array_col": adata.obs["array_col"].to_numpy(),
            "spot_ids": list(map(str, adata.obs_names)), "coords": coords,
            "mapping": mapping, "n_total": int(n_total)}


def _marker_sets(ref, visium_genes, n=N_MARKERS):
    R = ref.as_R_cpm(); types = [str(c) for c in ref.cell_types]
    genes = [str(g) for g in ref.gene_names]
    vis = set(visium_genes)
    mean_across = R.mean(axis=0) + 1e-9
    out = {}
    for i, t in enumerate(types):
        spec = R[i] / mean_across
        order = np.argsort(spec)[::-1]
        picks = [genes[j] for j in order if genes[j] in vis][:n]
        if picks:
            out[t] = picks
    return out


def _reconstruction_pearson(Y, props, ref, gene_names):
    from scipy.stats import pearsonr
    rg = {g: i for i, g in enumerate(map(str, ref.gene_names))}
    shared = [g for g in gene_names if g in rg]
    if len(shared) < 10:
        return float("nan")
    rows = [rg[g] for g in shared]
    cidx = [gene_names.index(g) for g in shared]
    R = np.asarray(ref.as_R_cpm())
    P = props.reindex(columns=list(ref.cell_types)).fillna(0.0).to_numpy(float)
    obs = Y[:, cidx]
    obs_cpm = obs / np.clip(obs.sum(1, keepdims=True), 1, None) * 1e6
    recon = P @ R[:, rows]
    vals = []
    for s in range(obs_cpm.shape[0]):
        o, r = np.log1p(obs_cpm[s]), np.log1p(recon[s])
        if o.std() > 1e-9 and r.std() > 1e-9:
            vals.append(pearsonr(o, r)[0])
    return float(np.mean(vals)) if vals else float("nan")


def _evidence(props, sc, marker_expr, marker_sets):
    coords = sc["coords"]
    nbr, w = SM.knn_weights(coords, 6)
    mi = props.apply(lambda c: SM.morans_i(c.to_numpy(float), nbr, w))
    gc = props.apply(lambda c: SM.gearys_c(c.to_numpy(float), nbr, w))
    mr = RANK.marker_recovery_score(props, marker_expr, marker_sets)
    ent = M.proportion_entropy(props)
    return {
        "marker_recovery_mean": float(np.nanmean(mr.to_numpy())),
        "mean_morans_i": float(np.nanmean(mi.to_numpy())),
        "mean_gearys_c": float(np.nanmean(gc.to_numpy())),
        "mean_entropy": float(ent.mean()),
        "mean_dominant_fraction": float(props.max(axis=1).mean()),
        "near_zero_fraction": M.near_zero_fraction(props),
        "recon_pearson": _reconstruction_pearson(sc["Y"], props, sc["ref"], sc["gene_names"]),
        "n_types_detected_mean": float((props > 0.01).sum(axis=1).mean()),
        "output_complete": bool(props.notna().all().all()),
    }


def main(argv=None) -> int:
    ap = argparse.ArgumentParser()
    ap.add_argument("--run-real-data", action="store_true")
    args = ap.parse_args(argv)
    if not H.real_data_enabled(args.run_real_data):
        print(f"Refusing without --run-real-data (or {H.REAL_DATA_ENV}=1).", file=sys.stderr)
        return 2

    import tissueresolve as tr
    from tissueresolve.api import deconv_bulk
    RAW.mkdir(parents=True, exist_ok=True); EI.mkdir(parents=True, exist_ok=True)
    sc = _load_scenario()
    ref, mapping = sc["ref"], sc["mapping"]
    print(f"Real Visium: {len(sc['spot_ids'])}/{sc['n_total']} spots, "
          f"{len(sc['gene_names'])} genes, {len(list(ref.cell_types))} ref types. NO GROUND TRUTH.")

    # marker expression matrix (spots × marker genes, CPM log1p)
    marker_sets = _marker_sets(ref, sc["gene_names"])
    all_markers = sorted({g for v in marker_sets.values() for g in v})
    gidx = {g: i for i, g in enumerate(sc["gene_names"])}
    Ym = sc["Y"][:, [gidx[g] for g in all_markers]]
    Ym_cpm = np.log1p(Ym / np.clip(Ym.sum(1, keepdims=True), 1, None) * 1e6)
    marker_expr = pd.DataFrame(Ym_cpm, index=sc["spot_ids"], columns=all_markers)

    # export inputs for external methods (same spots)
    pd.DataFrame(sc["Y"].T.astype(int), index=sc["gene_names"], columns=sc["spot_ids"]).to_csv(
        EI / "spot_counts_genes_by_spots.tsv", sep="\t")
    pd.DataFrame({"spot": sc["spot_ids"], "row": sc["array_row"], "col": sc["array_col"]}).to_csv(
        EI / "spot_coords.tsv", sep="\t", index=False)

    methods, preds, rows = {}, {}, []

    def run(name, fn):
        t0 = time.perf_counter()
        try:
            with warnings.catch_warnings():
                warnings.simplefilter("ignore")
                p = fn()
        except Exception as exc:  # noqa: BLE001
            print(f"  {name:34s} FAILED: {exc}"); return
        rt = time.perf_counter() - t0
        p = p.reindex(index=sc["spot_ids"])
        preds[name] = p
        p.to_csv(RAW / f"{name}_proportions.tsv", sep="\t")
        ev = _evidence(p, sc, marker_expr, marker_sets)
        rows.append({"method": name, "runtime_seconds": round(rt, 1), "status": "executed", **ev})
        print(f"  {name:34s} {rt:5.1f}s marker_rec={ev['marker_recovery_mean']:.3f} "
              f"moran={ev['mean_morans_i']:.3f} recon_r={ev['recon_pearson']:.3f}")

    bulk_df = pd.DataFrame(sc["Y"].T, index=sc["gene_names"], columns=sc["spot_ids"])
    run("NNLS_per_spot", lambda: deconv_bulk(bulk_df, ref, solver="nnls",
                                             resolution_mode="none").deconv.proportions)
    run("TissueResolve_spatial_flat", lambda: tr.deconv_spatial(
        sc["Y"], ref, sc["array_row"], sc["array_col"], sc["lib"], sc["gene_names"],
        spot_ids=sc["spot_ids"], resolution_mode="none", run_neighbourhood=False).deconv.proportions)
    run("TissueResolve_spatial_hierarchical", lambda: tr.deconv_spatial(
        sc["Y"], ref, sc["array_row"], sc["array_col"], sc["lib"], sc["gene_names"],
        spot_ids=sc["spot_ids"], resolution_mode="hierarchical", hierarchy_mapping=mapping,
        run_neighbourhood=False).estimates.combined_fine)

    # cross-method concordance (fine + family) — NOT accuracy
    fine_preds = {k: v.reindex(columns=list(ref.cell_types)).fillna(0.0) for k, v in preds.items()}
    conc = M.pairwise_method_correlation(fine_preds)
    fam_preds = {k: M.aggregate_to_families(v, mapping) for k, v in preds.items()}
    conc_fam = M.pairwise_method_correlation(fam_preds)

    OUT.mkdir(parents=True, exist_ok=True)
    ev_df = pd.DataFrame(rows)
    with open(OUT / "spatial_real_metrics.tsv", "w") as fh:
        fh.write("# Real Visium breast cancer — NO GROUND TRUTH. Evidence metrics only "
                 "(marker recovery=proxy; concordance!=accuracy; high Moran's I not automatically better).\n")
        ev_df.to_csv(fh, sep="\t", index=False)
    pd.DataFrame([{"pair": k, "fine_pearson": v} for k, v in conc.items()]).to_csv(
        OUT / "method_concordance_fine.tsv", sep="\t", index=False)
    pd.DataFrame([{"pair": k, "family_pearson": v} for k, v in conc_fam.items()]).to_csv(
        OUT / "method_concordance_family.tsv", sep="\t", index=False)
    (OUT / "manifest.json").write_text(json.dumps({
        "dataset": "10x Visium V1_Breast_Cancer_Block_A_Section_1 (human_breast_cancer_1.h5ad)",
        "n_spots": len(sc["spot_ids"]), "n_spots_total": sc["n_total"], "max_spots": MAX_SPOTS,
        "seed": SEED, "n_ref_types": len(list(ref.cell_types)), "n_markers_per_type": N_MARKERS,
        "ground_truth": "NONE — evidence metrics only, not accuracy",
    }, indent=2), encoding="utf-8")

    _figures(ev_df, conc, preds, sc)
    print("\n=== evidence (NO GROUND TRUTH) ===")
    print(ev_df[["method", "marker_recovery_mean", "mean_morans_i", "mean_entropy",
                 "near_zero_fraction", "recon_pearson"]].round(3).to_string(index=False))
    print("\n=== fine concordance (method agreement, NOT accuracy) ===")
    for k, v in conc.items():
        print(f"  {k}: {v:.3f}")
    print(f"\nWrote -> {OUT}/")
    return 0


def _figures(ev_df, conc, preds, sc):
    import matplotlib; matplotlib.use("Agg")
    import matplotlib.pyplot as plt
    FIG = OUT / "figures"; FIG.mkdir(parents=True, exist_ok=True)

    def save(fig, name, data, cap):
        fig.tight_layout(); fig.savefig(FIG / f"{name}.png", dpi=150)
        fig.savefig(FIG / f"{name}.svg"); plt.close(fig)
        data.to_csv(FIG / f"{name}.data.tsv", sep="\t")
        (FIG / f"{name}.caption.txt").write_text(cap.strip() + "\n", encoding="utf-8")

    sub = ev_df.set_index("method")[["marker_recovery_mean", "mean_morans_i", "recon_pearson"]]
    fig, ax = plt.subplots(figsize=(8, 4.2)); sub.plot.bar(ax=ax)
    ax.set_title("Real Visium evidence (NOT accuracy): marker recovery / Moran's I / reconstruction")
    ax.set_ylabel("score"); plt.setp(ax.get_xticklabels(), rotation=15, ha="right")
    save(fig, "figL_real_evidence", sub,
         "Real Visium breast cancer evidence metrics (no ground truth): marker-recovery "
         "proxy, mean Moran's I (structure), reconstruction Pearson (self-consistency). "
         "NOT accuracy.")


if __name__ == "__main__":
    raise SystemExit(main())
