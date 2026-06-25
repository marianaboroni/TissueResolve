#!/usr/bin/env python
"""Full-workflow validation after the Resolution Decision Layer refactor.

Runs the reorganized TissueResolve pipeline end-to-end on breast + HLCA/lung
donor-held-out pseudobulk and writes the QC-first report plus the required output
tables, so the evidence-based order can be verified:

    reference QC → query compatibility → trusted-resolution decision →
    modality-aware weighting → broad deconvolution → fine deconvolution only where
    supported → soft gating → QC-first report.

No model defaults are changed; donor- and seed-disjoint (reference vs query); the
query/test donors are never used for gene selection or tuning.  Outputs (all
git-ignored) under benchmarks/outputs/full_workflow_<tissue>/:
  run_metadata.json, predictions.tsv, unresolved_mass.tsv,
  trusted_resolution.tsv, qc.tsv, bulk_metrics.tsv, warnings.json,
  reference_suitability.tsv, report.html (+ figure source data).

Usage:  python full_workflow_validation.py --run-real-data --tissue both
"""
from __future__ import annotations

import argparse
import json
import sys
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

OUT = REPO / "benchmarks" / "outputs"
SEED = 3                      # held-out test seed (disjoint from any calibration)
N_SAMPLES, CELLS, MIN_CELLS = 12, 500, 20

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
               "fine_col": "cell_type", "donor_col": "donor_id", "supplement": _BREAST_SUPPLEMENT},
    "lung": {"h5ad": REPO / "examples/second_tissue_lung/data/derived/hlca_subset.h5ad",
             "hmap": REPO / "examples/second_tissue_lung/data/derived/hlca_hierarchy.tsv",
             "fine_col": "cell_type_fine", "donor_col": "donor_id"},
}


def _load(cfg):
    import anndata as ad
    adata = ad.read_h5ad(cfg["h5ad"])
    if "feature_name" in adata.var.columns:
        adata.var_names = adata.var["feature_name"].astype(str).to_numpy()
        adata.var_names_make_unique()
    hmap = pd.read_csv(cfg["hmap"], sep="\t")
    mapping = dict(zip(hmap["fine_cell_type"].astype(str), hmap["broad_cell_type"].astype(str)))
    mapping.update(cfg.get("supplement", {}))
    return adata, mapping


def run_tissue(tissue: str) -> dict:
    from tissueresolve.api import deconv_bulk, generate_report
    from tissueresolve.reference.hierarchy import build_cell_type_hierarchy
    from tissueresolve.reference.suitability import compute_reference_suitability_score
    cfg = TISSUES[tissue]
    fine_col, donor_col = cfg["fine_col"], cfg["donor_col"]
    outdir = OUT / f"full_workflow_{tissue}"
    outdir.mkdir(parents=True, exist_ok=True)
    adata, mapping = _load(cfg)

    # donor-disjoint reference / query split (query donors never used for gene selection)
    ref_d, qry = SH.split_donors(adata, donor_col, ref_frac=0.5, seed=SEED)
    rmask = adata.obs[donor_col].astype(str).isin(set(ref_d)).to_numpy()
    with warnings.catch_warnings():
        warnings.simplefilter("ignore")
        ref = H.prepare_reference(adata[rmask].copy(), cell_type_col=fine_col,
                                  min_cells=MIN_CELLS, estimate_overdispersion=False).reference
    ref_types = [str(c) for c in ref.cell_types]
    mp = dict(build_cell_type_hierarchy(ref_types, mapping))

    # 1. reference QC (suitability)
    with warnings.catch_warnings():
        warnings.simplefilter("ignore")
        suit = compute_reference_suitability_score(ref, mapping=mp)
    pd.DataFrame([{"score": suit.overall_score, "classification": suit.classification}]).to_csv(
        outdir / "reference_suitability.tsv", sep="\t", index=False)

    # held-out-donor pseudobulk (known truth)
    tgt = SH.build_target_proportions(ref_types, N_SAMPLES, "imbalanced", seed=1000 + SEED)
    ds = SH.realize_pseudobulk(adata, tgt, celltype_col=fine_col, donor_col=donor_col,
                               query_donors=qry, seed=SEED, cells_per_sample=CELLS)
    truth = ds.true_mrna_proportions

    # 2-7. hierarchical deconvolution: query compat → resolution decision →
    #      modality weighting → broad → fine-where-supported → soft gating
    with warnings.catch_warnings():
        warnings.simplefilter("ignore")
        h = deconv_bulk(ds.counts, ref, solver="nnls",
                        resolution_mode="hierarchical", hierarchy_mapping=mp)
    est, md = h.estimates, h.run_metadata

    # ----- required output tables -----
    est.combined_fine.to_csv(outdir / "predictions.tsv", sep="\t")
    (est.unresolved_mass if est.unresolved_mass.shape[1] else
     pd.DataFrame(index=est.combined_fine.index)).to_csv(outdir / "unresolved_mass.tsv", sep="\t")
    tr = md.get("trusted_resolution", {})
    rd = md.get("resolution_decision", {})
    pd.DataFrame([{"broad_family": f, "trusted_resolution": s,
                   "n_supported": len(rd.get("supported_subtypes", {}).get(f, [])),
                   "n_unsupported": len(rd.get("unsupported_subtypes", {}).get(f, []))}
                  for f, s in sorted(tr.items())]).to_csv(
        outdir / "trusted_resolution.tsv", sep="\t", index=False)
    est.qc.to_csv(outdir / "qc.tsv", sep="\t", index=False)

    # ----- diagnostic instrumentation: raw W_hat BEFORE L1-row normalisation -----
    # (opt-in solver flag; does NOT change any default output) — lets us inspect
    # sparsity / degenerate solutions (raw NNLS coefficients, row sums, effective N).
    from tissueresolve.bulk.solver import WNNLSSolver
    panel = sorted(set(map(str, ref.gene_names)) & set(map(str, ds.counts.index)))
    with warnings.catch_warnings():
        warnings.simplefilter("ignore")
        raw_res = WNNLSSolver().solve(ds.counts, ref, panel, return_raw_weights=True)
    raw_res.run_metadata["raw_weights_pre_l1norm"].round(6).to_csv(
        outdir / "raw_weights_pre_l1norm.tsv", sep="\t")
    raw_res.run_metadata["raw_weight_sparsity"].round(4).to_csv(
        outdir / "raw_weight_sparsity.tsv", sep="\t")

    # bulk accuracy vs known truth (broad + fine + conditional)
    cols = [c for c in truth.columns if c in ref_types]
    p = est.combined_fine.reindex(index=truth.index, columns=cols).fillna(0.0)
    t = truth[cols]
    acc = M.accuracy_metrics(t, p)
    tfam, pfam = M.aggregate_to_families(t, mp), M.aggregate_to_families(
        est.combined_fine.reindex(index=truth.index), mp)
    bacc = M.accuracy_metrics(tfam, pfam)
    pd.DataFrame([{"level": "fine", **acc}, {"level": "broad", **bacc}]).to_csv(
        outdir / "bulk_metrics.tsv", sep="\t", index=False)

    # run metadata + warnings
    def _jsonable(o):
        if isinstance(o, dict):
            return {str(k): _jsonable(v) for k, v in o.items()}
        if isinstance(o, (list, tuple)):
            return [_jsonable(v) for v in o]
        if isinstance(o, (np.integer,)):
            return int(o)
        if isinstance(o, (np.floating,)):
            return float(o)
        if isinstance(o, pd.DataFrame):
            return "<DataFrame>"
        return o
    (outdir / "run_metadata.json").write_text(json.dumps(_jsonable(md), indent=2, default=str))
    warns = list(getattr(h, "warnings", []) or [])
    (outdir / "warnings.json").write_text(json.dumps(warns, indent=2, default=str))

    # 8. QC-first report (figure source data written alongside by the report builder)
    report_ok, report_err = True, ""
    try:
        with warnings.catch_warnings():
            warnings.simplefilter("ignore")
            generate_report(h, outdir / "report.html")
    except Exception as exc:           # report failures must be visible, not hidden
        report_ok, report_err = False, str(exc)

    summary = {
        "tissue": tissue, "n_ref_donors": len(ref_d), "n_query_donors": len(qry),
        "n_fine_types": len(ref_types), "n_families": int(len(set(mp.values()))),
        "reference_suitability": f"{suit.overall_score:.3f} ({suit.classification})",
        "modality": md.get("modality"), "prediction_unit": md.get("prediction_unit"),
        "gene_weighting_mode": md.get("gene_weighting_mode"),
        "hierarchical_gating": md.get("hierarchical_gating"),
        "n_full_fine": rd.get("n_full_fine"), "n_selected_fine": rd.get("n_selected_fine"),
        "n_broad_only": rd.get("n_broad_only"),
        "fine_trusted": md.get("fine_predictions_trusted"),
        "fine_diagnostic_only": md.get("fine_predictions_diagnostic_only"),
        "mass_error": est.metadata.get("mass_conservation_max_error"),
        "unresolved_mass_fraction": est.metadata.get("unresolved_mass_fraction"),
        "fine_pearson": round(acc["pearson"], 4), "fine_rmse": round(acc["rmse"], 4),
        "broad_pearson": round(bacc["pearson"], 4), "broad_rmse": round(bacc["rmse"], 4),
        "report_generated": report_ok, "report_error": report_err,
    }
    (outdir / "summary.json").write_text(json.dumps(summary, indent=2, default=str))
    print(f"\n=== [{tissue}] full-workflow summary ===")
    for k, v in summary.items():
        print(f"  {k}: {v}")
    return summary


def main(argv=None) -> int:
    ap = argparse.ArgumentParser()
    ap.add_argument("--run-real-data", action="store_true")
    ap.add_argument("--tissue", choices=["breast", "lung", "both"], default="both")
    args = ap.parse_args(argv)
    if not args.run_real_data:
        print("Refusing without --run-real-data.", file=sys.stderr); return 2
    tissues = ["breast", "lung"] if args.tissue == "both" else [args.tissue]
    summaries = []
    for tissue in tissues:
        if not Path(TISSUES[tissue]["h5ad"]).exists():
            print(f"Missing {TISSUES[tissue]['h5ad']}; skipping {tissue}.", file=sys.stderr); continue
        summaries.append(run_tissue(tissue))
    OUT.mkdir(parents=True, exist_ok=True)
    pd.DataFrame(summaries).to_csv(OUT / "full_workflow_summary.tsv", sep="\t", index=False)
    print(f"\nwrote {OUT/'full_workflow_summary.tsv'}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
