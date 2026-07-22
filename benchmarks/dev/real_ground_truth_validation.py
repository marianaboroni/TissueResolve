#!/usr/bin/env python
"""Real ground-truth validation harness (bulk↔flow, spatial↔Xenium/CosMx segmentation).

Consumes LOCAL data only (never downloads); refuses without --run-real-data; if the
required files are absent it prints the exact expected paths + an acquisition checklist
and writes a manifest template, then exits non-zero. When data is present it deconvolves
with TissueResolve, aggregates estimates to the ground-truth label level via a
label_map, scores against the real truth (cross-subject / cross-spot), and applies gates.
No results are ever fabricated. Outputs are gitignored.

See docs/dev/REAL_GROUND_TRUTH_VALIDATION_PLAN.md.

Usage:
  # bulk vs flow/CyTOF
  PYTHONPATH=src:. python benchmarks/dev/real_ground_truth_validation.py \
      --run-real-data --mode bulk_flow --data-dir examples/real_bulk_flow/data
  # spatial vs Xenium/CosMx segmentation-derived per-spot truth
  PYTHONPATH=src:. python benchmarks/dev/real_ground_truth_validation.py \
      --run-real-data --mode spatial_seg --data-dir examples/spatial_ground_truth/data
"""
from __future__ import annotations

import argparse
import json
import sys
from pathlib import Path

import numpy as np
import pandas as pd

REPO = Path(__file__).resolve().parents[2]
sys.path.insert(0, str(REPO / "examples" / "real_breast_cancer" / "scripts"))
sys.path.insert(0, str(REPO))
sys.path.insert(0, str(REPO / "src"))
import _harness as H  # noqa: E402
from benchmarks.shared import metrics as M  # noqa: E402

OUT = REPO / "benchmarks" / "outputs" / "real_ground_truth"

REQUIRED = {
    "bulk_flow": ["bulk_counts.tsv", "flow_truth.tsv", "reference.h5ad", "label_map.tsv"],
    "spatial_seg": ["spot_truth.tsv", "reference.h5ad", "label_map.tsv"],  # +visium.h5ad or counts+coords
}
GATES = {"min_pearson": 0.6, "min_ccc": 0.5, "max_abs_bias": 0.05,
         "spatial_min_pearson": 0.5}


def _acquisition_help(mode, data_dir):
    files = REQUIRED[mode] + (["visium.h5ad (or counts.tsv + coords.tsv)"] if mode == "spatial_seg" else [])
    print(f"\n[real ground-truth validation | mode={mode}] required data NOT found.\n"
          f"Expected under: {data_dir}/", file=sys.stderr)
    for f in files:
        print(f"    - {f}", file=sys.stderr)
    print("\nAcquisition checklist (see docs/dev/REAL_GROUND_TRUTH_VALIDATION_PLAN.md):",
          file=sys.stderr)
    if mode == "bulk_flow":
        print("  * Obtain a cohort with matched bulk RNA-seq + flow/CyTOF proportions\n"
              "    (e.g. Finotello 2019 / EPIC-Racle 2017 / CIBERSORT-Newman / ImmPort\n"
              "    SDY311/SDY420 — verify accessions).\n"
              "  * bulk_counts.tsv: genes x samples; flow_truth.tsv: samples x flow_labels\n"
              "    (fractions); label_map.tsv: flow_label <TAB> reference_cell_type(s '|'-sep).",
              file=sys.stderr)
    else:
        print("  * Obtain matched Visium + Xenium/CosMx (e.g. 10x breast Janesick 2023);\n"
              "    upstream, co-register and turn segmentation calls into per-spot\n"
              "    proportions -> spot_truth.tsv (spot x cell_type fractions).", file=sys.stderr)
    # write a manifest template (does not overwrite an existing one)
    data_dir.mkdir(parents=True, exist_ok=True)
    mp = data_dir / "dataset_manifest.json"
    if not mp.exists():
        mp.write_text(json.dumps({
            "name": "", "source": "", "accession": "", "url": "", "access_date": "",
            "platform": "", "n_samples_or_spots": None, "filters_applied": "",
            "package_versions": "", "notes": "fill in when you place the data"}, indent=2),
            encoding="utf-8")
        print(f"\nWrote manifest template -> {mp} (fill it in).", file=sys.stderr)


def _load_label_map(path):
    m = {}
    for line in Path(path).read_text(encoding="utf-8").splitlines():
        if not line.strip() or line.startswith("#"):
            continue
        parts = line.split("\t")
        if len(parts) < 2:
            continue
        m[parts[0].strip()] = [x.strip() for x in parts[1].split("|") if x.strip()]
    return m


def _aggregate_to_truth_labels(pred_df, label_map):
    """Sum predicted reference cell types into ground-truth labels via label_map."""
    out = {}
    for tlabel, refcts in label_map.items():
        cols = [c for c in refcts if c in pred_df.columns]
        out[tlabel] = pred_df[cols].sum(axis=1) if cols else 0.0
    agg = pd.DataFrame(out, index=pred_df.index)
    return agg


def _score_truth(truth, pred):
    common_cols = [c for c in truth.columns if c in pred.columns]
    common_idx = [i for i in truth.index if i in pred.index]
    t = truth.loc[common_idx, common_cols]
    p = pred.reindex(index=common_idx, columns=common_cols).fillna(0.0)
    rows = []
    for c in common_cols:
        tv, pv = t[c].to_numpy(float), p[c].to_numpy(float)
        rows.append({"label": c,
                     "pearson": M._safe_corr(tv, pv, "pearson"),
                     "ccc": M.concordance_correlation_coefficient(tv, pv),
                     "rmse": float(np.sqrt(np.mean((tv - pv) ** 2))),
                     "bias": float(np.mean(pv - tv)), "n": len(tv)})
    per = pd.DataFrame(rows)
    overall = {"mean_pearson": float(per["pearson"].mean()),
               "mean_ccc": float(per["ccc"].mean()),
               "rmse": float(np.sqrt(np.mean((t.to_numpy(float) - p.to_numpy(float)) ** 2))),
               "max_abs_bias": float(per["bias"].abs().max()),
               "n_labels": len(common_cols), "n_units": len(common_idx)}
    return per, overall


def run_bulk_flow(data_dir):
    import anndata as ad
    import tissueresolve as tr
    from tissueresolve.config import TissueResolveConfig
    bulk = pd.read_csv(data_dir / "bulk_counts.tsv", sep="\t", index_col=0)
    flow = pd.read_csv(data_dir / "flow_truth.tsv", sep="\t", index_col=0)
    label_map = _load_label_map(data_dir / "label_map.tsv")
    adata = ad.read_h5ad(data_dir / "reference.h5ad")
    with __import__("warnings").catch_warnings():
        __import__("warnings").simplefilter("ignore")
        ref = H.prepare_reference(adata, min_cells=30).reference
    results = {}
    for method in ("wNNLS", "poisson_glm_experimental"):
        cfg = TissueResolveConfig(); cfg.bulk_solver.method = method
        pred = tr.deconv_bulk(bulk, ref, config=cfg, resolution_mode="none").deconv.proportions
        agg = _aggregate_to_truth_labels(pred, label_map)
        per, overall = _score_truth(flow, agg)
        results[method] = {"per_label": per, "overall": overall}
    return results


def run_spatial_seg(data_dir):
    import anndata as ad
    import tissueresolve as tr
    from tissueresolve.config import TissueResolveConfig
    truth = pd.read_csv(data_dir / "spot_truth.tsv", sep="\t", index_col=0)
    label_map = _load_label_map(data_dir / "label_map.tsv")
    ref = H.prepare_reference(ad.read_h5ad(data_dir / "reference.h5ad"), min_cells=30).reference
    vis = ad.read_h5ad(data_dir / "visium.h5ad")
    Y = np.asarray(vis.X.todense()) if hasattr(vis.X, "todense") else np.asarray(vis.X)
    ar = vis.obs["array_row"].to_numpy() if "array_row" in vis.obs else np.arange(vis.n_obs)
    ac = vis.obs["array_col"].to_numpy() if "array_col" in vis.obs else np.zeros(vis.n_obs)
    lib = Y.sum(1)
    res = tr.deconv_spatial(Y, ref, ar, ac, lib, list(vis.var_names),
                            spot_ids=list(vis.obs_names), resolution_mode="none")
    agg = _aggregate_to_truth_labels(res.deconv.proportions, label_map)
    per, overall = _score_truth(truth, agg)
    return {"TissueResolve_spatial": {"per_label": per, "overall": overall}}


def main(argv=None) -> int:
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument("--run-real-data", action="store_true")
    ap.add_argument("--mode", required=True, choices=["bulk_flow", "spatial_seg"])
    ap.add_argument("--data-dir", required=True)
    args = ap.parse_args(argv)
    if not H.real_data_enabled(args.run_real_data):
        print(f"Refusing without --run-real-data (or {H.REAL_DATA_ENV}=1).", file=sys.stderr)
        return 2
    data_dir = Path(args.data_dir)
    needed = list(REQUIRED[args.mode])
    if args.mode == "spatial_seg" and not (data_dir / "visium.h5ad").exists():
        needed = needed + ["visium.h5ad"]
    missing = [f for f in REQUIRED[args.mode] if not (data_dir / f).exists()]
    if args.mode == "spatial_seg" and not (data_dir / "visium.h5ad").exists():
        missing.append("visium.h5ad")
    if missing:
        _acquisition_help(args.mode, data_dir)
        return 3

    OUT.mkdir(parents=True, exist_ok=True)
    results = run_bulk_flow(data_dir) if args.mode == "bulk_flow" else run_spatial_seg(data_dir)

    # write outputs + gates
    all_per, overall_rows, gate_rows = [], [], []
    for method, r in results.items():
        per = r["per_label"].copy(); per.insert(0, "method", method)
        all_per.append(per)
        overall_rows.append({"method": method, **r["overall"]})
        o = r["overall"]
        if args.mode == "bulk_flow":
            passed = (o["mean_pearson"] >= GATES["min_pearson"] and o["mean_ccc"] >= GATES["min_ccc"]
                      and o["max_abs_bias"] <= GATES["max_abs_bias"])
        else:
            passed = o["mean_pearson"] >= GATES["spatial_min_pearson"]
        gate_rows.append({"method": method, "mode": args.mode, "mean_pearson": round(o["mean_pearson"], 3),
                          "mean_ccc": round(o["mean_ccc"], 3), "max_abs_bias": round(o["max_abs_bias"], 3),
                          "gate_passed": bool(passed)})
    pd.concat(all_per, ignore_index=True).to_csv(OUT / f"{args.mode}_per_label_metrics.tsv", sep="\t", index=False)
    pd.DataFrame(overall_rows).to_csv(OUT / f"{args.mode}_overall_metrics.tsv", sep="\t", index=False)
    pd.DataFrame(gate_rows).to_csv(OUT / f"{args.mode}_gates.tsv", sep="\t", index=False)
    manifest = data_dir / "dataset_manifest.json"
    (OUT / f"{args.mode}_run_manifest.json").write_text(json.dumps({
        "mode": args.mode, "data_dir": str(data_dir),
        "dataset_manifest": json.loads(manifest.read_text()) if manifest.exists() else None,
        "gates": GATES}, indent=2), encoding="utf-8")
    print(f"\nWrote -> {OUT}/")
    print(pd.DataFrame(gate_rows).to_string(index=False))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
