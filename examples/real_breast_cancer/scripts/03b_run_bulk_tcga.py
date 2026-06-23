#!/usr/bin/env python
"""
Run TissueResolve bulk deconvolution on **real TCGA-BRCA TNBC** bulk RNA-seq.

Real bulk tumours have **no ground-truth cell proportions**, so this is a
*no-ground-truth* benchmark: it reports cross-method **concordance**, per-sample
**reconstruction QC**, **runtime**, gene overlap, and the mean predicted
composition / dominant cell types — never accuracy.  (For ground-truth accuracy
use the pseudobulk harness, scripts 02–03.)

Inputs
------
- ``data/bulk_tcga_tnbc/tcga_tnbc_counts.tsv``  (from 00b_download_tcga_tnbc_bulk.py)
- the saved single-cell reference signature (outputs/reference/breast_cancer_reference)

Outputs (under ``outputs/bulk_tcga/``)
--------------------------------------
- ``tcga_bulk_estimated_proportions.tsv``  (TissueResolve auto solver)
- ``tcga_bulk_method_benchmark.tsv``       (per-method runtime, mean recon R²,
                                            concordance vs the auto solver)
- ``tcga_bulk_mean_composition.tsv``       (mean RNA-derived proportion / type)
- ``tcga_bulk_dominant_type_counts.tsv``
- ``warnings.json``, ``methods.txt``, ``report.html``

Usage
-----
    python scripts/03b_run_bulk_tcga.py --run-real-data
"""
from __future__ import annotations

import argparse
import sys
import time
import warnings
from pathlib import Path

import numpy as np
import pandas as pd

sys.path.insert(0, str(Path(__file__).resolve().parent))
import _harness as H  # noqa: E402

COUNTS_TSV = H.DATA_DIR / "bulk_tcga_tnbc" / "tcga_tnbc_counts.tsv"
OUT_DIR = H.OUTPUTS_DIR / "bulk_tcga"

SOLVERS = ["auto", "nnls", "weighted_nnls"]


def _recon_pearson(bulk_df: pd.DataFrame, props: pd.DataFrame, ref) -> float:
    """Mean per-sample reconstruction concordance (no ground truth needed).

    Reconstructs each sample's expression as ``sum_k proportion_k · reference_k``
    (reference profiles are CPM) and correlates it with the observed profile
    (CPM, log1p) over genes shared with the reference.  Higher = the estimated
    composition explains the observed bulk better (a self-consistency check).
    """
    from scipy.stats import pearsonr
    ref_genes = list(map(str, ref.gene_names))
    gidx = {g: i for i, g in enumerate(ref_genes)}
    shared = [g for g in map(str, bulk_df.index) if g in gidx]
    if len(shared) < 10:
        return float("nan")
    rows = [gidx[g] for g in shared]
    R = np.asarray(ref.R_cpm)               # (K cell types × G genes), CPM
    cts = list(ref.cell_types)
    P = props.reindex(columns=cts).fillna(0.0).to_numpy(float)  # samples × K
    obs = bulk_df.loc[shared].to_numpy(float).T                 # samples × shared
    obs_cpm = obs / np.clip(obs.sum(1, keepdims=True), 1, None) * 1e6
    recon = P @ R[:, rows]                                       # samples × shared (CPM)
    vals = []
    for s in range(obs_cpm.shape[0]):
        o, r = np.log1p(obs_cpm[s]), np.log1p(recon[s])
        if np.std(o) > 1e-9 and np.std(r) > 1e-9:
            vals.append(pearsonr(o, r)[0])
    return float(np.mean(vals)) if vals else float("nan")


def _flatten_concordance(a: pd.DataFrame, b: pd.DataFrame) -> dict:
    from scipy.stats import pearsonr, spearmanr
    samples = [s for s in a.index if s in set(b.index)]
    types = [c for c in a.columns if c in set(b.columns)]
    av = a.loc[samples, types].to_numpy(float).ravel()
    bv = b.loc[samples, types].to_numpy(float).ravel()
    if np.std(av) < 1e-12 or np.std(bv) < 1e-12:
        return {"pearson": float("nan"), "spearman": float("nan")}
    return {"pearson": float(pearsonr(av, bv)[0]),
            "spearman": float(spearmanr(av, bv)[0])}


def main(argv=None) -> int:
    ap = argparse.ArgumentParser(description=__doc__,
                                 formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--run-real-data", action="store_true")
    args = ap.parse_args(argv)

    if not H.real_data_enabled(args.run_real_data):
        print(f"Refusing to run without --run-real-data (or {H.REAL_DATA_ENV}=1).",
              file=sys.stderr)
        return 2
    if not COUNTS_TSV.exists():
        print(f"Missing {COUNTS_TSV}.  Run 00b_download_tcga_tnbc_bulk.py first.",
              file=sys.stderr)
        return 2
    if not H.SAVED_REFERENCE_DIR.exists():
        print(f"Missing reference {H.SAVED_REFERENCE_DIR}.  Run 01_prepare_reference.py.",
              file=sys.stderr)
        return 2

    from tissueresolve.results import ReferenceSignature
    from tissueresolve.api import deconv_bulk

    OUT_DIR.mkdir(parents=True, exist_ok=True)
    ref = ReferenceSignature.load(H.SAVED_REFERENCE_DIR)
    bulk_raw = pd.read_csv(COUNTS_TSV, sep="\t", index_col=0, comment="#")
    bulk, orient = H.orient_bulk_genes_by_samples(bulk_raw, ref.gene_names)
    overlap = H.gene_overlap(bulk.index, ref.gene_names)
    print(f"TCGA TNBC bulk: {bulk.shape[1]} samples × {bulk.shape[0]} genes; "
          f"gene overlap with reference: {overlap['n_shared']} shared.")

    results, bench_rows = {}, []
    for solver in SOLVERS:
        t0 = time.perf_counter()
        with warnings.catch_warnings():
            warnings.simplefilter("ignore")
            res = deconv_bulk(bulk, ref, solver=solver, resolution_mode="none")
        rt = time.perf_counter() - t0
        results[solver] = res
        recon = _recon_pearson(bulk, res.deconv.proportions, ref)
        bench_rows.append({"method": f"TissueResolve_{solver}", "solver": solver,
                           "runtime_seconds": round(rt, 3),
                           "mean_recon_pearson": round(recon, 4)})
        print(f"  {solver:14s} {rt:6.2f}s  mean reconstruction Pearson={recon:.3f}")

    auto = results["auto"]
    auto_props = auto.deconv.proportions
    # concordance of each method vs auto (no ground truth → method agreement)
    for row in bench_rows:
        conc = _flatten_concordance(auto_props, results[row["solver"]].deconv.proportions)
        row["concordance_vs_auto_pearson"] = round(conc["pearson"], 4)
        row["concordance_vs_auto_spearman"] = round(conc["spearman"], 4)

    # --- outputs ---
    H.write_tsv(auto_props.round(4), OUT_DIR / "tcga_bulk_estimated_proportions.tsv",
                comment=["TissueResolve auto-solver; RNA-derived mRNA proportions "
                         "(NOT absolute cell fractions). Real TCGA TNBC bulk — no "
                         "ground truth."])
    H.write_tsv(pd.DataFrame(bench_rows).set_index("method"),
                OUT_DIR / "tcga_bulk_method_benchmark.tsv",
                comment=["No-ground-truth benchmark: runtime, mean reconstruction "
                         "Pearson (observed vs reference-reconstructed CPM, a "
                         "self-consistency check), and cross-method concordance vs "
                         "the auto solver. NOT accuracy."])
    mean_comp = auto_props.mean(0).sort_values(ascending=False)
    H.write_tsv(mean_comp.to_frame("mean_proportion").round(4),
                OUT_DIR / "tcga_bulk_mean_composition.tsv")
    dom = auto_props.idxmax(axis=1).value_counts()
    H.write_tsv(dom.to_frame("n_samples_dominant"),
                OUT_DIR / "tcga_bulk_dominant_type_counts.tsv")

    import json
    warns = [
        {"severity": "info", "category": "estimate_type",
         "message": "Bulk estimates are RNA-derived mRNA proportions, not absolute "
                    "cell fractions."},
        {"severity": "info", "category": "ground_truth",
         "message": "Real TCGA TNBC bulk has no ground-truth cell proportions; "
                    "this is a concordance/QC/runtime benchmark, not accuracy."},
        {"severity": "info", "category": "gene_overlap",
         "message": f"{overlap['n_shared']} genes shared between TCGA bulk and the "
                    f"reference (of {overlap['n_reference']} reference genes)."},
    ]
    (OUT_DIR / "warnings.json").write_text(json.dumps(warns, indent=2), encoding="utf-8")

    from tissueresolve.report import methods_text as mt
    from tissueresolve.report import orchestration
    (OUT_DIR / "methods.txt").write_text(mt.compose_bulk_methods(auto), encoding="utf-8")
    try:
        orchestration.generate_report("bulk", auto, out=OUT_DIR / "report.html")
    except Exception as exc:  # noqa: BLE001
        print(f"  warning: report generation failed: {exc}", file=sys.stderr)

    print("\nMean composition (top 8):")
    print(mean_comp.head(8).round(3).to_string())
    print(f"\nOutputs → {OUT_DIR}/")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
