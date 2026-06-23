#!/usr/bin/env python
"""TCGA bulk prediction sparsity / complexity audit (PART 1 + PART 17).

Diagnostic only — does NOT modify any core algorithm.  Re-runs TissueResolve on
the real TCGA-BRCA TNBC bulk cohort across every internal solver and the
hierarchical path, preserves each method's raw proportion output, and computes a
composition-complexity panel per sample/method so we can determine *where* the
"only 4-5 populations per sample" observation comes from:

    solver behaviour | post-processing | hierarchy gating | report display | true concentration

Outputs (under benchmarks/outputs/)
-----------------------------------
- tcga_prediction_complexity_audit.tsv   per (method, sample) complexity metrics
- tcga_method_complexity_summary.tsv     per-method mean/median summary
- raw/<method>_proportions.tsv           preserved raw method output
- raw/hierarchical_*.tsv                  family / conditional / fine / combined / unresolved

Usage:  python benchmarks/audit/tcga_sparsity_audit.py --run-real-data
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
sys.path.insert(0, str(REPO / "examples" / "real_breast_cancer" / "scripts"))
import _harness as H  # noqa: E402

COUNTS_TSV = H.DATA_DIR / "bulk_tcga_tnbc" / "tcga_tnbc_counts.tsv"
OUT_DIR = REPO / "benchmarks" / "outputs"
RAW_DIR = OUT_DIR / "raw"

FLAT_SOLVERS = ["nnls", "weighted_nnls", "marker_nnls", "ridge_nnls", "auto"]
THRESHOLDS = [0.0, 1e-4, 1e-3, 5e-3, 1e-2]
SEED = 0


# ----------------------------------------------------------------------------- metrics
def _gini(x: np.ndarray) -> float:
    """Gini coefficient of a non-negative composition vector."""
    x = np.sort(np.asarray(x, float))
    n = x.size
    s = x.sum()
    if n == 0 or s <= 0:
        return float("nan")
    idx = np.arange(1, n + 1)
    return float((np.sum((2 * idx - n - 1) * x)) / (n * s))


def _entropy_nats(p: np.ndarray) -> float:
    p = np.asarray(p, float)
    p = p[p > 0]
    s = p.sum()
    if s <= 0:
        return float("nan")
    p = p / s
    return float(-np.sum(p * np.log(p)))


def complexity_row(vec: np.ndarray, n_total: int) -> dict:
    """Composition-complexity metrics for one sample's proportion vector."""
    v = np.asarray(vec, float)
    v = np.where(np.isfinite(v), v, 0.0)
    out = {"n_reference_populations": n_total, "mass_sum": float(v.sum())}
    for t in THRESHOLDS:
        key = "n_gt_0" if t == 0.0 else f"n_gt_{t:g}"
        out[key] = int(np.sum(v > t))
    s = v.sum()
    p = v / s if s > 0 else v
    H_ = _entropy_nats(p)
    out["shannon_entropy_nats"] = H_
    out["effective_n_populations"] = float(np.exp(H_)) if np.isfinite(H_) else float("nan")
    out["gini"] = _gini(v)
    out["dominant_fraction"] = float(np.max(p)) if s > 0 else float("nan")
    srt = np.sort(p)[::-1]
    for k in (1, 3, 5):
        out[f"top{k}_cumulative"] = float(srt[:k].sum())
    return out


# ----------------------------------------------------------------------------- run
def main(argv=None) -> int:
    ap = argparse.ArgumentParser(description=__doc__,
                                 formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--run-real-data", action="store_true")
    ap.add_argument("--solvers", nargs="*", default=FLAT_SOLVERS)
    args = ap.parse_args(argv)

    if not H.real_data_enabled(args.run_real_data):
        print(f"Refusing without --run-real-data (or {H.REAL_DATA_ENV}=1).", file=sys.stderr)
        return 2
    if not COUNTS_TSV.exists():
        print(f"Missing {COUNTS_TSV}. Run 00b_download_tcga_tnbc_bulk.py first.", file=sys.stderr)
        return 2
    if not H.SAVED_REFERENCE_DIR.exists():
        print(f"Missing reference {H.SAVED_REFERENCE_DIR}.", file=sys.stderr)
        return 2

    np.random.seed(SEED)
    from tissueresolve.results import ReferenceSignature
    from tissueresolve.api import deconv_bulk

    RAW_DIR.mkdir(parents=True, exist_ok=True)
    ref = ReferenceSignature.load(H.SAVED_REFERENCE_DIR)
    n_ref_types = len(list(ref.cell_types))
    bulk_raw = pd.read_csv(COUNTS_TSV, sep="\t", index_col=0, comment="#")
    bulk, orient = H.orient_bulk_genes_by_samples(bulk_raw, ref.gene_names)
    overlap = H.gene_overlap(bulk.index, ref.gene_names)
    print(f"TCGA TNBC: {bulk.shape[1]} samples x {bulk.shape[0]} genes; "
          f"{overlap['n_shared']} shared with reference ({n_ref_types} populations).")

    rows = []  # per (method, sample) complexity
    summaries = []

    # ---- FLAT solvers (resolution_mode='none' = raw normalized solver output) ----
    for solver in args.solvers:
        t0 = time.perf_counter()
        try:
            with warnings.catch_warnings():
                warnings.simplefilter("ignore")
                res = deconv_bulk(bulk, ref, solver=solver, resolution_mode="none")
            props = res.deconv.proportions
        except Exception as exc:  # noqa: BLE001
            print(f"  FLAT {solver:14s} FAILED: {exc}")
            summaries.append({"method": f"flat_{solver}", "status": f"failed: {exc}"})
            continue
        rt = time.perf_counter() - t0
        method = f"flat_{solver}"
        props.round(6).to_csv(RAW_DIR / f"{method}_proportions.tsv", sep="\t")
        for samp, vec in props.iterrows():
            r = complexity_row(vec.to_numpy(float), n_ref_types)
            r.update({"method": method, "sample": samp, "level": "flat",
                      "unresolved_mass": 0.0})
            rows.append(r)
        sub = pd.DataFrame([r for r in rows if r["method"] == method])
        summaries.append(_summary(method, "flat", sub, rt))
        print(f"  FLAT {solver:14s} {rt:6.1f}s  mean n_gt_0={sub['n_gt_0'].mean():.1f} "
              f"mean eff_N={sub['effective_n_populations'].mean():.1f}")

    # ---- HIERARCHICAL path ----
    t0 = time.perf_counter()
    try:
        with warnings.catch_warnings():
            warnings.simplefilter("ignore")
            hres = deconv_bulk(bulk, ref, solver="auto", resolution_mode="hierarchical")
        est = hres.estimates
        fam = est.family_proportions
        comb = est.combined_fine
        fine = est.fine_proportions
        cond = est.conditional_proportions
        unres = est.unresolved_mass
        fam.round(6).to_csv(RAW_DIR / "hierarchical_family_proportions.tsv", sep="\t")
        comb.round(6).to_csv(RAW_DIR / "hierarchical_combined_proportions.tsv", sep="\t")
        fine.round(6).to_csv(RAW_DIR / "hierarchical_fine_proportions.tsv", sep="\t")
        cond.round(6).to_csv(RAW_DIR / "hierarchical_conditional_proportions.tsv", sep="\t")
        if unres is not None:
            pd.DataFrame(unres).round(6).to_csv(RAW_DIR / "hierarchical_unresolved_mass.tsv", sep="\t")
        rt = time.perf_counter() - t0
        n_fam = fam.shape[1]
        n_fine = fine.shape[1]
        # unresolved mass per sample (sum of unresolved_* columns in combined, or est.unresolved_mass)
        if isinstance(unres, pd.DataFrame):
            unres_per = unres.sum(axis=1)
        elif isinstance(unres, pd.Series):
            unres_per = unres
        else:
            ucols = [c for c in comb.columns if str(c).startswith("unresolved")]
            unres_per = comb[ucols].sum(axis=1) if ucols else pd.Series(0.0, index=comb.index)

        for level, df, ntot in [("hier_family", fam, n_fam),
                                ("hier_fine", fine, n_fine),
                                ("hier_combined", comb, comb.shape[1])]:
            for samp, vec in df.iterrows():
                r = complexity_row(vec.to_numpy(float), ntot)
                r.update({"method": f"hierarchical_{level}", "sample": samp,
                          "level": level,
                          "unresolved_mass": float(unres_per.get(samp, np.nan))})
                rows.append(r)
            sub = pd.DataFrame([rr for rr in rows if rr["method"] == f"hierarchical_{level}"])
            summaries.append(_summary(f"hierarchical_{level}", level, sub, rt))
            print(f"  HIER {level:14s}        mean n_gt_0={sub['n_gt_0'].mean():.2f}/{ntot} "
                  f"mean eff_N={sub['effective_n_populations'].mean():.2f}")
        # resolvability summary
        resolv = getattr(est, "resolvability", None)
        unres_fams = getattr(est, "unresolved_families", None)
        meta = {"n_families": int(n_fam), "n_fine": int(n_fine),
                "unresolved_families": list(unres_fams) if unres_fams is not None else None,
                "mean_unresolved_mass": float(unres_per.mean()),
                "resolvability": _jsonable(resolv)}
        (RAW_DIR / "hierarchical_resolvability.json").write_text(
            json.dumps(meta, indent=2, default=str), encoding="utf-8")
        print(f"  HIER unresolved families: {meta['unresolved_families']}")
        print(f"  HIER mean unresolved mass: {meta['mean_unresolved_mass']:.3f}")
    except Exception as exc:  # noqa: BLE001
        import traceback
        traceback.print_exc()
        print(f"  HIERARCHICAL FAILED: {exc}")
        summaries.append({"method": "hierarchical", "status": f"failed: {exc}"})

    # ---- write tables ----
    audit = pd.DataFrame(rows)
    front = ["method", "sample", "level", "n_reference_populations"]
    cols = front + [c for c in audit.columns if c not in front]
    audit = audit[cols]
    OUT_DIR.mkdir(parents=True, exist_ok=True)
    _write_tsv(audit, OUT_DIR / "tcga_prediction_complexity_audit.tsv",
               "PART 17 TCGA per-sample composition-complexity audit. Diagnostic only.")
    summ = pd.DataFrame(summaries)
    _write_tsv(summ, OUT_DIR / "tcga_method_complexity_summary.tsv",
               "Per-method complexity summary across TCGA samples. Diagnostic only.")
    print(f"\nWrote {OUT_DIR/'tcga_prediction_complexity_audit.tsv'}")
    print(f"Wrote {OUT_DIR/'tcga_method_complexity_summary.tsv'}")
    print("\n=== METHOD SUMMARY ===")
    show = [c for c in ["method", "mean_n_gt_0", "mean_n_gt_0.01", "mean_effective_n",
                        "mean_dominant_fraction", "mean_top5_cumulative",
                        "mean_unresolved_mass"] if c in summ.columns]
    print(summ[show].to_string(index=False))
    return 0


def _summary(method, level, sub, rt):
    return {
        "method": method, "level": level, "runtime_seconds": round(rt, 1),
        "n_samples": int(len(sub)),
        "mean_n_reference_populations": float(sub["n_reference_populations"].mean()),
        "mean_n_gt_0": float(sub["n_gt_0"].mean()),
        "mean_n_gt_0.0001": float(sub["n_gt_0.0001"].mean()),
        "mean_n_gt_0.001": float(sub["n_gt_0.001"].mean()),
        "mean_n_gt_0.005": float(sub["n_gt_0.005"].mean()),
        "mean_n_gt_0.01": float(sub["n_gt_0.01"].mean()),
        "mean_effective_n": float(sub["effective_n_populations"].mean()),
        "median_effective_n": float(sub["effective_n_populations"].median()),
        "mean_shannon_entropy": float(sub["shannon_entropy_nats"].mean()),
        "mean_gini": float(sub["gini"].mean()),
        "mean_dominant_fraction": float(sub["dominant_fraction"].mean()),
        "mean_top1_cumulative": float(sub["top1_cumulative"].mean()),
        "mean_top3_cumulative": float(sub["top3_cumulative"].mean()),
        "mean_top5_cumulative": float(sub["top5_cumulative"].mean()),
        "mean_unresolved_mass": float(sub["unresolved_mass"].mean()),
    }


def _jsonable(obj):
    if obj is None:
        return None
    if isinstance(obj, pd.DataFrame):
        return obj.to_dict(orient="index")
    if isinstance(obj, pd.Series):
        return obj.to_dict()
    try:
        json.dumps(obj)
        return obj
    except TypeError:
        return str(obj)


def _write_tsv(df, path, comment):
    with open(path, "w", encoding="utf-8") as fh:
        fh.write(f"# {comment}\n")
        df.to_csv(fh, sep="\t", index=False)


if __name__ == "__main__":
    raise SystemExit(main())
