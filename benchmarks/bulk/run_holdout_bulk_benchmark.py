#!/usr/bin/env python
"""Phase-2 ground-truth bulk benchmark — TissueResolve only (no external installs).

End-to-end, donor-disjoint, count-level. For each scenario:

1. Split the 41 single-cell donors into disjoint *reference* / *query* sets.
2. Build the TissueResolve reference signature from **reference donors only**.
3. Generate count-level pseudobulk from **query donors** (held out) with known
   cell-fraction AND mRNA-proportion ground truth.
4. Run TissueResolve flat solvers (nnls/weighted/ridge/auto) + hierarchical.
5. Score against the **mRNA-proportion truth** (what an RNA deconvolver estimates,
   CLAUDE.md bulk rule 1): fine-level + broad-level accuracy & compositional
   distances, plus composition-complexity preservation and rare-pop detection.

Outputs (benchmarks/outputs/holdout_bulk/)
- bulk_fine_metrics.tsv, bulk_broad_metrics.tsv, bulk_complexity_metrics.tsv
- bulk_rare_population_metrics.tsv
- raw/<scenario>__<method>_pred.tsv, raw/<scenario>_truth_mrna.tsv
- manifest.json  (donors, versions, params, reference summary)

This is validation tooling; it does NOT modify any core algorithm.

Usage:  python benchmarks/bulk/run_holdout_bulk_benchmark.py --run-real-data
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
sys.path.insert(0, str(REPO))
import _harness as H  # noqa: E402
from benchmarks.shared import synthetic_holdout as SH  # noqa: E402
from benchmarks.shared import metrics as M  # noqa: E402

OUT = REPO / "benchmarks" / "outputs" / "holdout_bulk"
RAW = OUT / "raw"
HIERARCHY_TSV = REPO / "examples" / "real_breast_cancer" / "config" / \
    "breast_cancer_cell_type_hierarchy.tsv"

FLAT_SOLVERS = ["nnls", "weighted_nnls", "ridge_nnls", "auto"]
SEED = 0
N_SAMPLES = 12
CELLS_PER_SAMPLE = 600
MIN_CELLS = 30  # lower than shipped (50) since only ~half the donors are used


def _versions() -> dict:
    import scipy, sklearn  # noqa
    import tissueresolve
    return {"tissueresolve": getattr(tissueresolve, "__version__", "?"),
            "numpy": np.__version__, "pandas": pd.__version__,
            "scipy": scipy.__version__}


def _complexity_pair(truth: pd.DataFrame, pred: pd.DataFrame) -> dict:
    """Composition-complexity preservation: predicted vs true (mean over samples)."""
    ct, cp = M.composition_complexity(truth), M.composition_complexity(pred)
    return {
        "true_effective_n": float(ct["effective_n_populations"].mean()),
        "pred_effective_n": float(cp["effective_n_populations"].mean()),
        "effective_n_abs_error": float(np.mean(np.abs(
            cp["effective_n_populations"].to_numpy() - ct["effective_n_populations"].to_numpy()))),
        "true_n_gt_0.01": float(ct["n_gt_0.01"].mean()),
        "pred_n_gt_0.01": float(cp["n_gt_0.01"].mean()),
        "true_gini": float(ct["gini"].mean()),
        "pred_gini": float(cp["gini"].mean()),
        "true_dominant": float(ct["dominant_fraction"].mean()),
        "pred_dominant": float(cp["dominant_fraction"].mean()),
    }


def _align_fine(truth: pd.DataFrame, pred: pd.DataFrame) -> tuple:
    """Restrict to common columns; missing pred columns filled 0 (fair to method)."""
    cols = [c for c in truth.columns if not str(c).startswith("unresolved_")]
    p = pred.reindex(index=truth.index, columns=cols).fillna(0.0)
    return truth[cols], p


def main(argv=None) -> int:
    ap = argparse.ArgumentParser(description=__doc__,
                                 formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--run-real-data", action="store_true")
    ap.add_argument("--scenarios", nargs="*",
                    default=["balanced", "imbalanced", "rare", "similar_subtypes",
                             "missing_population"])
    args = ap.parse_args(argv)
    if not H.real_data_enabled(args.run_real_data):
        print(f"Refusing without --run-real-data (or {H.REAL_DATA_ENV}=1).", file=sys.stderr)
        return 2
    if not H.REFERENCE_H5AD.exists():
        print(f"Missing reference h5ad {H.REFERENCE_H5AD}.", file=sys.stderr)
        return 2

    import anndata as ad
    from tissueresolve.api import deconv_bulk
    from tissueresolve.reference.hierarchy import (
        load_hierarchy_mapping, build_cell_type_hierarchy)

    RAW.mkdir(parents=True, exist_ok=True)
    print(f"Loading {H.REFERENCE_H5AD} ...")
    adata = ad.read_h5ad(H.REFERENCE_H5AD)
    # Harmonise gene ids to symbols so pseudobulk counts and the reference (which
    # prepare_reference builds on feature_name) share the same gene namespace.
    # CELLxGENE var_names are numeric soma_joinid; feature_name holds symbols.
    if "feature_name" in adata.var.columns:
        fn = adata.var["feature_name"].astype(str)
        if fn.nunique() != adata.n_vars:
            raise SystemExit("duplicate gene symbols in feature_name; this runner "
                             "assumes unique symbols (none in this dataset).")
        adata.var_names = fn.to_numpy()
    print(f"  {adata.n_obs} cells x {adata.n_vars} genes (genes={list(adata.var_names[:3])}...)")

    # one shared donor split (seed 0): reference donors vs disjoint query donors
    ref_donors, query_donors = SH.split_donors(adata, "donor_id", ref_frac=0.5, seed=SEED)
    assert set(ref_donors).isdisjoint(query_donors)
    print(f"  reference donors: {len(ref_donors)}  query donors: {len(query_donors)} (disjoint)")

    # build the shared reference from REFERENCE donors only (held-out query)
    ref_mask = adata.obs["donor_id"].astype(str).isin(set(ref_donors)).to_numpy()
    t0 = time.perf_counter()
    with warnings.catch_warnings():
        warnings.simplefilter("ignore")
        prep = H.prepare_reference(adata[ref_mask].copy(), min_cells=MIN_CELLS,
                                   estimate_overdispersion=True)
    ref = prep.reference
    ref_types = [str(c) for c in ref.cell_types]
    print(f"  reference: {ref.n_genes} genes x {len(ref_types)} cell types "
          f"({time.perf_counter()-t0:.1f}s)")

    raw_map = load_hierarchy_mapping(HIERARCHY_TSV)
    mapping = build_cell_type_hierarchy(ref_types, raw_map)
    families = sorted(set(mapping.values()))

    # scenario-specific params drawn from the reference's own types
    rare_type = next((t for t in ["regulatory T cell", "natural killer cell",
                                   "conventional dendritic cell"] if t in ref_types),
                     ref_types[-1])
    # a same-family similar pair
    fam_members: dict[str, list] = {}
    for t in ref_types:
        fam_members.setdefault(mapping.get(t, t), []).append(t)
    sim_fam = max(fam_members, key=lambda f: len(fam_members[f]))
    similar_pair = tuple(fam_members[sim_fam][:2])
    missing_type = next((t for t in ref_types if mapping.get(t) == "Epithelial"), ref_types[0])

    # deterministic per-scenario target seed (NOT Python hash(), which varies per process)
    SCEN_SEED = {"balanced": 11, "imbalanced": 22, "rare": 33,
                 "similar_subtypes": 44, "missing_population": 55}

    # dedicated reference EXCLUDING the missing type's cells (built once, if needed)
    ref_missing = None
    if "missing_population" in args.scenarios:
        miss_mask = ref_mask & (adata.obs["cell_type"].astype(str) != missing_type).to_numpy()
        with warnings.catch_warnings():
            warnings.simplefilter("ignore")
            ref_missing = H.prepare_reference(adata[miss_mask].copy(), min_cells=MIN_CELLS,
                                              estimate_overdispersion=True).reference
        print(f"  missing-pop reference excludes {missing_type!r}: "
              f"{len(list(ref_missing.cell_types))} types")

    fine_rows, broad_rows, cx_rows, rare_rows = [], [], [], []

    for scen in args.scenarios:
        kw = {}
        ref_scen = ref
        if scen == "rare":
            kw = {"rare_type": rare_type, "rare_level": 0.01}
        elif scen == "similar_subtypes":
            kw = {"similar_pair": similar_pair}
        elif scen == "missing_population":
            ref_scen = ref_missing  # reference lacks missing_type; mixture still contains it
            kw = {"missing_type": missing_type}

        targets = SH.build_target_proportions(ref_types, N_SAMPLES, scen,
                                               seed=SCEN_SEED.get(scen, SEED), **kw)
        ds = SH.realize_pseudobulk(adata, targets, celltype_col="cell_type",
                                   donor_col="donor_id", query_donors=query_donors,
                                   seed=SEED, cells_per_sample=CELLS_PER_SAMPLE)
        truth = ds.true_mrna_proportions  # the RNA-based target
        truth.to_csv(RAW / f"{scen}_truth_mrna.tsv", sep="\t")
        ds.true_cell_fractions.to_csv(RAW / f"{scen}_truth_cellfrac.tsv", sep="\t")
        counts = ds.counts
        print(f"\n[{scen}] {counts.shape[1]} samples, {counts.shape[0]} genes; "
              f"true eff-N={M.composition_complexity(truth)['effective_n_populations'].mean():.1f}")

        runs = {}
        for solver in FLAT_SOLVERS:
            try:
                with warnings.catch_warnings():
                    warnings.simplefilter("ignore")
                    res = deconv_bulk(counts, ref_scen, solver=solver, resolution_mode="none")
                runs[f"flat_{solver}"] = res.deconv.proportions
            except Exception as exc:  # noqa: BLE001
                print(f"    flat_{solver} FAILED: {exc}")
        try:
            with warnings.catch_warnings():
                warnings.simplefilter("ignore")
                hres = deconv_bulk(counts, ref_scen, solver="auto",
                                   resolution_mode="hierarchical", hierarchy_mapping=mapping)
            runs["hierarchical"] = hres.estimates.combined_fine
        except Exception as exc:  # noqa: BLE001
            print(f"    hierarchical FAILED: {exc}")

        for method, pred in runs.items():
            pred.to_csv(RAW / f"{scen}__{method}_pred.tsv", sep="\t")
            tf, pf = _align_fine(truth, pred)
            acc = M.accuracy_metrics(tf, pf)
            comp = M.compositional_metrics(tf, pf)
            fine_rows.append({"scenario": scen, "method": method, **acc,
                              "jsd": comp["jsd_mean"], "tv": comp["tv_mean"],
                              "aitchison": comp["aitchison_mean"], "ccc": comp["ccc"]})
            # broad: aggregate truth + pred to families
            tfam = M.aggregate_to_families(tf, mapping)
            pfam = M.aggregate_to_families(pred, mapping)
            bacc = M.accuracy_metrics(tfam, pfam)
            bcomp = M.compositional_metrics(tfam, pfam)
            broad_rows.append({"scenario": scen, "method": method, **bacc,
                               "jsd": bcomp["jsd_mean"], "ccc": bcomp["ccc"],
                               "dominant_family_accuracy": M.dominant_accuracy(tfam, pfam)})
            # complexity preservation
            cx_rows.append({"scenario": scen, "method": method, **_complexity_pair(tf, pf)})
            # rare detection
            if scen == "rare" and rare_type in pf.columns:
                est = pf[rare_type].to_numpy(float)
                tru = tf[rare_type].to_numpy(float)
                rare_rows.append({
                    "scenario": scen, "method": method, "rare_type": rare_type,
                    "true_level": float(np.mean(tru)),
                    "mean_estimated": float(np.mean(est)),
                    "detection_rate_gt_1e-3": float(np.mean(est > 1e-3)),
                    "detection_rate_gt_5e-3": float(np.mean(est > 5e-3)),
                    "abundance_bias": float(np.mean(est - tru))})
            print(f"    {method:18s} fine r={acc['pearson']:.3f} rmse={acc['rmse']:.3f} "
                  f"jsd={comp['jsd_mean']:.3f} | broad r={bacc['pearson']:.3f}")

    # write tables
    OUT.mkdir(parents=True, exist_ok=True)
    def _w(rows, name, comment):
        df = pd.DataFrame(rows)
        with open(OUT / name, "w", encoding="utf-8") as fh:
            fh.write(f"# {comment}\n")
            df.to_csv(fh, sep="\t", index=False)
        return df
    _w(fine_rows, "bulk_fine_metrics.tsv",
       "TissueResolve-only held-out-donor bulk benchmark. Fine-level accuracy vs mRNA-proportion truth.")
    _w(broad_rows, "bulk_broad_metrics.tsv",
       "Broad-family accuracy (truth+pred aggregated to families) vs mRNA-proportion truth.")
    cx = _w(cx_rows, "bulk_complexity_metrics.tsv",
            "Composition-complexity preservation: predicted vs true effective-N / richness / Gini.")
    if rare_rows:
        _w(rare_rows, "bulk_rare_population_metrics.tsv",
           f"Rare-population detection for {rare_type} at true level ~0.01.")

    manifest = {
        "reference_donors": ref_donors, "query_donors": query_donors,
        "donor_disjoint": True, "n_reference_types": len(ref_types),
        "min_cells": MIN_CELLS, "n_samples_per_scenario": N_SAMPLES,
        "cells_per_sample": CELLS_PER_SAMPLE, "seed": SEED,
        "scenarios": args.scenarios, "families": families,
        "rare_type": rare_type, "similar_pair": list(similar_pair),
        "missing_type": missing_type, "versions": _versions(),
        "ground_truth": "mRNA-proportion (true_mrna_proportions); cell-fraction also saved",
    }
    (OUT / "manifest.json").write_text(json.dumps(manifest, indent=2), encoding="utf-8")
    print(f"\nWrote metrics + manifest -> {OUT}/")
    if not cx.empty:
        print("\n=== complexity preservation (pred vs true effective-N) ===")
        print(cx[["scenario", "method", "true_effective_n", "pred_effective_n",
                  "effective_n_abs_error"]].to_string(index=False))
    else:
        print("WARNING: no successful runs — all method×scenario combinations failed.")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
