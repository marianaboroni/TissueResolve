#!/usr/bin/env python
"""Multi-split held-out-donor bulk benchmark with statistics (PART 14).

Hardens the single-split results: repeats the donor-disjoint benchmark across
several random donor splits, treats each (split, scenario) as a paired replicate,
and reports bootstrap 95% CIs, paired Wilcoxon tests, and bootstrap rank
stability.  Also adds a FAIR hierarchical re-score (resolvable-only fine metrics)
so honest abstention is not counted as a wrong fine prediction.

TissueResolve-only, no external installs, deterministic.

Outputs (benchmarks/outputs/holdout_bulk_multisplit/)
- metrics_long.tsv         (split, scenario, method, metric, value)
- metrics_ci.tsv           (method, metric, point, lo, hi, n)
- pairwise_wilcoxon.tsv    (metric, method_a, method_b, statistic, p, effect, n)
- rank_stability.tsv       (method, mean_rank, prob_best) for fine Pearson
- manifest.json

Usage:  python benchmarks/bulk/run_holdout_bulk_multisplit.py --run-real-data
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
from benchmarks.shared import stats as ST  # noqa: E402

OUT = REPO / "benchmarks" / "outputs" / "holdout_bulk_multisplit"
HIERARCHY_TSV = REPO / "examples" / "real_breast_cancer" / "config" / \
    "breast_cancer_cell_type_hierarchy.tsv"

FLAT_SOLVERS = ["nnls", "weighted_nnls", "ridge_nnls", "auto"]
SCENARIOS = ["balanced", "imbalanced", "rare", "similar_subtypes", "missing_population"]
SCEN_OFFSET = {s: i * 100 for i, s in enumerate(SCENARIOS)}
N_SPLITS = 5
N_SAMPLES = 12
CELLS_PER_SAMPLE = 600
MIN_CELLS = 30


def _align_fine(truth, pred):
    cols = [c for c in truth.columns if not str(c).startswith("unresolved_")]
    return truth[cols], pred.reindex(index=truth.index, columns=cols).fillna(0.0)


def run_one_split(adata, raw_map, split_seed: int) -> list[dict]:
    from tissueresolve.api import deconv_bulk
    from tissueresolve.reference.hierarchy import build_cell_type_hierarchy

    ref_donors, query_donors = SH.split_donors(adata, "donor_id",
                                               ref_frac=0.5, seed=split_seed)
    ref_mask = adata.obs["donor_id"].astype(str).isin(set(ref_donors)).to_numpy()
    with warnings.catch_warnings():
        warnings.simplefilter("ignore")
        ref = H.prepare_reference(adata[ref_mask].copy(), min_cells=MIN_CELLS,
                                  estimate_overdispersion=True).reference
    ref_types = [str(c) for c in ref.cell_types]
    mapping = build_cell_type_hierarchy(ref_types, raw_map)

    fam_members: dict[str, list] = {}
    for t in ref_types:
        fam_members.setdefault(mapping.get(t, t), []).append(t)
    sim_fam = max(fam_members, key=lambda f: len(fam_members[f]))
    similar_pair = tuple(fam_members[sim_fam][:2])
    rare_type = next((t for t in ["regulatory T cell", "natural killer cell"]
                      if t in ref_types), ref_types[-1])
    missing_type = next((t for t in ref_types if mapping.get(t) == "Epithelial"),
                        ref_types[0])

    ref_missing = None
    if "missing_population" in SCENARIOS:
        mm = ref_mask & (adata.obs["cell_type"].astype(str) != missing_type).to_numpy()
        with warnings.catch_warnings():
            warnings.simplefilter("ignore")
            ref_missing = H.prepare_reference(adata[mm].copy(), min_cells=MIN_CELLS,
                                              estimate_overdispersion=True).reference

    rows = []
    for scen in SCENARIOS:
        kw, ref_scen = {}, ref
        if scen == "rare":
            kw = {"rare_type": rare_type, "rare_level": 0.01}
        elif scen == "similar_subtypes":
            kw = {"similar_pair": similar_pair}
        elif scen == "missing_population":
            kw, ref_scen = {"missing_type": missing_type}, ref_missing
        tseed = 10 + SCEN_OFFSET[scen] + split_seed * 1000
        targets = SH.build_target_proportions(ref_types, N_SAMPLES, scen, seed=tseed, **kw)
        ds = SH.realize_pseudobulk(adata, targets, celltype_col="cell_type",
                                   donor_col="donor_id", query_donors=query_donors,
                                   seed=split_seed * 7 + 1, cells_per_sample=CELLS_PER_SAMPLE)
        truth = ds.true_mrna_proportions

        runs = {}
        for solver in FLAT_SOLVERS:
            try:
                with warnings.catch_warnings():
                    warnings.simplefilter("ignore")
                    runs[f"flat_{solver}"] = deconv_bulk(
                        ds.counts, ref_scen, solver=solver,
                        resolution_mode="none").deconv.proportions
            except Exception:  # noqa: BLE001
                pass
        hres = None
        try:
            with warnings.catch_warnings():
                warnings.simplefilter("ignore")
                hres = deconv_bulk(ds.counts, ref_scen, solver="auto",
                                   resolution_mode="hierarchical", hierarchy_mapping=mapping)
            runs["hierarchical"] = hres.estimates.combined_fine
        except Exception:  # noqa: BLE001
            pass

        def _emit(method, fine_acc, fine_comp, pred_fine, truth_fine):
            tfam = M.aggregate_to_families(truth_fine, mapping)
            pfam = M.aggregate_to_families(pred_fine, mapping)
            bacc = M.accuracy_metrics(tfam, pfam)
            cxt = M.composition_complexity(truth_fine)["effective_n_populations"].mean()
            cxp = M.composition_complexity(pred_fine)["effective_n_populations"].mean()
            base = {"split": split_seed, "scenario": scen, "method": method}
            for metric, val in [
                ("fine_pearson", fine_acc["pearson"]), ("fine_rmse", fine_acc["rmse"]),
                ("fine_jsd", fine_comp["jsd_mean"]), ("fine_ccc", fine_comp["ccc"]),
                ("broad_pearson", bacc["pearson"]),
                ("dominant_family_accuracy", M.dominant_accuracy(tfam, pfam)),
                ("complexity_abs_error", abs(cxp - cxt))]:
                rows.append({**base, "metric": metric, "value": float(val)})

        for method, pred in runs.items():
            tf, pf = _align_fine(truth, pred)
            _emit(method, M.accuracy_metrics(tf, pf), M.compositional_metrics(tf, pf), pf, tf)

        # FAIR hierarchical: fine restricted to families it did NOT abstain on
        if hres is not None:
            unresolved = set(getattr(hres.estimates, "unresolved_families", []) or [])
            fair = M.fine_metrics_resolvable_only(truth, hres.estimates.combined_fine,
                                                  mapping, unresolved)
            if fair.get("n_cell_types", 0) > 0:
                keep = [c for c in truth.columns
                        if mapping.get(c, c) not in unresolved
                        and not str(c).startswith("unresolved_")]
                tf = truth[keep]
                pf = hres.estimates.combined_fine.reindex(index=tf.index, columns=keep).fillna(0.0)
                _emit("hierarchical_resolvable", M.accuracy_metrics(tf, pf),
                      M.compositional_metrics(tf, pf), pf, tf)
    return rows


def main(argv=None) -> int:
    ap = argparse.ArgumentParser(description=__doc__,
                                 formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--run-real-data", action="store_true")
    ap.add_argument("--n-splits", type=int, default=N_SPLITS)
    args = ap.parse_args(argv)
    if not H.real_data_enabled(args.run_real_data):
        print(f"Refusing without --run-real-data (or {H.REAL_DATA_ENV}=1).", file=sys.stderr)
        return 2

    import anndata as ad
    from tissueresolve.reference.hierarchy import load_hierarchy_mapping

    adata = ad.read_h5ad(H.REFERENCE_H5AD)
    if "feature_name" in adata.var.columns:
        fn = adata.var["feature_name"].astype(str)
        if fn.nunique() == adata.n_vars:
            adata.var_names = fn.to_numpy()
    raw_map = load_hierarchy_mapping(HIERARCHY_TSV)

    all_rows = []
    for s in range(args.n_splits):
        print(f"[split {s}] ...", flush=True)
        all_rows.extend(run_one_split(adata, raw_map, s))
    long = pd.DataFrame(all_rows)
    OUT.mkdir(parents=True, exist_ok=True)

    def _w(df, name, comment):
        with open(OUT / name, "w", encoding="utf-8") as fh:
            fh.write(f"# {comment}\n")
            df.to_csv(fh, sep="\t", index=False)
    _w(long, "metrics_long.tsv", "Per (split, scenario, method) scalar metrics.")

    # bootstrap CIs per (method, metric) over split×scenario replicates
    ci_rows = []
    for (method, metric), g in long.groupby(["method", "metric"]):
        c = ST.bootstrap_ci(g["value"].to_numpy(), seed=0)
        ci_rows.append({"method": method, "metric": metric, **c})
    ci = pd.DataFrame(ci_rows)
    _w(ci, "metrics_ci.tsv", "Bootstrap 95% CIs per method/metric over split x scenario replicates.")

    # paired Wilcoxon on fine_pearson and complexity_abs_error
    wil_rows = []
    for metric in ["fine_pearson", "complexity_abs_error", "broad_pearson"]:
        piv = long[long.metric == metric].pivot_table(
            index=["split", "scenario"], columns="method", values="value")
        methods = list(piv.columns)
        for i in range(len(methods)):
            for j in range(i + 1, len(methods)):
                a, b = methods[i], methods[j]
                r = ST.paired_wilcoxon(piv[a].to_numpy(), piv[b].to_numpy())
                wil_rows.append({"metric": metric, "method_a": a, "method_b": b, **r})
    _w(pd.DataFrame(wil_rows), "pairwise_wilcoxon.tsv",
       "Paired Wilcoxon signed-rank across split x scenario replicates.")

    # rank stability on fine_pearson
    piv = long[long.metric == "fine_pearson"].pivot_table(
        index=["split", "scenario"], columns="method", values="value")
    rs = ST.rank_stability(piv, higher_is_better=True, seed=0)
    _w(rs.reset_index(), "rank_stability.tsv",
       "Bootstrap rank distribution on fine Pearson (higher better).")

    (OUT / "manifest.json").write_text(json.dumps({
        "n_splits": args.n_splits, "scenarios": SCENARIOS, "n_samples": N_SAMPLES,
        "cells_per_sample": CELLS_PER_SAMPLE, "min_cells": MIN_CELLS,
        "replicate_unit": "split x scenario", "ground_truth": "mRNA-proportion",
    }, indent=2), encoding="utf-8")

    print(f"\nWrote -> {OUT}/")
    print("\n=== fine Pearson (mean [95% CI]) ===")
    fp = ci[ci.metric == "fine_pearson"].sort_values("point", ascending=False)
    for _, r in fp.iterrows():
        print(f"  {r['method']:24s} {r['point']:.3f} [{r['lo']:.3f}, {r['hi']:.3f}]  n={int(r['n'])}")
    print("\n=== rank stability (fine Pearson) ===")
    print(rs.round(3).to_string())
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
