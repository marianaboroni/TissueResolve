#!/usr/bin/env python
"""PHASE 1 — Hierarchical-failure diagnosis via oracle experiments (DIAGNOSIS ONLY).

Does NOT modify or add any algorithm.  Uses existing public APIs + oracle
injections (true broad mass, oracle/global/within-family gene panels, gating
toggles, spatial-lambda sweep) to localise WHY the hierarchical path underperforms
the flat path on bulk + spatial ground truth.

Bulk oracle experiments (held-out-donor synthetic, broad+fine mRNA truth):
  E1 oracle-broad      : inject TRUE family mass into the broad→fine combine
  E2 oracle-gene-panel : within-family panel = reference within-family DE (oracle
                         for synthetic) vs global panel vs pipeline default
  E3 pre- vs post-gating: all-subtype combine vs gated/unresolved vs no-gating
Spatial oracle experiment:
  E4 smoothing sweep   : lambda_spatial in {0, weak, default, strong}

Outputs:
  benchmarks/outputs/hierarchical_oracle_experiments.tsv

Usage:  python benchmarks/diagnostics/hierarchical_oracle.py --run-real-data
"""
from __future__ import annotations

import argparse
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
from benchmarks.shared import synthetic_spatial as SSP  # noqa: E402
from benchmarks.shared import metrics as M  # noqa: E402
from benchmarks.shared import spatial_metrics as SM  # noqa: E402

OUT = REPO / "benchmarks" / "outputs"
HIER = REPO / "examples" / "real_breast_cancer" / "config" / "breast_cancer_cell_type_hierarchy.tsv"
BULK_SEEDS = [0, 1, 2]
BULK_SCENARIOS = ["imbalanced", "similar_subtypes"]   # fine-relevant scenarios
SPATIAL_SEEDS = [0, 1]
N_SAMPLES = 12
CELLS_PER_SAMPLE = 600
MIN_CELLS = 30
SCEN_SEED = {"balanced": 11, "imbalanced": 22, "rare": 33, "similar_subtypes": 44}


def _within_family_de_panels(ref, mapping, n_top=40):
    """Oracle discriminative panel: top within-family DE genes from the REFERENCE
    (not the query) — for synthetic data generated from these profiles these ARE
    the truly discriminative genes. Respects rule 11 (no test-donor info)."""
    R = ref.as_R_cpm(); genes = [str(g) for g in ref.gene_names]
    types = [str(c) for c in ref.cell_types]
    fam_members: dict = {}
    for i, t in enumerate(types):
        fam_members.setdefault(mapping.get(t, t), []).append((i, t))
    panels = {}
    for fam, members in fam_members.items():
        if len(members) < 2:
            continue
        idx = [i for i, _ in members]
        sub = R[idx]                                   # members × genes (CPM)
        # within-family discriminability: variance across members / mean (CV-like)
        mu = sub.mean(0) + 1e-6
        score = sub.std(0) / mu                        # genes that separate subtypes
        top = np.argsort(score)[::-1][:n_top]
        panels[fam] = [genes[j] for j in top]
    return panels


def _fine_metrics(truth, pred, mapping):
    cols = [c for c in truth.columns if not str(c).startswith("unresolved_")]
    t = truth[cols]; p = pred.reindex(index=t.index, columns=cols).fillna(0.0)
    acc = M.accuracy_metrics(t, p)
    # within-family conditional accuracy: renormalise within each family, compare
    fam_of = {c: mapping.get(c, c) for c in cols}
    cond_t, cond_p = t.copy() * 0.0, p.copy() * 0.0
    for fam in set(fam_of.values()):
        mem = [c for c in cols if fam_of[c] == fam]
        st, sp = t[mem].sum(1), p[mem].sum(1)
        for c in mem:
            cond_t[c] = np.where(st > 0, t[c] / st.replace(0, np.nan), 0.0)
            cond_p[c] = np.where(sp > 0, p[c] / sp.replace(0, np.nan), 0.0)
    cond = M.accuracy_metrics(cond_t.fillna(0), cond_p.fillna(0))
    cx_t = M.composition_complexity(t)["effective_n_populations"].mean()
    cx_p = M.composition_complexity(p)["effective_n_populations"].mean()
    return {"fine_pearson": acc["pearson"], "fine_rmse": acc["rmse"],
            "within_family_cond_pearson": cond["pearson"],
            "true_eff_n": float(cx_t), "pred_eff_n": float(cx_p),
            "complexity_abs_error": float(abs(cx_p - cx_t)),
            "mass_sum_mean": float(p.sum(1).mean())}


def run_bulk(adata, raw_map, rows):
    from tissueresolve.api import deconv_bulk
    from tissueresolve.reference.hierarchy import (
        build_cell_type_hierarchy, combine_family_and_conditional_estimates)
    for s in BULK_SEEDS:
        ref_donors, query = SH.split_donors(adata, "donor_id", ref_frac=0.5, seed=s)
        rmask = adata.obs["donor_id"].astype(str).isin(set(ref_donors)).to_numpy()
        with warnings.catch_warnings():
            warnings.simplefilter("ignore")
            ref = H.prepare_reference(adata[rmask].copy(), min_cells=MIN_CELLS,
                                      estimate_overdispersion=True).reference
        ref_types = [str(c) for c in ref.cell_types]
        mapping = build_cell_type_hierarchy(ref_types, raw_map)
        oracle_panels = _within_family_de_panels(ref, mapping)
        global_panel = ref_types  # placeholder; global = full gene set (panel=None below)
        for scen in BULK_SCENARIOS:
            kw = {"similar_pair": tuple([t for t in ref_types
                  if mapping.get(t) == max({mapping.get(x, x) for x in ref_types},
                  key=lambda f: sum(mapping.get(x, x) == f for x in ref_types))][:2])} \
                  if scen == "similar_subtypes" else {}
            tgt = SH.build_target_proportions(ref_types, N_SAMPLES, scen,
                                              seed=SCEN_SEED[scen], **kw)
            ds = SH.realize_pseudobulk(adata, tgt, celltype_col="cell_type",
                                       donor_col="donor_id", query_donors=query,
                                       seed=s, cells_per_sample=CELLS_PER_SAMPLE)
            truth = ds.true_mrna_proportions
            true_fam = M.aggregate_to_families(truth, mapping)
            counts = ds.counts

            def emit(method, pred, extra=None):
                m = _fine_metrics(truth, pred, mapping)
                if extra:
                    m.update(extra)
                rows.append({"experiment": "bulk", "seed": s, "scenario": scen,
                             "method": method, **m})

            with warnings.catch_warnings():
                warnings.simplefilter("ignore")
                # flat baseline
                flat = deconv_bulk(counts, ref, solver="nnls", resolution_mode="none").deconv.proportions
                emit("flat_nnls", flat)
                # standard hierarchical (current defaults)
                h = deconv_bulk(counts, ref, solver="auto", resolution_mode="hierarchical",
                                hierarchy_mapping=mapping)
                est = h.estimates
                emit("hier_standard", est.combined_fine)
                # E3 pre-gating: all-subtype combine (no abstention)
                pre = combine_family_and_conditional_estimates(
                    est.family_proportions, est.conditional_proportions, mapping)
                emit("hier_pregating", pre)
                # E3 no-gating run (allow_unresolved=False)
                hng = deconv_bulk(counts, ref, solver="auto", resolution_mode="hierarchical",
                                  hierarchy_mapping=mapping, allow_unresolved=False,
                                  allow_partial_resolution=False)
                emit("hier_no_gating", hng.estimates.combined_fine)
                # E1 oracle broad mass: inject TRUE family mass into combine
                ora = combine_family_and_conditional_estimates(
                    true_fam.reindex(columns=est.family_proportions.columns).fillna(0.0),
                    est.conditional_proportions, mapping)
                emit("hier_oracle_broad", ora)
                # E2 oracle gene panel (within-family DE) + global (full genes)
                hop = deconv_bulk(counts, ref, solver="auto", resolution_mode="hierarchical",
                                  hierarchy_mapping=mapping, allow_unresolved=False,
                                  allow_partial_resolution=False, family_gene_panels=oracle_panels)
                emit("hier_oracle_panel_nogating", hop.estimates.combined_fine)
            # attach unresolved-aware info to the standard gated row
            ufams = set(getattr(est, "unresolved_families", []) or [])
            umass = float(est.unresolved_mass.to_numpy().sum() /
                          max(est.combined_fine.to_numpy().sum(), 1e-9)) \
                if getattr(est, "unresolved_mass", None) is not None else 0.0
            scope = [r for r in rows if r["seed"] == s and r["scenario"] == scen
                     and r["experiment"] == "bulk"]
            bym = {r["method"]: r for r in scope}
            bym["hier_standard"]["unresolved_mass_frac"] = umass
            bym["hier_standard"]["n_unresolved_families"] = len(ufams)
            print(f"  bulk seed{s}/{scen}: flat={bym['flat_nnls']['fine_pearson']:.3f} "
                  f"hier={bym['hier_standard']['fine_pearson']:.3f} "
                  f"pregate={bym['hier_pregating']['fine_pearson']:.3f} "
                  f"oracle_broad={bym['hier_oracle_broad']['fine_pearson']:.3f} "
                  f"oracle_panel={bym['hier_oracle_panel_nogating']['fine_pearson']:.3f}", flush=True)


def run_spatial(adata, raw_map, rows):
    import tissueresolve as tr
    from tissueresolve.config import TissueResolveConfig
    from tissueresolve.reference.hierarchy import build_cell_type_hierarchy
    lambdas = {"none_0.0": 0.0, "weak_0.02": 0.02, "default_0.1": 0.1, "strong_0.5": 0.5}
    for s in SPATIAL_SEEDS:
        ref_donors, query = SH.split_donors(adata, "donor_id", ref_frac=0.5, seed=s)
        rmask = adata.obs["donor_id"].astype(str).isin(set(ref_donors)).to_numpy()
        with warnings.catch_warnings():
            warnings.simplefilter("ignore")
            ref = H.prepare_reference(adata[rmask].copy(), min_cells=MIN_CELLS,
                                      estimate_overdispersion=True).reference
        ref_types = [str(c) for c in ref.cell_types]
        mapping = build_cell_type_hierarchy(ref_types, raw_map)
        rare = next((t for t in ["regulatory T cell", "natural killer cell"] if t in ref_types), ref_types[-1])
        sc = SSP.generate_spatial_scenario(adata, ref_types, query, mapping, n_side=14,
                                           cells_per_spot=40, seed=s, rare_type=rare)
        truth, coords, domains = sc["truth"], sc["coords"], sc["domain_labels"]
        for name, lam in lambdas.items():
            cfg = TissueResolveConfig()
            cfg.spatial_solver.lambda_spatial = lam
            with warnings.catch_warnings():
                warnings.simplefilter("ignore")
                try:
                    res = tr.deconv_spatial(sc["Y"], ref, sc["array_row"], sc["array_col"],
                                            sc["lib_sizes"], sc["gene_names"], spot_ids=sc["spot_ids"],
                                            resolution_mode="none", config=cfg, run_neighbourhood=False)
                    pred = res.deconv.proportions
                    lam_used = res.deconv.lambda_spatial
                except Exception as exc:  # noqa: BLE001
                    print(f"  spatial seed{s} lambda={name} FAILED: {exc}"); continue
            cols = list(truth.columns)
            p = pred.reindex(index=truth.index, columns=cols).fillna(0.0)
            acc = M.accuracy_metrics(truth, p); comp = M.compositional_metrics(truth, p)
            fid = SM.spatial_fidelity_metrics(truth, p, coords, domain_labels=domains, k=6)
            rows.append({"experiment": "spatial", "seed": s, "scenario": f"lambda_{name}",
                         "method": f"spatial_lambda_{name}", "lambda_used": float(lam_used),
                         "fine_pearson": acc["pearson"], "fine_rmse": acc["rmse"],
                         "fine_jsd": comp["jsd_mean"], "local_rmse": fid["local_rmse"],
                         "oversmoothing_score": fid["oversmoothing_score"],
                         "morans_i_mae": fid["morans_i_mae"],
                         "domain_recovery_ari": fid.get("domain_recovery_ari", float("nan"))})
            print(f"  spatial seed{s} {name} (lam={lam_used}): fine_r={acc['pearson']:.3f} "
                  f"oversmooth={fid['oversmoothing_score']:.2f} local_rmse={fid['local_rmse']:.3f}", flush=True)


def main(argv=None) -> int:
    ap = argparse.ArgumentParser(); ap.add_argument("--run-real-data", action="store_true")
    ap.add_argument("--skip-spatial", action="store_true")
    args = ap.parse_args(argv)
    if not H.real_data_enabled(args.run_real_data):
        print(f"Refusing without --run-real-data.", file=sys.stderr); return 2
    import anndata as ad
    from tissueresolve.reference.hierarchy import load_hierarchy_mapping
    adata = ad.read_h5ad(H.REFERENCE_H5AD)
    if "feature_name" in adata.var.columns:
        fn = adata.var["feature_name"].astype(str)
        if fn.nunique() == adata.n_vars:
            adata.var_names = fn.to_numpy()
    raw_map = load_hierarchy_mapping(HIER)
    rows = []
    print("E1-E3 bulk oracle experiments ..."); run_bulk(adata, raw_map, rows)
    if not args.skip_spatial:
        print("E4 spatial smoothing sweep ..."); run_spatial(adata, raw_map, rows)
    df = pd.DataFrame(rows)
    OUT.mkdir(parents=True, exist_ok=True)
    with open(OUT / "hierarchical_oracle_experiments.tsv", "w") as fh:
        fh.write("# PHASE 1 hierarchical-failure diagnosis (oracle experiments). DIAGNOSIS ONLY; no algorithm changed.\n")
        df.to_csv(fh, sep="\t", index=False)
    print(f"\nWrote {OUT/'hierarchical_oracle_experiments.tsv'}")
    # summaries
    b = df[df.experiment == "bulk"]
    if not b.empty:
        print("\n=== BULK: mean fine Pearson by method (over seeds×scenarios) ===")
        print(b.groupby("method")[["fine_pearson", "within_family_cond_pearson",
              "complexity_abs_error"]].mean().round(3).to_string())
    sp = df[df.experiment == "spatial"]
    if not sp.empty:
        print("\n=== SPATIAL: mean by lambda ===")
        print(sp.groupby("scenario")[["fine_pearson", "oversmoothing_score", "local_rmse",
              "domain_recovery_ari"]].mean().round(3).to_string())
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
