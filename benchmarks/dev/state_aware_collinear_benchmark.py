#!/usr/bin/env python
"""Test whether state-aware 3-level deconvolution improves COLLINEAR fine-type
recovery vs the standard 2-level hierarchical path (breast pseudobulk gold-truth).

Setup (collinear states become "states" within a well-separated "cell type"):
  state      = fine cell type (the collinear subtypes, e.g. the 11 T/NK subtypes)
  cell type  = broad family (T/NK, Myeloid, Endothelial, ...)   <- where collinear states live
  broad      = super-compartment (Immune / Stromal / Epithelial)

Compares, on the SAME donor-disjoint pseudobulk and SAME fine truth:
  - baseline : 2-level hierarchical  (family -> fine)            [deconv_bulk hierarchical]
  - state_aware : 3-level            (compartment -> family -> state)

Primary question: does the extra level improve conditional WITHIN-FAMILY RMSE for
collinear families (T/NK, Myeloid, Endothelial) without hurting fine accuracy,
rare recall, or spillover? Honest gates; experimental; no defaults changed.

Usage:
  PYTHONPATH=src:. python benchmarks/dev/state_aware_collinear_benchmark.py \
      --run-real-data --seeds 0 1 2 3 4
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
sys.path.insert(0, str(REPO / "src"))
import _harness as H  # noqa: E402
from benchmarks.shared import synthetic_holdout as SH  # noqa: E402
from benchmarks.shared import metrics as M  # noqa: E402
from benchmarks.spatial.run_weak_smoothing_grid import _load_dataset  # noqa: E402

OUT = REPO / "benchmarks" / "outputs" / "state_aware_collinear"
SCENARIOS = ["similar_subtypes", "imbalanced", "rare"]
N_SAMPLES = 12
CELLS_PER_SAMPLE = 600
# super-compartment grouping of breast broad families
FAMILY_TO_COMPARTMENT = {
    "T/NK": "Immune", "Myeloid": "Immune", "B/Plasma": "Immune",
    "Endothelial": "Stromal", "Mural": "Stromal", "Stromal/Fibroblast": "Stromal",
    "Adipocyte": "Stromal", "Epithelial": "Epithelial",
}


def _conditional_rmse(truth, pred, mapping, only_family=None):
    """RMSE of conditional P(state|family) over families with >=2 fine members."""
    fam_members = {}
    for c in truth.columns:
        fam = str(mapping.get(c, c))
        if only_family is not None and fam != only_family:
            continue
        fam_members.setdefault(fam, []).append(c)
    sq = []
    for fam, members in fam_members.items():
        members = [m for m in members if m in pred.columns]
        if len(members) < 2:
            continue
        ts = truth[members].sum(axis=1).replace(0, np.nan)
        ps = pred[members].sum(axis=1).replace(0, np.nan)
        for m in members:
            tc = (truth[m] / ts).fillna(0.0).to_numpy()
            pc = (pred[m] / ps).fillna(0.0).to_numpy()
            sq.append((tc - pc) ** 2)
    return float(np.sqrt(np.mean(np.concatenate(sq)))) if sq else float("nan")


def _score(truth, pred, mapping, rare, collinear_families):
    cols = [c for c in truth.columns if not str(c).startswith("unresolved_")]
    t = truth[cols]
    p = pred.reindex(index=t.index, columns=cols).fillna(0.0)
    acc = M.accuracy_metrics(t, p)
    tfam, pfam = M.aggregate_to_families(t, mapping), M.aggregate_to_families(p, mapping)
    bacc = M.accuracy_metrics(tfam, pfam)
    absent = t.to_numpy(float) < 1e-9
    spill = float(np.where(absent, p.to_numpy(float), 0.0).sum(axis=1).mean())
    fp_sub = float(((p.to_numpy(float) > 0.01) & absent).mean())
    effn_p = float(np.mean([M.effective_n_populations(p.iloc[i].to_numpy()) for i in range(p.shape[0])]))
    effn_t = float(np.mean([M.effective_n_populations(t.iloc[i].to_numpy()) for i in range(t.shape[0])]))
    # rare recall/precision
    if rare in t.columns and rare in p.columns:
        tp = t[rare].to_numpy() > 0.02; pp = p[rare].to_numpy() > 0.02
        rsens = float((tp & pp).sum() / max((tp).sum(), 1))
        rprec = float((tp & pp).sum() / max((pp).sum(), 1)) if pp.sum() else float("nan")
    else:
        rsens = rprec = float("nan")
    # unresolved mass actually returned by the method (mass not in cols)
    unres = float(pred.reindex(index=t.index).filter(regex="^unresolved_").sum(axis=1).mean()) \
        if any(str(c).startswith("unresolved_") for c in pred.columns) else 0.0
    out = {"fine_pearson": acc["pearson"], "fine_rmse": acc["rmse"],
           "broad_pearson": bacc["pearson"], "conditional_rmse": _conditional_rmse(t, p, mapping),
           "rare_sensitivity": rsens, "rare_precision": rprec,
           "pairwise_spillover": spill, "false_positive_subtype_rate": fp_sub,
           "effective_n_pred": effn_p, "effective_n_truth": effn_t, "unresolved_mass": unres}
    for fam in collinear_families:
        out[f"cond_rmse[{fam}]"] = _conditional_rmse(t, p, mapping, only_family=fam)
    return out


def main(argv=None) -> int:
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument("--run-real-data", action="store_true")
    ap.add_argument("--seeds", type=int, nargs="+", default=[0, 1, 2, 3, 4])
    args = ap.parse_args(argv)
    if not H.real_data_enabled(args.run_real_data):
        print(f"Refusing without --run-real-data (or {H.REAL_DATA_ENV}=1).", file=sys.stderr)
        return 2

    import tissueresolve as tr
    from tissueresolve.config import TissueResolveConfig
    from tissueresolve.reference.three_level_hierarchy import build_three_level_hierarchy
    from tissueresolve.bulk.state_aware_hierarchical import run_state_aware_hierarchical_bulk

    D = _load_dataset("breast")
    adata, ref, mapping = D["adata"], D["ref"], D["mapping"]   # mapping: fine -> family
    ct_col, donor_col = D["celltype_col"], D["donor_col"]
    ref_types = [str(c) for c in ref.cell_types]
    rare = D["rare_type"]
    q_donors = D["query_donors"]

    # 3-level maps restricted to ref types present
    state_to_celltype = {t: str(mapping.get(t, t)) for t in ref_types}
    celltype_to_broad = {fam: FAMILY_TO_COMPARTMENT.get(fam, "Other")
                         for fam in set(state_to_celltype.values())}
    hierarchy = build_three_level_hierarchy(state_to_celltype, celltype_to_broad)
    fam_members = {}
    for t, f in state_to_celltype.items():
        fam_members.setdefault(f, []).append(t)
    collinear_families = sorted([f for f, m in fam_members.items() if len(m) >= 3])
    sim_pair = tuple(sorted(fam_members.get("T/NK", ref_types))[:2])
    print(f"[breast] {len(ref_types)} states; families={ {k:len(v) for k,v in fam_members.items()} }")
    print(f"collinear families (>=3 states): {collinear_families}; similar_pair={sim_pair}; rare={rare}")

    OUT.mkdir(parents=True, exist_ok=True)
    rows, rt_rows = [], []
    cfg = TissueResolveConfig()
    for seed in args.seeds:
        for scen in SCENARIOS:
            kw = {}
            if scen == "rare":
                kw = {"rare_type": rare, "rare_level": 0.01}
            elif scen == "similar_subtypes":
                kw = {"similar_pair": sim_pair}
            targets = SH.build_target_proportions(ref_types, N_SAMPLES, scen, seed=seed, **kw)
            ds = SH.realize_pseudobulk(adata, targets, celltype_col=ct_col, donor_col=donor_col,
                                       query_donors=q_donors, seed=seed,
                                       cells_per_sample=CELLS_PER_SAMPLE)
            counts, truth = ds.counts, ds.true_mrna_proportions

            # --- baseline: 2-level hierarchical (family -> fine) ---
            for method in ("hierarchical_2level", "state_aware_3level"):
                t0 = time.perf_counter()
                try:
                    with warnings.catch_warnings():
                        warnings.simplefilter("ignore")
                        if method == "hierarchical_2level":
                            res = tr.deconv_bulk(counts, ref, config=cfg, solver="auto",
                                                 resolution_mode="hierarchical",
                                                 hierarchy_mapping=mapping)
                            pred = res.estimates.combined_fine
                        else:
                            sa = run_state_aware_hierarchical_bulk(counts, ref, hierarchy, config=cfg)
                            pred = sa.state_proportions
                except Exception as exc:  # noqa: BLE001
                    rt_rows.append({"seed": seed, "scenario": scen, "method": method,
                                    "status": f"failed: {exc}"})
                    print(f"  seed={seed} {scen} {method} FAILED: {exc}")
                    continue
                rt = float(time.perf_counter() - t0)
                m = _score(truth, pred, mapping, rare, collinear_families)
                rows.append({"seed": seed, "scenario": scen, "method": method, **m})
                rt_rows.append({"seed": seed, "scenario": scen, "method": method,
                                "runtime_seconds": round(rt, 2), "status": "ok"})
                print(f"  seed={seed} {scen:<16} {method:<20} {rt:5.1f}s fine={m['fine_pearson']:.3f} "
                      f"cond={m['conditional_rmse']:.3f} cond[T/NK]={m.get('cond_rmse[T/NK]', float('nan')):.3f} "
                      f"unres={m['unresolved_mass']:.3f}")

    df = pd.DataFrame(rows)
    df.to_csv(OUT / "per_scenario_metrics.tsv", sep="\t", index=False)
    pd.DataFrame(rt_rows).to_csv(OUT / "runtime_metrics.tsv", sep="\t", index=False)
    mcols = [c for c in df.columns if c not in ("seed", "scenario", "method")]
    by_method = df.groupby("method")[mcols].mean(numeric_only=True).reset_index()
    by_method.to_csv(OUT / "by_method_metrics.tsv", sep="\t", index=False)

    # gates: state_aware vs hierarchical_2level (means)
    gates = {}
    if {"hierarchical_2level", "state_aware_3level"} <= set(by_method["method"].values):
        b = by_method.set_index("method").loc["hierarchical_2level"]
        v = by_method.set_index("method").loc["state_aware_3level"]
        gates = {
            "conditional_rmse_improved": bool(v.conditional_rmse < b.conditional_rmse),
            "cond_rmse_T_NK_improved": bool(v.get("cond_rmse[T/NK]", np.nan) < b.get("cond_rmse[T/NK]", np.nan)),
            "fine_pearson_not_worse_1pct": bool(v.fine_pearson >= b.fine_pearson - 0.01),
            "rare_sens_not_worse": bool(v.rare_sensitivity >= b.rare_sensitivity - 1e-9),
            "spillover_not_increased": bool(v.pairwise_spillover <= b.pairwise_spillover + 1e-9),
            "effn_closer_or_equal": bool(abs(v.effective_n_pred - v.effective_n_truth)
                                         <= abs(b.effective_n_pred - b.effective_n_truth) + 1e-9),
        }
    (OUT / "gates.json").write_text(json.dumps(gates, indent=2), encoding="utf-8")
    (OUT / "manifest.json").write_text(json.dumps({
        "dataset": "breast", "seeds": args.seeds, "scenarios": SCENARIOS,
        "levels": "state=fine, cell_type=family, broad=compartment",
        "collinear_families": collinear_families, "rare_type": rare,
        "note": "experimental state-aware (opt-in); default behaviour unchanged",
    }, indent=2), encoding="utf-8")
    print(f"\nWrote -> {OUT}/")
    show = ["method", "fine_pearson", "broad_pearson", "conditional_rmse"] + \
           [f"cond_rmse[{f}]" for f in collinear_families] + \
           ["rare_sensitivity", "pairwise_spillover", "unresolved_mass"]
    show = [c for c in show if c in by_method.columns]
    print(by_method[show].round(3).to_string(index=False))
    print("\nGATES (state_aware vs 2-level):", json.dumps(gates))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
