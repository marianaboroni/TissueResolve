#!/usr/bin/env python
"""Strategies C (union control) and D (hierarchical multi-panel) vs B (flat fine_global),
Stage 1.5 §8-10. Donor-disjoint, Poisson, real-data. Reuses the EXISTING hierarchical
inference (broad->fine->resolution-gate->unresolved, mass-conserving) — no new solver.

  B_fine_global_flat   : get_solver('poisson', genes=fine_global).solve  (flat)
  C_union_flat         : flat Poisson on broad+sibling+rare union (negative control)
  D_hierarchical       : deconv_bulk(solver='poisson', resolution_mode='hierarchical')
                         with ref.selected_genes = fine_global (broad+fine+gating+unresolved)

Reports conditional within-family RMSE, mass conservation, resolved coverage, paired-by-mixture.
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
sys.path.insert(0, str(REPO / "src"))
import _harness as H  # noqa: E402
from benchmarks.spatial.run_weak_smoothing_grid import _load_dataset  # noqa: E402

OUT = REPO / "benchmarks" / "results" / "signatures"


def _cond(truth, pred, mp):
    fam = {}
    for c in truth.columns:
        fam.setdefault(str(mp.get(c, c)), []).append(c)
    multi = [c for m in fam.values() if len(m) > 1 for c in m]
    if not multi:
        return np.nan

    def cd(df):
        o = pd.DataFrame(0.0, index=df.index, columns=truth.columns)
        for m in fam.values():
            m = [x for x in m if x in df.columns]
            if len(m) < 2:
                continue
            s = df[m].sum(1).replace(0, np.nan)
            for x in m:
                o[x] = (df[x] / s).fillna(0)
        return o
    tc, pc = cd(truth), cd(pred.reindex(columns=truth.columns).fillna(0))
    return float(np.sqrt(np.mean((tc[multi].to_numpy(float) - pc[multi].to_numpy(float)) ** 2)))


def _rmse(truth, pred):
    p = pred.reindex(index=truth.index, columns=truth.columns).fillna(0)
    return float(np.sqrt(np.mean((p.to_numpy(float) - truth.to_numpy(float)) ** 2)))


def run(dataset, seeds, rows):
    import tissueresolve as tr
    from tissueresolve.solver import get_solver
    from tissueresolve.reference.signature_optimizer import ReferenceSignatureOptimizer
    D = _load_dataset(dataset)
    ad, ctc, dc, ref, types = (D["adata"], D["celltype_col"], D["donor_col"], D["ref"], D["ref_types"])
    mp = {t: D["mapping"].get(t, t) for t in types}
    train = ad[ad.obs[dc].astype(str).isin(set(D["ref_donors"])).to_numpy()].copy()
    with warnings.catch_warnings():
        warnings.simplefilter("ignore")
        m = ReferenceSignatureOptimizer(ctc, dc, mapping=mp, seed=0).optimize(train)
    fine = m.fine_global_genes
    union = list(dict.fromkeys(list(m.broad_genes)
                 + [g for gl in m.sibling_genes.values() for g in gl]
                 + [g for gl in m.rare_confirmation.values() for g in gl]))
    import copy
    ref_fine = copy.deepcopy(ref); ref_fine.selected_genes = list(fine)
    qmask = ad.obs[dc].astype(str).isin(set(D["query_donors"])).to_numpy()
    for seed in seeds:
        counts, truth, _ = H.generate_pseudobulk(ad[qmask].copy(), ctc,
                                                 n_per_regime=2, n_cells=200, seed=seed)
        with warnings.catch_warnings():
            warnings.simplefilter("ignore")
            preds = {
                "B_fine_global_flat": get_solver("poisson", genes=fine).solve(counts, ref).proportions,
                "C_union_flat": get_solver("poisson", genes=union).solve(counts, ref).proportions,
                "D_hierarchical": tr.deconv_bulk(counts, ref_fine, solver="poisson",
                                                 resolution_mode="hierarchical",
                                                 hierarchy_mapping=mp, n_bootstrap=0).deconv.proportions,
            }
        for name, p in preds.items():
            # mass conservation: rows (incl unresolved_*) sum to 1
            mass = float(p.sum(axis=1).mean())
            fine_cols = [c for c in p.columns if not str(c).startswith("unresolved_")]
            unresolved = float(p[[c for c in p.columns if str(c).startswith("unresolved_")]].sum(axis=1).mean()) \
                if any(str(c).startswith("unresolved_") for c in p.columns) else 0.0
            rows.append(dict(dataset=dataset, seed=seed, strategy=name,
                             RMSE=_rmse(truth, p[fine_cols]), cond_RMSE=_cond(truth, p[fine_cols], mp),
                             mass_sum=mass, unresolved_mass=unresolved,
                             resolved_coverage=1.0 - unresolved))


def main(argv=None):
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument("--run-real-data", action="store_true")
    ap.add_argument("--datasets", nargs="+", default=["breast", "lung"])
    ap.add_argument("--seeds", type=int, nargs="+", default=[0, 1, 2])
    args = ap.parse_args(argv)
    if not H.real_data_enabled(args.run_real_data):
        print(f"Refusing without --run-real-data (or {H.REAL_DATA_ENV}=1).", file=sys.stderr)
        return 2
    OUT.mkdir(parents=True, exist_ok=True)
    rows = []
    for ds in args.datasets:
        run(ds, args.seeds, rows)
    df = pd.DataFrame(rows)
    df.to_csv(OUT / "strategies_cd_rows.tsv", sep="\t", index=False)
    summ = df.groupby(["dataset", "strategy"])[["RMSE", "cond_RMSE", "mass_sum",
                                                "unresolved_mass", "resolved_coverage"]].mean().round(4)
    summ.to_csv(OUT / "strategies_cd_summary.tsv", sep="\t")
    (OUT / "strategies_cd_manifest.json").write_text(json.dumps(
        {"datasets": args.datasets, "seeds": args.seeds, "solver": "poisson",
         "D_reuses": "deconv_bulk hierarchical (broad->fine->gate->unresolved)"}, indent=2))
    print("\n=== STRATEGIES B/C/D (Poisson, donor-disjoint) ===")
    print(summ.to_string())
    print(f"\nWrote {OUT}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
