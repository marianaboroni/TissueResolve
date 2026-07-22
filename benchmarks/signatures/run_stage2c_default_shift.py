#!/usr/bin/env python
"""Stage 2C — validate a FIXED default shift_scale (no per-run calibration split).

Stage 2B showed a calibration-fit donor-shift inflation (~2) gives ~nominal coverage. For the
certificate to be usable as a default core output *without* requiring the user to hold out a
calibration split, a single fixed default must give acceptable held-out coverage across tissues.

This checks predicted-interval coverage (nominal 90%) at fixed shift_scale in {1,2,3} on breast+lung
(donor-disjoint reference vs held-out), full-depth and 1000-gene panel. No fitting; no default changed.

Usage: PYTHONPATH=src:. python benchmarks/signatures/run_stage2c_default_shift.py --run-real-data
"""
from __future__ import annotations

import argparse
import json
import sys
import warnings
from collections import defaultdict
from pathlib import Path

import numpy as np
import pandas as pd

REPO = Path(__file__).resolve().parents[2]
sys.path.insert(0, str(REPO / "examples" / "real_breast_cancer" / "scripts"))
sys.path.insert(0, str(REPO))
sys.path.insert(0, str(REPO / "src"))

from tissueresolve.reference.identifiability_calibration import (  # noqa: E402
    calibrated_identifiability_certificate as ccert)

OUT = REPO / "benchmarks" / "results" / "signatures" / "stage2c"
NOMINAL, Z, N_SIM = 0.90, 1.645, 150
CONDITIONS = [("full_depth", 1.0, None), ("panel_1000", 1.0, 1000)]


def _pool(query, ctc, seeds, H):
    cl, tl = [], []
    for s in seeds:
        c, t, _ = H.generate_pseudobulk(query, ctc, n_per_regime=3, n_cells=200, seed=s)
        cl.append(c); tl.append(t)
    counts = pd.concat(cl, axis=1); counts.columns = [f"m{i}" for i in range(counts.shape[1])]
    truth = pd.concat(tl, axis=0).reset_index(drop=True); truth.index = counts.columns
    return counts, truth


def _thin(counts, f, rng):
    return counts if f >= 1 else pd.DataFrame(rng.binomial(counts.to_numpy().astype(np.int64), f),
                                              index=counts.index, columns=counts.columns)


def _coverage(truth, pred, pred_sd, types):
    ins = []
    for k in types:
        tv = truth[k].to_numpy(float) if k in truth else np.zeros(len(truth))
        pv = pred[k].to_numpy(float) if k in pred else np.zeros(len(truth))
        sd = float(pred_sd.get(k, np.nan))
        if sd == sd:
            ins.extend((np.abs(pv - tv) <= Z * sd + 1e-9).tolist())
    return round(float(np.mean(ins)), 4) if ins else np.nan


def run(dataset, seeds, rows):
    import _harness as H
    from benchmarks.spatial.run_weak_smoothing_grid import _load_dataset
    from tissueresolve.reference.build import ReferenceBuilder
    from tissueresolve.config import ReferenceConfig
    from tissueresolve.solver import PoissonGLMSolver

    D = _load_dataset(dataset)
    ad, ctc, dc = D["adata"], D["celltype_col"], D["donor_col"]
    train = set(map(str, D["ref_donors"])); test = set(map(str, D["query_donors"]))
    train_ad = ad[ad.obs[dc].astype(str).isin(train).to_numpy()].copy()
    with warnings.catch_warnings():
        warnings.simplefilter("ignore")
        ref = ReferenceBuilder(ReferenceConfig(celltype_col=ctc, donor_col=dc, min_cells=25)).build_from_adata(train_ad)
    types = [str(c) for c in ref.cell_types]; genes = list(ref.gene_names)
    q_test = ad[ad.obs[dc].astype(str).isin(test).to_numpy()].copy()
    counts_full, truth = _pool(q_test, ctc, seeds, H)
    truth = truth.reindex(columns=types).fillna(0.0)
    panel_1000 = list(np.random.default_rng(12345).choice(genes, size=min(1000, len(genes)), replace=False))

    for cname, depth, n_panel in CONDITIONS:
        rng = np.random.default_rng(abs(hash((dataset, cname))) % (2**32))
        panel = panel_1000 if n_panel else genes
        counts = _thin(counts_full, depth, rng)
        det = [g for g in panel if (counts.loc[panel] > 0).mean(axis=1).get(g, 0) >= 0.10]
        lib = float(counts.loc[det].sum(axis=0).median())
        with warnings.catch_warnings():
            warnings.simplefilter("ignore")
            pred = PoissonGLMSolver(genes=(panel if n_panel else None)).solve(
                counts, ref).proportions.reindex(columns=types).fillna(0.0)
        for s in (1.0, 2.0, 3.0):
            cert = ccert(ref, query_detectable_genes=det, library_size=lib, n_sim=N_SIM, shift_scale=s, seed=0)
            sd = dict(zip(cert.per_type.cell_type, cert.per_type.pred_sd))
            rows.append(dict(dataset=dataset, condition=cname, shift_scale=s,
                             coverage=_coverage(truth, pred, sd, types)))
            print(f"  {dataset:7s} {cname:11s} s={s}: coverage={rows[-1]['coverage']}")


def main(argv=None):
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument("--run-real-data", action="store_true")
    ap.add_argument("--datasets", nargs="+", default=["breast", "lung"])
    ap.add_argument("--seeds", type=int, nargs="+", default=[0, 1, 2, 3, 4])
    args = ap.parse_args(argv)
    import _harness as H
    if not H.real_data_enabled(args.run_real_data):
        print(f"Refusing without --run-real-data (or {H.REAL_DATA_ENV}=1).", file=sys.stderr); return 2
    OUT.mkdir(parents=True, exist_ok=True)
    rows = []
    for ds in args.datasets:
        run(ds, args.seeds, rows)
    df = pd.DataFrame(rows); df.to_csv(OUT / "default_shift_coverage.tsv", sep="\t", index=False)
    piv = df.pivot_table(index="shift_scale", values="coverage", aggfunc="mean").round(4)
    best = (df.groupby("shift_scale").coverage.apply(lambda v: float(np.mean(np.abs(v - NOMINAL))))).idxmin()
    print("\n=== mean coverage by fixed shift_scale (target 0.90) ===\n", piv.to_string())
    print(f"closest-to-nominal fixed default: shift_scale={best}")
    (OUT / "default_shift_manifest.json").write_text(json.dumps(
        {"datasets": args.datasets, "seeds": args.seeds, "nominal": NOMINAL,
         "recommended_default_shift_scale": float(best),
         "mean_coverage_by_s": piv.coverage.to_dict()}, indent=2))
    print(f"Wrote {OUT}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
