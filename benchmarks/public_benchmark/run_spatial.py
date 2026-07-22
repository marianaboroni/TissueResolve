#!/usr/bin/env python
"""Public benchmark — SPATIAL (synthetic, known ground truth). Same reference, same synthetic spots,
same scoring for every method. Runs TissueResolve default / weak_smoothing / edge_aware presets +
NNLS-per-spot baseline, scores broad + fine + local + boundary + rare-niche + absent-mass (via the
existing spatial scorer) plus an identifiability-stratified layer at spot depth. External spatial
tools (RCTD/CARD/cell2location) are attempted separately (run_external_spatial.*) and recorded
honestly. No default changed; TR presets are opt-in config only.

Usage: PYTHONPATH=src:. python benchmarks/public_benchmark/run_spatial.py --run-real-data \
       --datasets breast lung --seeds 0 1
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
sys.path.insert(0, str(REPO)); sys.path.insert(0, str(REPO / "src"))

from benchmarks.spatial.run_weak_smoothing_grid import _load_dataset, _score, DATASET_CFG  # noqa: E402
from benchmarks.shared import synthetic_spatial as SSP  # noqa: E402

OUT = REPO / "benchmarks" / "results" / "public_benchmark" / "spatial"
PRESETS = {"TissueResolve_default": dict(lambda_spatial=0.1, edge_aware=False),
           "TissueResolve_weak_smoothing": dict(lambda_spatial=0.02, edge_aware=False),
           "TissueResolve_edge_aware": dict(lambda_spatial=0.05, edge_aware=True)}


def _spatial_identifiability(ref, Y, gene_names, types):
    """Certificate at spot depth: pseudo-query = summed spots, library_size = median spot depth."""
    import tissueresolve as tr
    dfq = pd.DataFrame(np.asarray(Y).T, index=[str(g) for g in gene_names],
                       columns=[f"s{i}" for i in range(np.asarray(Y).shape[0])])
    with warnings.catch_warnings():
        warnings.simplefilter("ignore")
        cert = tr.bulk_identifiability(dfq, ref)
    return dict(zip(cert.per_type.cell_type, cert.per_type.recoverability))


def _ident_stratified(truth, pred, class_by):
    cols = [c for c in truth.columns if not str(c).startswith("unresolved_")]
    t = truth[cols]; p = pred.reindex(index=t.index, columns=cols).fillna(0.0)
    from collections import defaultdict
    by = defaultdict(list)
    for c in cols:
        by[class_by.get(c, "NOT_TESTABLE")].append(
            float(np.sqrt(np.mean((p[c].to_numpy(float) - t[c].to_numpy(float)) ** 2))))
    return {f"err_{k}": (float(np.mean(v)) if v else np.nan)
            for k, v in [(k, by[k]) for k in
            ("RESOLVABLE", "WEAKLY_RESOLVABLE", "GROUP_ONLY", "UNRESOLVABLE", "NOT_TESTABLE")]}


def run(dataset, seeds, rows, status_rows):
    import tissueresolve as tr
    from tissueresolve.config import TissueResolveConfig
    from benchmarks.spatial.methods.nnls_spot_baseline import NNLSSpotBaseline

    D = _load_dataset(dataset)
    types = [str(c) for c in D["ref_types"]]
    for seed in seeds:
        sc = SSP.generate_spatial_scenario(
            D["adata"], D["ref_types"], D["query_donors"], D["mapping"],
            celltype_col=D["celltype_col"], donor_col=D["donor_col"],
            n_side=12, cells_per_spot=8, seed=seed, rare_type=D["rare_type"],
            domain_families=D["domain_families"])
        truth, coords, domains = sc["truth"], sc["coords"], sc["domain_labels"]
        class_by = _spatial_identifiability(D["ref"], sc["Y"], sc["gene_names"], types)
        # TissueResolve presets
        for name, cfgkw in PRESETS.items():
            tcfg = TissueResolveConfig()
            tcfg.spatial_solver.lambda_spatial = cfgkw["lambda_spatial"]
            tcfg.spatial_solver.edge_aware = cfgkw["edge_aware"]
            t0 = time.perf_counter()
            try:
                with warnings.catch_warnings():
                    warnings.simplefilter("ignore")
                    res = tr.deconv_spatial(sc["Y"], D["ref"], sc["array_row"], sc["array_col"],
                                            sc["lib_sizes"], sc["gene_names"], spot_ids=sc["spot_ids"],
                                            config=tcfg, resolution_mode="none", run_neighbourhood=False)
                pred = res.deconv.proportions
                rt = time.perf_counter() - t0
                m = _score(truth, pred, coords, domains, D["mapping"], D["rare_type"])
                rows.append(dict(dataset=dataset, seed=seed, method=name, runtime_s=round(rt, 2),
                                 **m, **_ident_stratified(truth, pred, class_by)))
                status_rows.append(dict(dataset=dataset, method=name, status="executed",
                                        runtime_s=round(rt, 2), note=""))
            except Exception as e:
                status_rows.append(dict(dataset=dataset, method=name, status="failed_run",
                                        runtime_s=round(time.perf_counter() - t0, 2),
                                        note=f"{type(e).__name__}: {str(e)[:150]}"))
        # NNLS-per-spot baseline
        t0 = time.perf_counter()
        try:
            with warnings.catch_warnings():
                warnings.simplefilter("ignore")
                r = NNLSSpotBaseline().run({"reference": D["ref"], "Y": sc["Y"],
                                            "gene_names": sc["gene_names"], "normalization_status": "counts"})
            if r.status != "success":
                raise RuntimeError(r.skip_reason or r.status)
            pred = r.predictions
            rt = time.perf_counter() - t0
            m = _score(truth, pred, coords, domains, D["mapping"], D["rare_type"])
            rows.append(dict(dataset=dataset, seed=seed, method="NNLS_spot_baseline",
                             runtime_s=round(rt, 2), **m, **_ident_stratified(truth, pred, class_by)))
            status_rows.append(dict(dataset=dataset, method="NNLS_spot_baseline", status="executed",
                                    runtime_s=round(rt, 2), note=""))
        except Exception as e:
            status_rows.append(dict(dataset=dataset, method="NNLS_spot_baseline", status="failed_run",
                                    runtime_s=round(time.perf_counter() - t0, 2),
                                    note=f"{type(e).__name__}: {str(e)[:150]}"))
        print(f"[{dataset} seed={seed}] spots={len(truth)} types={len(types)} "
              f"classes={pd.Series(class_by).value_counts().to_dict()}")


def main(argv=None):
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument("--run-real-data", action="store_true")
    ap.add_argument("--datasets", nargs="+", default=["breast", "lung"])
    ap.add_argument("--seeds", type=int, nargs="+", default=[0, 1])
    args = ap.parse_args(argv)
    import _harness as H
    if not H.real_data_enabled(args.run_real_data):
        print(f"Refusing without --run-real-data (or {H.REAL_DATA_ENV}=1).", file=sys.stderr); return 2
    OUT.mkdir(parents=True, exist_ok=True)
    rows, srows = [], []
    t0 = time.perf_counter()
    for ds in args.datasets:
        if not Path(DATASET_CFG.get(ds, {}).get("h5ad", "x")).exists():
            print(f"[skip] {ds}: missing h5ad"); continue
        run(ds, args.seeds, rows, srows)
    pd.DataFrame(rows).to_csv(OUT / "metrics_by_scenario.tsv", sep="\t", index=False)
    pd.DataFrame(srows).to_csv(OUT / "method_status.tsv", sep="\t", index=False)
    (OUT / "run_manifest.json").write_text(json.dumps(
        {"datasets": args.datasets, "seeds": args.seeds, "presets": list(PRESETS),
         "runtime_s": round(time.perf_counter() - t0, 1)}, indent=2, default=str))
    m = pd.DataFrame(rows)
    if len(m):
        cols = [c for c in ("broad_rmse", "broad_pearson", "fine_rmse", "fine_pearson", "cond_rmse",
                            "rare_niche_recall", "absent_mass", "runtime_s") if c in m.columns]
        print("\n=== SPATIAL SUMMARY (mean) ===\n", m.groupby("method")[cols].mean().round(4).to_string())
    print(f"\nWrote {OUT}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
