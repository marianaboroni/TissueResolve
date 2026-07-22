#!/usr/bin/env python
"""Public benchmark — BULK. Same donor-disjoint reference, same held-out test mixtures, same ground
truth and scoring for every method. Runs TissueResolve variants + internal/external NNLS in-process,
scores them with the two-resolution + identifiability suite, and EXPORTS single-cell references +
mixtures so external R tools (MuSiC/Bisque/BayesPrism) can be run on identical inputs (run_external.sh)
and scored by the same code (merge_and_report.py).

Donor-disjoint: reference on TRAIN donors (built WITH donor_col so the identifiability certificate has
donor_cv); held-out mixtures from TEST donors. No test donor is used for gene/threshold/marker/model
selection. No TissueResolve default is changed.

Usage: PYTHONPATH=src:. python benchmarks/public_benchmark/run_bulk.py --run-real-data \
       --datasets breast lung --seeds 0 1 2
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

from benchmarks.public_benchmark import scoring as SC  # noqa: E402

OUT = REPO / "benchmarks" / "results" / "public_benchmark" / "bulk"
EXPORT = OUT / "external_inputs"
MAX_CELLS_PER_TYPE = 80
SCENARIOS = [("base", 1.0, None), ("low_depth", 0.03, None), ("reduced_overlap", 1.0, 1000)]


def _props(res):
    for attr in ("deconv",):
        d = getattr(res, attr, None)
        if d is not None and hasattr(d, "proportions"):
            return d.proportions
    return res.proportions


def _thin(counts, f, rng):
    return counts if f >= 1 else pd.DataFrame(
        rng.binomial(counts.to_numpy().astype(np.int64), f), index=counts.index, columns=counts.columns)


def _inprocess_methods():
    """name -> callable(bulk, ref, mapping) -> (proportions DataFrame). All share the same inputs."""
    import tissueresolve as tr
    from benchmarks.bulk.methods.nnls_baseline import NNLSBaseline
    from benchmarks.bulk.methods.extra_baselines import RidgeNNLSBaseline, CorrelationMatcherBaseline

    def _tr(bulk, ref, mapping, **kw):
        with warnings.catch_warnings():
            warnings.simplefilter("ignore")
            return _props(tr.deconv_bulk(bulk, ref, n_bootstrap=0, **kw))

    def _reg(cls, bulk, ref):
        with warnings.catch_warnings():
            warnings.simplefilter("ignore")
            r = cls().run({"bulk": bulk, "reference": ref, "normalization_status": "counts"})
        if r.status != "success":
            raise RuntimeError(r.skip_reason or r.status)
        return r.predictions

    return {
        "TissueResolve_wNNLS": lambda b, r, m: _tr(b, r, m, resolution_mode="flat"),
        "TissueResolve_PoissonGLM": lambda b, r, m: _tr(b, r, m, resolution_mode="flat", solver="poisson"),
        "TissueResolve_hierarchical": lambda b, r, m: _tr(b, r, m, resolution_mode="hierarchical",
                                                          hierarchy_mapping=m),
        "NNLS_baseline": lambda b, r, m: _reg(NNLSBaseline, b, r),
        "RidgeNNLS_baseline": lambda b, r, m: _reg(RidgeNNLSBaseline, b, r),
        "CorrelationMatcher_baseline": lambda b, r, m: _reg(CorrelationMatcherBaseline, b, r),
    }


def _export_reference(train_ad, ctc, dc, dest, rng):
    ct = train_ad.obs[ctc].astype(str).to_numpy()
    donor = train_ad.obs[dc].astype(str).to_numpy()
    keep = []
    for t in np.unique(ct):
        pool = np.where(ct == t)[0]
        keep.extend((rng.choice(pool, MAX_CELLS_PER_TYPE, replace=False)
                     if pool.size > MAX_CELLS_PER_TYPE else pool).tolist())
    keep = np.sort(np.array(keep))
    X = train_ad.X[keep]
    X = np.asarray(X.todense()) if hasattr(X, "todense") else np.asarray(X)
    cid = [f"c{i}" for i in keep]
    dest.mkdir(parents=True, exist_ok=True)
    pd.DataFrame(X.T.astype(int), index=list(map(str, train_ad.var_names)), columns=cid).to_csv(
        dest / "reference_counts_genes_by_cells.tsv", sep="\t")
    pd.DataFrame({"cell_id": cid, "cellType": ct[keep], "SubjectName": donor[keep]}).to_csv(
        dest / "reference_cell_metadata.tsv", sep="\t", index=False)


def run(dataset, seeds, metrics_rows, ct_rows, status_rows, cert_rows):
    import tissueresolve as tr
    import _harness as H
    from benchmarks.spatial.run_weak_smoothing_grid import _load_dataset
    from benchmarks.signatures import validation as V
    from tissueresolve.reference.build import ReferenceBuilder
    from tissueresolve.config import ReferenceConfig

    D = _load_dataset(dataset)
    ad, ctc, dc, mapping = D["adata"], D["celltype_col"], D["donor_col"], D["mapping"]
    train = set(map(str, D["ref_donors"])); test = set(map(str, D["query_donors"]))
    V.validate_donor_disjoint(train, test)
    train_ad = ad[ad.obs[dc].astype(str).isin(train).to_numpy()].copy()
    with warnings.catch_warnings():
        warnings.simplefilter("ignore")
        ref = ReferenceBuilder(ReferenceConfig(celltype_col=ctc, donor_col=dc, min_cells=25)).build_from_adata(train_ad)
    types = [str(c) for c in ref.cell_types]
    mp = {t: str(mapping.get(t, t)) for t in types}
    q_test = ad[ad.obs[dc].astype(str).isin(test).to_numpy()].copy()

    # held-out mixtures (easy/medium/hard regimes → balanced..imbalanced), known truth
    cl, tl = [], []
    for s in seeds:
        c, t, _ = H.generate_pseudobulk(q_test, ctc, n_per_regime=3, n_cells=200, seed=s)
        cl.append(c); tl.append(t)
    counts_full = pd.concat(cl, axis=1); counts_full.columns = [f"m{i}" for i in range(counts_full.shape[1])]
    truth = pd.concat(tl, axis=0).reset_index(drop=True); truth.index = counts_full.columns
    truth = truth.reindex(columns=types).fillna(0.0)

    genes = list(ref.gene_names)
    panel = list(np.random.default_rng(12345).choice(genes, size=min(1000, len(genes)), replace=False))
    _export_reference(train_ad, ctc, dc, EXPORT / dataset / "reference", np.random.default_rng(1))
    methods = _inprocess_methods()

    for scen, depth, n_panel in SCENARIOS:
        rng = np.random.default_rng(abs(hash((dataset, scen))) % (2**32))
        counts = _thin(counts_full, depth, rng)
        gpanel = panel if n_panel else None
        cbulk = counts.loc[panel] if n_panel else counts
        # identifiability certificate (shared evaluation axis) — donor-aware, default shift
        with warnings.catch_warnings():
            warnings.simplefilter("ignore")
            cert = tr.bulk_identifiability(cbulk, ref)
        class_by = dict(zip(cert.per_type.cell_type, cert.per_type.recoverability))
        clusters = [g["group"].split(";") for _, g in cert.clusters.iterrows()] if len(cert.clusters) else []
        for _, r in cert.per_type.iterrows():
            cert_rows.append(dict(dataset=dataset, scenario=scen, **r.to_dict()))
        # export mixtures for external tools
        ed = EXPORT / dataset / scen; ed.mkdir(parents=True, exist_ok=True)
        cbulk.to_csv(ed / "bulk_counts_genes_by_samples.tsv", sep="\t")
        truth.to_csv(ed / "truth_mrna.tsv", sep="\t")
        # run in-process methods on identical inputs
        for name, fn in methods.items():
            t0 = time.perf_counter()
            try:
                pred_raw = fn(cbulk, ref, mp)             # keep unresolved_* cols for fair broad scoring
                pred = pred_raw.reindex(index=truth.index).fillna(0.0)
                rt = time.perf_counter() - t0
                status_rows.append(dict(dataset=dataset, scenario=scen, method=name,
                                        status="executed", runtime_s=round(rt, 3), note=""))
                pred.to_csv(ed / f"pred_{name}.tsv", sep="\t")
                row = SC.score_all(truth, pred, mp, clusters, class_by)
                metrics_rows.append(dict(dataset=dataset, scenario=scen, method=name,
                                         runtime_s=round(rt, 3), n_samples=len(truth), **row))
                for cr in SC.per_celltype_rows(truth, pred, class_by, mp):
                    ct_rows.append(dict(dataset=dataset, scenario=scen, method=name, **cr))
            except Exception as e:  # honest failure, never fabricate
                status_rows.append(dict(dataset=dataset, scenario=scen, method=name,
                                        status="failed_run", runtime_s=round(time.perf_counter() - t0, 3),
                                        note=f"{type(e).__name__}: {str(e)[:160]}"))
        print(f"[{dataset}/{scen}] genes={cbulk.shape[0]} samples={truth.shape[1]}types "
              f"clusters={len(clusters)} classes={pd.Series(class_by).value_counts().to_dict()}")


def main(argv=None):
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument("--run-real-data", action="store_true")
    ap.add_argument("--datasets", nargs="+", default=["breast", "lung"])
    ap.add_argument("--seeds", type=int, nargs="+", default=[0, 1, 2])
    args = ap.parse_args(argv)
    import _harness as H
    if not H.real_data_enabled(args.run_real_data):
        print(f"Refusing without --run-real-data (or {H.REAL_DATA_ENV}=1).", file=sys.stderr); return 2
    OUT.mkdir(parents=True, exist_ok=True)
    mrows, ctrows, srows, certrows = [], [], [], []
    t0 = time.perf_counter()
    for ds in args.datasets:
        run(ds, args.seeds, mrows, ctrows, srows, certrows)
    pd.DataFrame(mrows).to_csv(OUT / "metrics_by_scenario.tsv", sep="\t", index=False)
    pd.DataFrame(ctrows).to_csv(OUT / "metrics_by_celltype.tsv", sep="\t", index=False)
    pd.DataFrame(srows).to_csv(OUT / "method_status.tsv", sep="\t", index=False)
    pd.DataFrame(certrows).to_csv(OUT / "identifiability_certificate.tsv", sep="\t", index=False)
    (OUT / "run_manifest.json").write_text(json.dumps(
        {"datasets": args.datasets, "seeds": args.seeds, "scenarios": [s[0] for s in SCENARIOS],
         "runtime_s": round(time.perf_counter() - t0, 1), "max_cells_per_type": MAX_CELLS_PER_TYPE},
        indent=2, default=str))
    m = pd.DataFrame(mrows)
    print("\n=== IN-PROCESS SUMMARY (mean over scenarios/datasets) ===")
    if len(m):
        show = m.groupby("method")[["broad_rmse", "broad_pearson", "fine_rmse", "fine_pearson",
                                    "condfam_rmse", "swap_rmse", "rare_recall", "absent_mass",
                                    "overconfident_error_rate", "runtime_s"]].mean().round(4)
        print(show.to_string())
    print(f"\nWrote {OUT}\nExports for external tools: {EXPORT}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
