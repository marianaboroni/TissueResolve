#!/usr/bin/env python
"""P1 — unified external bulk comparison on identical donor-held-out inputs.

Scores ALL methods on the SAME exported breast pseudobulk scenarios (donor-disjoint,
seed 0) used by the external-tool pipeline, on the focus metrics: conditional
within-family RMSE, paired rare precision/recall, false-positive subtype rate,
absent subtype mass, effective-N, runtime.

Methods:
  * external: MuSiC, BisqueRNA  (predictions already produced by
    benchmarks/bulk/run_external_holdout.R — loaded here, never refit/fabricated)
  * TissueResolve: wNNLS (default), Poisson GLM, NB GLM, external plain-NNLS control
    (run here on the identical reference + bulk)

BayesPrism / DWLS / SCDC are NOT installed in this environment → recorded as
deferred (with install pointers), never fabricated. CARD is spatial and excluded.

Reuses the exported inputs under benchmarks/outputs/holdout_bulk/external_inputs/.

Usage:
  PYTHONPATH=src:. python benchmarks/diagnostics/external_bulk_comparison.py --run-real-data
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
from benchmarks.diagnostics.nb_bulk_solver_benchmark import _score  # noqa: E402

EXT = REPO / "benchmarks" / "outputs" / "holdout_bulk" / "external_inputs"
OUT = REPO / "benchmarks" / "outputs" / "nb_bulk_solver" / "external_comparison"
HIERARCHY_TSV = (REPO / "examples" / "real_breast_cancer" / "config" /
                 "breast_cancer_cell_type_hierarchy.tsv")
MIN_CELLS = 30
EXTERNAL = ["MuSiC", "BisqueRNA"]
TR_METHODS = ["wNNLS", "poisson_glm_experimental", "nb_glm_experimental", "external_nnls_control"]
DEFERRED = ["BayesPrism", "DWLS", "SCDC"]


def _build_ref(adata, ref_mask, drop_type=None):
    import warnings as _w
    mask = ref_mask.copy()
    if drop_type is not None:
        mask = mask & (adata.obs["cell_type"].astype(str) != drop_type).to_numpy()
    with _w.catch_warnings():
        _w.simplefilter("ignore")
        return H.prepare_reference(adata[mask].copy(), min_cells=MIN_CELLS,
                                   estimate_overdispersion=True).reference


def _tr_predict(method, bulk, ref):
    import tissueresolve as tr
    from tissueresolve.config import TissueResolveConfig
    if method == "external_nnls_control":
        from scipy.optimize import nnls
        genes = sorted(set(bulk.index) & set(ref.gene_names))
        Phi = (ref.as_R_cpm()[:, [list(ref.gene_names).index(g) for g in genes]].T / 1e6)
        B = bulk.loc[genes].to_numpy(float)
        Bn = B / np.where(B.sum(0, keepdims=True) > 0, B.sum(0, keepdims=True), 1.0)
        Pn = Phi / np.where(Phi.sum(0, keepdims=True) > 0, Phi.sum(0, keepdims=True), 1.0)
        out = np.zeros((bulk.shape[1], len(ref.cell_types)))
        for n in range(bulk.shape[1]):
            theta, _ = nnls(Pn, Bn[:, n]); s = theta.sum()
            out[n] = theta / s if s > 0 else np.full(len(ref.cell_types), 1.0 / len(ref.cell_types))
        return pd.DataFrame(out, index=bulk.columns, columns=list(ref.cell_types))
    cfg = TissueResolveConfig()
    cfg.bulk_solver.method = method
    return tr.deconv_bulk(bulk, ref, config=cfg, resolution_mode="none").deconv.proportions


def main(argv=None) -> int:
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument("--run-real-data", action="store_true")
    args = ap.parse_args(argv)
    if not H.real_data_enabled(args.run_real_data):
        print(f"Refusing without --run-real-data (or {H.REAL_DATA_ENV}=1).", file=sys.stderr)
        return 2
    if not (EXT / "export_manifest.json").exists():
        print(f"Missing exported external inputs at {EXT}. Run "
              "benchmarks/bulk/export_external_inputs.py + run_external_holdout.R first.",
              file=sys.stderr)
        return 3

    import anndata as ad
    from tissueresolve.reference.hierarchy import load_hierarchy_mapping, build_cell_type_hierarchy

    manifest = json.loads((EXT / "export_manifest.json").read_text())
    rare = manifest["rare_type"]
    missing_type = manifest.get("missing_type")
    OUT.mkdir(parents=True, exist_ok=True)

    adata = ad.read_h5ad(H.REFERENCE_H5AD)
    if "feature_name" in adata.var.columns:
        fn = adata.var["feature_name"].astype(str)
        if fn.nunique() == adata.n_vars:
            adata.var_names = fn.to_numpy()
    ref_donors = set(manifest["reference_donors"])
    ref_mask = adata.obs["donor_id"].astype(str).isin(ref_donors).to_numpy()
    ref_full = _build_ref(adata, ref_mask)
    ref_missing = _build_ref(adata, ref_mask, drop_type=missing_type) if missing_type else ref_full
    ref_types = [str(c) for c in ref_full.cell_types]
    mapping = build_cell_type_hierarchy(ref_types, load_hierarchy_mapping(HIERARCHY_TSV))

    rows, runtime_rows = [], []
    scen_dirs = sorted([d.name for d in EXT.iterdir() if d.is_dir()
                        and (d / "truth_mrna.tsv").exists()])
    for scen in scen_dirs:
        truth = pd.read_csv(EXT / scen / "truth_mrna.tsv", sep="\t", index_col=0)
        bulk = pd.read_csv(EXT / scen / "bulk_counts_genes_by_samples.tsv", sep="\t", index_col=0)
        ref = ref_missing if scen == "missing_population" else ref_full
        # external predictions (already produced by R; load, never refit)
        for method in EXTERNAL:
            pf = EXT / scen / f"{method}_pred.tsv"
            if not pf.exists():
                runtime_rows.append({"scenario": scen, "method": method,
                                     "status": "skipped: no prediction file"})
                continue
            pred = pd.read_csv(pf, sep="\t", index_col=0)
            m = _score(truth, pred, mapping, rare)
            rows.append({"scenario": scen, "method": method, "source": "external", **m})
        # TissueResolve methods on the identical reference + bulk
        for method in TR_METHODS:
            t0 = time.perf_counter()
            try:
                with warnings.catch_warnings():
                    warnings.simplefilter("ignore")
                    pred = _tr_predict(method, bulk, ref)
            except Exception as exc:  # noqa: BLE001
                runtime_rows.append({"scenario": scen, "method": method,
                                     "status": f"failed: {exc}"})
                print(f"  {scen} {method} FAILED: {exc}")
                continue
            rt = float(time.perf_counter() - t0)
            m = _score(truth, pred, mapping, rare)
            rows.append({"scenario": scen, "method": method, "source": "tissueresolve", **m})
            runtime_rows.append({"scenario": scen, "method": method,
                                 "runtime_seconds": round(rt, 2), "status": "ok"})
            print(f"  {scen:<18} {method:<26} fine={m['fine_pearson']:.3f} "
                  f"cond={m['conditional_rmse']:.3f} rareR={m['rare_sensitivity']} "
                  f"rareP={m['rare_precision']} effN={m['effective_n_pred']:.1f}/{m['effective_n_truth']:.1f}")

    per_scen = pd.DataFrame(rows)
    per_scen.to_csv(OUT / "per_scenario_metrics.tsv", sep="\t", index=False)
    pd.DataFrame(runtime_rows).to_csv(OUT / "runtime_metrics.tsv", sep="\t", index=False)
    focus = ["fine_pearson", "broad_pearson", "conditional_rmse", "rare_sensitivity",
             "rare_precision", "rare_fpr", "false_positive_subtype_rate",
             "absent_subtype_mass", "effective_n_pred", "effective_n_truth"]
    by_method = per_scen.groupby("method")[focus].mean(numeric_only=True).reset_index()
    by_method.to_csv(OUT / "by_method_focus_metrics.tsv", sep="\t", index=False)

    status = [{"method": m, "kind": "external", "status": "executed (predictions reused)"}
              for m in EXTERNAL]
    status += [{"method": m, "kind": "tissueresolve",
                "status": "executed"} for m in TR_METHODS]
    status += [{"method": m, "kind": "external", "status":
                "deferred: not installed in this environment (see "
                "benchmarks/envs/install_external_tools.sh)"} for m in DEFERRED]
    pd.DataFrame(status).to_csv(OUT / "method_status.tsv", sep="\t", index=False)

    (OUT / "manifest.json").write_text(json.dumps({
        "scenarios": scen_dirs, "seed": manifest.get("seed", 0), "dataset": "breast",
        "external_executed": EXTERNAL, "external_deferred": DEFERRED,
        "tissueresolve_methods": TR_METHODS, "rare_type": rare,
        "note": "identical donor-disjoint reference + bulk for all methods; "
                "external preds reused from run_external_holdout.R; CARD excluded (spatial)",
    }, indent=2), encoding="utf-8")
    print(f"\nWrote -> {OUT}/")
    print(by_method.round(3).to_string(index=False))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
