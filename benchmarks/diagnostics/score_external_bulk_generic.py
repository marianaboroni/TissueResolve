#!/usr/bin/env python
"""Score a generic external-bulk base dir (exported by export_external_bulk_generic.py
+ run_external_bulk_generic.R): all methods on identical inputs, focus metrics.

Rebuilds the TissueResolve reference from the SAME exported single cells MuSiC/Bisque
used (fair apples-to-apples), runs the TR solvers (wNNLS / Poisson GLM / NB GLM /
external NNLS) on each scenario's bulk, loads the external predictions, and scores
everything on the focus metrics, aggregated per method (mean over scenario×seed).

Usage:
  PYTHONPATH=src:. python benchmarks/diagnostics/score_external_bulk_generic.py \
      --run-real-data --base benchmarks/outputs/external_bulk/lung
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
from benchmarks.diagnostics.nb_bulk_solver_benchmark import _score  # noqa: E402

EXTERNAL = ["MuSiC", "BisqueRNA"]
TR_METHODS = ["wNNLS", "poisson_glm_experimental", "nb_glm_experimental", "external_nnls_control"]
DEFERRED = ["BayesPrism", "DWLS", "SCDC"]


def _rebuild_reference(base: Path):
    """Rebuild a TissueResolve reference from the exported subsampled cells."""
    import anndata as ad
    rc = pd.read_csv(base / "reference" / "reference_counts_genes_by_cells.tsv",
                     sep="\t", index_col=0)                       # genes × cells
    meta = pd.read_csv(base / "reference" / "reference_cell_metadata.tsv", sep="\t")
    X = rc.to_numpy(float).T                                       # cells × genes
    obs = pd.DataFrame({"cell_type": meta["cellType"].astype(str).to_numpy(),
                        "donor_id": meta["SubjectName"].astype(str).to_numpy()},
                       index=meta["cell_id"].astype(str).to_numpy())
    adata = ad.AnnData(X=X, obs=obs, var=pd.DataFrame(index=rc.index.astype(str)))
    with warnings.catch_warnings():
        warnings.simplefilter("ignore")
        ref = H.prepare_reference(adata, min_cells=1, cell_type_col="cell_type",
                                  estimate_overdispersion=True).reference
    return ref


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
    cfg = TissueResolveConfig(); cfg.bulk_solver.method = method
    return tr.deconv_bulk(bulk, ref, config=cfg, resolution_mode="none").deconv.proportions


def main(argv=None) -> int:
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument("--run-real-data", action="store_true")
    ap.add_argument("--base", required=True)
    args = ap.parse_args(argv)
    if not H.real_data_enabled(args.run_real_data):
        print(f"Refusing without --run-real-data (or {H.REAL_DATA_ENV}=1).", file=sys.stderr)
        return 2
    base = Path(args.base)
    manifest = json.loads((base / "manifest.json").read_text())
    mapping = {str(k): str(v) for k, v in manifest["mapping"].items()}
    rare = manifest["rare_type"]
    ref = _rebuild_reference(base)
    print(f"[{manifest['dataset']}] rebuilt reference: {ref.n_cell_types} types, "
          f"{ref.n_genes} genes; rare={rare}")

    rows, runtime_rows = [], []
    scen_dirs = [d for d in sorted(base.iterdir()) if (d / "truth_mrna.tsv").exists()]
    for d in scen_dirs:
        scen = d.name
        truth = pd.read_csv(d / "truth_mrna.tsv", sep="\t", index_col=0)
        bulk = pd.read_csv(d / "bulk_counts_genes_by_samples.tsv", sep="\t", index_col=0)
        for method in EXTERNAL:
            pf = d / f"{method}_pred.tsv"
            if not pf.exists():
                continue
            pred = pd.read_csv(pf, sep="\t", index_col=0)
            rows.append({"scenario": scen, "method": method, "source": "external",
                         **_score(truth, pred, mapping, rare)})
        for method in TR_METHODS:
            t0 = time.perf_counter()
            try:
                with warnings.catch_warnings():
                    warnings.simplefilter("ignore")
                    pred = _tr_predict(method, bulk, ref)
            except Exception as exc:  # noqa: BLE001
                runtime_rows.append({"scenario": scen, "method": method,
                                     "status": f"failed: {exc}"})
                continue
            rt = float(time.perf_counter() - t0)
            rows.append({"scenario": scen, "method": method, "source": "tissueresolve",
                         **_score(truth, pred, mapping, rare)})
            runtime_rows.append({"scenario": scen, "method": method,
                                 "runtime_seconds": round(rt, 2), "status": "ok"})
        print(f"  scored {scen}")

    per_scen = pd.DataFrame(rows)
    out = base / "comparison"
    out.mkdir(parents=True, exist_ok=True)
    per_scen.to_csv(out / "per_scenario_metrics.tsv", sep="\t", index=False)
    pd.DataFrame(runtime_rows).to_csv(out / "runtime_metrics.tsv", sep="\t", index=False)
    focus = ["fine_pearson", "broad_pearson", "conditional_rmse", "rare_sensitivity",
             "rare_precision", "rare_fpr", "false_positive_subtype_rate",
             "absent_subtype_mass", "effective_n_pred", "effective_n_truth"]
    by_method = per_scen.groupby("method")[focus].mean(numeric_only=True).reset_index()
    by_method.to_csv(out / "by_method_focus_metrics.tsv", sep="\t", index=False)
    ext_status = base / "external_method_status.tsv"
    status_note = pd.read_csv(ext_status, sep="\t") if ext_status.exists() else pd.DataFrame()
    status_note.to_csv(out / "external_method_status.tsv", sep="\t", index=False)
    pd.DataFrame([{"method": m, "status": "deferred: not installed (compilation/deps)"}
                  for m in DEFERRED]).to_csv(out / "deferred_methods.tsv", sep="\t", index=False)
    print(f"\nWrote -> {out}/")
    print(by_method.round(3).to_string(index=False))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
