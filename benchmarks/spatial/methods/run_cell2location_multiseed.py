#!/usr/bin/env python
"""cell2location across all seed<s>/external_inputs/ from run_spatial_multiseed.py (CPU).

Run with the c2l_py39 venv python.  Reduced epochs for CPU (recorded).
"""
from __future__ import annotations

import time
from pathlib import Path

import numpy as np
import pandas as pd

REPO = Path(__file__).resolve().parents[3]
BASE = REPO / "benchmarks" / "outputs" / "spatial_multiseed"
REF_EPOCHS = 250
SP_EPOCHS = 4000          # reduced for CPU multi-seed; default 30000


def run_one(ei: Path, predf: Path):
    import anndata as ad
    from cell2location.models import RegressionModel, Cell2location
    sc = pd.read_csv(ei / "reference" / "reference_counts_genes_by_cells.tsv", sep="\t", index_col=0)
    meta = pd.read_csv(ei / "reference" / "reference_cell_metadata.tsv", sep="\t")
    spots = pd.read_csv(ei / "spatial_spot_counts.tsv", sep="\t", index_col=0)
    shared = sc.index.intersection(spots.index)
    sc, spots = sc.loc[shared], spots.loc[shared]
    aref = ad.AnnData(sc.T.values.astype("float32"),
                      obs=pd.DataFrame({"cellType": meta["cellType"].values},
                                       index=meta["cell_id"].values),
                      var=pd.DataFrame(index=shared.astype(str)))
    ast = ad.AnnData(spots.T.values.astype("float32"),
                     obs=pd.DataFrame(index=spots.columns.astype(str)),
                     var=pd.DataFrame(index=shared.astype(str)))
    RegressionModel.setup_anndata(adata=aref, labels_key="cellType")
    rm = RegressionModel(aref); rm.train(max_epochs=REF_EPOCHS)
    aref = rm.export_posterior(aref, sample_kwargs={"num_samples": 200})
    means = aref.varm["means_per_cluster_mu_fg"] if "means_per_cluster_mu_fg" in aref.varm else \
        aref.var[[c for c in aref.var.columns if "means_per_cluster" in c]]
    inf = pd.DataFrame(means, index=aref.var_names)
    inf.columns = [c.replace("means_per_cluster_mu_fg_", "") for c in inf.columns]
    common = ast.var_names.intersection(inf.index)
    ast = ast[:, common].copy(); inf = inf.loc[common]
    Cell2location.setup_anndata(adata=ast)
    mod = Cell2location(ast, cell_state_df=inf, N_cells_per_location=30, detection_alpha=20)
    mod.train(max_epochs=SP_EPOCHS, batch_size=None, train_size=1)
    ast = mod.export_posterior(ast, sample_kwargs={"num_samples": 200})
    ab = ast.obsm["q05_cell_abundance_w_sf"]
    if not isinstance(ab, pd.DataFrame):
        ab = pd.DataFrame(np.asarray(ab), index=ast.obs_names,
                          columns=[str(c).split("_")[-1] for c in ast.uns["mod"]["factor_names"]])
    ab.columns = [str(c).split("_")[-1] if "abundance" in str(c) else str(c) for c in ab.columns]
    props = ab.div(ab.sum(axis=1), axis=0)
    props.round(6).to_csv(predf, sep="\t")
    return props.shape[1]


def main() -> int:
    import cell2location
    ver = cell2location.__version__
    rows = []
    for sdir in sorted(BASE.glob("seed*")):
        ei = sdir / "external_inputs"
        predf = sdir / "cell2location_pred.tsv"
        t0 = time.time()
        try:
            nct = run_one(ei, predf)
            st = "executed"; err = ""
        except Exception as exc:  # noqa: BLE001
            import traceback; traceback.print_exc()
            st = "failed"; err = str(exc)[:200]; nct = 0
        rt = time.time() - t0
        rows.append({"seed": sdir.name, "method": "cell2location", "status": st,
                     "version": ver, "runtime_seconds": round(rt, 1),
                     "ref_epochs": REF_EPOCHS, "spatial_epochs": SP_EPOCHS, "error": err})
        print(f"{sdir.name} cell2location {st} {rt:.0f}s", flush=True)
    pd.DataFrame(rows).to_csv(BASE / "cell2location_multiseed_status.tsv", sep="\t", index=False)
    print("Wrote cell2location_multiseed_status.tsv")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
