#!/usr/bin/env python
"""cell2location on the synthetic-spatial scenario (CPU). Honest status.

Two-stage: (1) regression model estimates per-cell-type reference signatures from
the sc reference; (2) Cell2location maps them to spots.  Epochs are reduced for a
CPU benchmark (recorded in status).  Outputs spot×celltype proportions
(q05 cell-abundance, row-normalised).  Run with the c2l_py39 venv's python.

Usage:  benchmarks/envs/c2l_py39/bin/python benchmarks/spatial/methods/run_cell2location_synthetic.py
"""
from __future__ import annotations

import json
import time
from pathlib import Path

import numpy as np
import pandas as pd

REPO = Path(__file__).resolve().parents[3]
BASE = REPO / "benchmarks" / "outputs" / "spatial_real_visium" / "external_inputs"
REF_EPOCHS = 250
SP_EPOCHS = 6000          # reduced for CPU; default is 30000


def _status(st, ver, msg, runtime, predf=""):
    pd.DataFrame([{"method": "cell2location", "status": st, "version": ver,
                   "runtime_seconds": round(runtime, 1), "ref_epochs": REF_EPOCHS,
                   "spatial_epochs": SP_EPOCHS, "error": str(msg)[:200], "output": predf}]
                 ).to_csv(BASE / "cell2location_status.tsv", sep="\t", index=False)
    print("cell2location:", st, ver, msg)


def main() -> int:
    t0 = time.time()
    try:
        import anndata as ad
        import cell2location
        from cell2location.models import RegressionModel, Cell2location
        ver = cell2location.__version__
    except Exception as exc:  # noqa: BLE001
        _status("skipped", "", f"import failed: {exc}", time.time() - t0)
        return 0
    try:
        sc_counts = pd.read_csv(BASE / "reference" / "reference_counts_genes_by_cells.tsv",
                                sep="\t", index_col=0)              # genes × cells
        meta = pd.read_csv(BASE / "reference" / "reference_cell_metadata.tsv", sep="\t")
        spots = pd.read_csv(BASE / "spot_counts_genes_by_spots.tsv",
                            sep="\t", index_col=0)                   # genes × spots
        shared = sc_counts.index.intersection(spots.index)
        sc_counts, spots = sc_counts.loc[shared], spots.loc[shared]

        adata_ref = ad.AnnData(sc_counts.T.values.astype("float32"),
                               obs=pd.DataFrame({"cellType": meta["cellType"].values},
                                                index=meta["cell_id"].values),
                               var=pd.DataFrame(index=shared.astype(str)))
        adata_st = ad.AnnData(spots.T.values.astype("float32"),
                              obs=pd.DataFrame(index=spots.columns.astype(str)),
                              var=pd.DataFrame(index=shared.astype(str)))

        # stage 1 — reference signatures
        RegressionModel.setup_anndata(adata=adata_ref, labels_key="cellType")
        rm = RegressionModel(adata_ref)
        rm.train(max_epochs=REF_EPOCHS)
        adata_ref = rm.export_posterior(adata_ref,
                                        sample_kwargs={"num_samples": 200})
        means = adata_ref.varm["means_per_cluster_mu_fg"] \
            if "means_per_cluster_mu_fg" in adata_ref.varm else \
            adata_ref.var[[c for c in adata_ref.var.columns if "means_per_cluster" in c]]
        inf_aver = pd.DataFrame(means, index=adata_ref.var_names)
        inf_aver.columns = [c.replace("means_per_cluster_mu_fg_", "") for c in inf_aver.columns]

        # stage 2 — spatial mapping
        common = adata_st.var_names.intersection(inf_aver.index)
        adata_st = adata_st[:, common].copy()
        inf_aver = inf_aver.loc[common]
        Cell2location.setup_anndata(adata=adata_st)
        mod = Cell2location(adata_st, cell_state_df=inf_aver,
                            N_cells_per_location=30, detection_alpha=20)
        mod.train(max_epochs=SP_EPOCHS, batch_size=None, train_size=1)
        adata_st = mod.export_posterior(adata_st,
                                        sample_kwargs={"num_samples": 200})
        ab = adata_st.obsm["q05_cell_abundance_w_sf"]
        ab = pd.DataFrame(np.asarray(ab), index=adata_st.obs_names,
                          columns=[c.replace("q05cell_abundance_w_sf_", "")
                                   for c in adata_st.uns["mod"]["factor_names"]]) \
            if not isinstance(ab, pd.DataFrame) else ab
        ab.columns = [str(c).split("_")[-1] if "abundance" in str(c) else str(c) for c in ab.columns]
        props = ab.div(ab.sum(axis=1), axis=0)
        predf = str(BASE / "cell2location_pred.tsv")
        props.round(6).to_csv(predf, sep="\t")
        _status("executed", ver, f"{props.shape[1]} cell types", time.time() - t0, predf)
        return 0
    except Exception as exc:  # noqa: BLE001
        import traceback; traceback.print_exc()
        _status("failed", ver, exc, time.time() - t0)
        return 0


if __name__ == "__main__":
    raise SystemExit(main())
