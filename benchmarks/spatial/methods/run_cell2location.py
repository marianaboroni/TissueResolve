#!/usr/bin/env python
"""
Real cell2location spatial deconvolution runner.

Trains the cell2location regression model on the harmonised single-cell
reference, then the spatial model on the exported Visium counts, and writes
per-spot cell-type proportions. Runs on CPU here (no GPU); --fast uses few
epochs (approximate — recorded in metadata). Real Visium has no ground truth,
so downstream comparison is concordance/structure, not accuracy. Fails
gracefully without aborting the benchmark.
"""
from __future__ import annotations
import argparse, importlib.util, json, time, warnings
from pathlib import Path

OUT = Path("benchmarks/outputs/spatial")
PREP = Path("benchmarks/outputs/prepared_inputs")


def _meta(status, **kw):
    d = OUT / "method_metadata"; d.mkdir(parents=True, exist_ok=True)
    base = dict(method="cell2location", version="", executed=False, imported=False,
                exported_only=False, status=status, runtime_seconds=0.0,
                input_normalization_used="counts", reference_level_used="fine",
                command_run="run_cell2location.py", warnings="", error_message="",
                output_path="")
    base.update(kw)
    (d / "cell2location.json").write_text(json.dumps(base, indent=2))


def main(argv=None) -> int:
    ap = argparse.ArgumentParser()
    ap.add_argument("--fast", action="store_true")
    args = ap.parse_args(argv)
    t0 = time.perf_counter()

    if importlib.util.find_spec("cell2location") is None:
        _meta("skipped", error_message="cell2location not installed (see benchmarks/envs/)")
        print("cell2location: skipped (not installed)"); return 0

    sc_counts = PREP / "reference" / "reference_counts_genes_by_cells.tsv"
    sc_meta = PREP / "reference" / "reference_cell_metadata.tsv"
    sp_counts = PREP / "spatial" / "spatial_counts_genes_by_spots.tsv"
    for f in (sc_counts, sc_meta, sp_counts):
        if not f.exists():
            _meta("failed", error_message=f"missing prepared input: {f}")
            print(f"cell2location: failed - missing {f}"); return 0

    try:
        import numpy as np, pandas as pd, anndata as ad
        import cell2location
        from cell2location.models import RegressionModel, Cell2location
        import torch
        gpu = bool(torch.cuda.is_available())
        acc = "gpu" if gpu else "cpu"
        reg_epochs = 50 if args.fast else 250
        c2l_epochs = 200 if args.fast else 3000

        ref = pd.read_csv(sc_counts, sep="\t", index_col=0)          # genes × cells
        meta = pd.read_csv(sc_meta, sep="\t")
        ad_ref = ad.AnnData(X=ref.T.to_numpy(dtype="float32"),
                            obs=pd.DataFrame({"cell_type": meta["cellType"].to_numpy(),
                                              "batch": meta["SubjectName"].astype(str).to_numpy()},
                                             index=list(ref.columns)),
                            var=pd.DataFrame(index=list(ref.index)))
        # regression model → per-cell-type reference signatures
        RegressionModel.setup_anndata(ad_ref, batch_key="batch", labels_key="cell_type")
        rm = RegressionModel(ad_ref)
        with warnings.catch_warnings():
            warnings.simplefilter("ignore")
            rm.train(max_epochs=reg_epochs, accelerator=acc)
            ad_ref = rm.export_posterior(ad_ref, sample_kwargs={"num_samples": 100,
                                          "batch_size": 2500, "accelerator": acc})
        # per-cluster signature
        if "means_per_cluster_mu_fg" in ad_ref.varm:
            inf = ad_ref.varm["means_per_cluster_mu_fg"].copy()
        else:
            inf = ad_ref.var[[c for c in ad_ref.var.columns
                              if c.startswith("means_per_cluster_mu_fg")]].copy()
        inf.columns = [c.replace("means_per_cluster_mu_fg_", "") for c in inf.columns]

        # spatial
        sp = pd.read_csv(sp_counts, sep="\t", index_col=0)            # genes × spots
        ad_sp = ad.AnnData(X=sp.T.to_numpy(dtype="float32"),
                           obs=pd.DataFrame(index=list(sp.columns)),
                           var=pd.DataFrame(index=list(sp.index)))
        shared = [g for g in ad_sp.var_names if g in set(inf.index)]
        ad_sp = ad_sp[:, shared].copy(); inf_s = inf.loc[shared]
        Cell2location.setup_anndata(ad_sp)
        cm = Cell2location(ad_sp, cell_state_df=inf_s, N_cells_per_location=30,
                           detection_alpha=20)
        with warnings.catch_warnings():
            warnings.simplefilter("ignore")
            cm.train(max_epochs=c2l_epochs, accelerator=acc)
            ad_sp = cm.export_posterior(ad_sp, sample_kwargs={"num_samples": 100,
                                        "batch_size": ad_sp.n_obs, "accelerator": acc})
        key = "q05_cell_abundance_w_sf"
        abund = ad_sp.obsm[key] if key in ad_sp.obsm else ad_sp.obsm[list(ad_sp.obsm)[0]]
        abund.columns = [str(c).split("q05_cell_abundance_w_sf_")[-1] for c in abund.columns]
        props = abund.div(abund.sum(axis=1), axis=0).fillna(0.0)
        OUT.joinpath("predictions").mkdir(parents=True, exist_ok=True)
        predf = OUT / "predictions" / "cell2location.tsv"
        props.round(6).to_csv(predf, sep="\t")
        _meta("executed", executed=True, version=cell2location.__version__,
              runtime_seconds=time.perf_counter() - t0,
              command_run=f"cell2location RegressionModel+Cell2location ({acc}, "
                          f"reg={reg_epochs}/c2l={c2l_epochs} epochs)",
              warnings=f"CPU/fast approximate (gpu={gpu}); no ground truth",
              output_path=str(predf),
              input_normalization_used="counts (abundance normalised per spot)")
        print(f"cell2location: executed -> {predf}")
    except Exception as exc:  # noqa: BLE001
        _meta("failed", runtime_seconds=time.perf_counter() - t0,
              error_message=f"{type(exc).__name__}: {exc}")
        print(f"cell2location: failed - {exc}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
