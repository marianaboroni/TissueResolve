#!/usr/bin/env python
"""Isolated-env runner for Rectangle (rectanglepy) — benchmark comparison only.

Rectangle requires Python >=3.10 and cannot import into the TissueResolve
project venv (3.9), so — like cell2location — it runs in its own environment
(``benchmarks/envs/rectangle_py311``) invoked via subprocess. This script has NO
TissueResolve dependency; it consumes already-prepared inputs and writes
predictions + a metadata JSON, then the py39 side reads them back for scoring.

Fair-comparison protocol: Rectangle is given the SAME reference cells and the
SAME bulk samples as TissueResolve, but is entitled to its OWN gene/marker
selection (that is part of the method). It is run twice per call:
  * ``correct_mrna_bias=True``  — Rectangle's default; cell-fraction-comparable.
  * ``correct_mrna_bias=False`` — RNA-proportion-comparable (matches how
    TissueResolve reports mRNA proportions). Both are recorded; the scorer
    picks the quantity matching each scenario's ground-truth type.

Rectangle also emits a single scalar ``Unknown`` column (per-sample residual);
it is preserved in the output (never dropped) so unknown-content behaviour can
be compared honestly against TissueResolve's diagnostics.

Inputs are files (not Python objects) so the process boundary stays clean:
  --ref-h5ad     single-cell reference AnnData (raw counts preferred)
  --cell-type-col obs column with the cell-type / state labels
  --bulk-tpm     TSV, samples (rows) x genes (cols), TPM
  --out-dir      where predictions_*.tsv + rectangle_metadata.json are written

Usage (from repo root):
  benchmarks/envs/rectangle_py311/bin/python benchmarks/dev/run_rectangle.py \
      --ref-h5ad path/to/ref.h5ad --cell-type-col cell_type \
      --bulk-tpm bulk_tpm.tsv --out-dir out/ --gene-col feature_name
"""
from __future__ import annotations

import argparse
import json
import time
import warnings
from pathlib import Path


def _write_meta(out_dir: Path, status: str, **kw) -> None:
    base = dict(method="Rectangle", status=status, executed=False,
                rectanglepy_version="", runtime_seconds=0.0, n_cell_types=0,
                n_bulk_samples=0, n_shared_genes=0, correct_mrna_bias_variants=[],
                optimize_cutoffs=False, command_run="run_rectangle.py",
                warnings="", error_message="")
    base.update(kw)
    (out_dir / "rectangle_metadata.json").write_text(json.dumps(base, indent=2))


def _resolve_counts(adata):
    """Best-available raw counts, mirroring _harness.resolve_counts_matrix."""
    if "counts" in getattr(adata, "layers", {}):
        return adata.layers["counts"], "layers['counts']"
    raw = getattr(adata, "raw", None)
    if raw is not None and getattr(raw, "X", None) is not None \
            and raw.X.shape[1] == adata.n_vars:
        return raw.X, "raw.X"
    return adata.X, "X"


def main(argv=None) -> int:
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument("--ref-h5ad", required=True)
    ap.add_argument("--cell-type-col", default="cell_type")
    ap.add_argument("--bulk-tpm", required=True, help="TSV: samples x genes, TPM")
    ap.add_argument("--out-dir", required=True)
    ap.add_argument("--gene-col", default=None,
                    help="var column with gene symbols (else var_names)")
    ap.add_argument("--counts-layer", default=None)
    ap.add_argument("--cells-per-type", type=int, default=0,
                    help="deterministic subsample cap per cell type (0 = all)")
    ap.add_argument("--seed", type=int, default=0)
    ap.add_argument("--optimize-cutoffs", action="store_true",
                    help="Rectangle marker-cutoff grid search (slow; off by default)")
    ap.add_argument("--bootstraps", type=int, default=7)
    args = ap.parse_args(argv)

    out_dir = Path(args.out_dir)
    out_dir.mkdir(parents=True, exist_ok=True)
    t0 = time.perf_counter()

    try:
        import importlib.util
        if importlib.util.find_spec("rectanglepy") is None:
            _write_meta(out_dir, "skipped",
                        error_message="rectanglepy not installed in this env "
                                       "(see benchmarks/envs/install_rectangle.sh)")
            print("Rectangle: skipped (not installed)")
            return 0

        import numpy as np
        import pandas as pd
        import anndata as ad
        import rectanglepy as rp

        rp_version = rp.version("rectanglepy")

        # --- reference ----------------------------------------------------
        adata = ad.read_h5ad(args.ref_h5ad)
        if args.counts_layer:
            counts = adata.layers[args.counts_layer]
            counts_source = f"layers['{args.counts_layer}']"
        else:
            counts, counts_source = _resolve_counts(adata)
        # gene symbols on the columns so they can meet the bulk TPM genes
        if args.gene_col and args.gene_col in adata.var.columns:
            genes = adata.var[args.gene_col].astype(str).to_numpy()
        else:
            genes = adata.var_names.astype(str).to_numpy()

        labels = adata.obs[args.cell_type_col].astype(str).to_numpy()

        # optional deterministic per-type subsample (keeps runtime bounded)
        keep = np.arange(adata.n_obs)
        if args.cells_per_type and args.cells_per_type > 0:
            rng = np.random.default_rng(args.seed)
            sel = []
            for ct in pd.unique(labels):
                idx = np.where(labels == ct)[0]
                if len(idx) > args.cells_per_type:
                    idx = rng.choice(idx, args.cells_per_type, replace=False)
                sel.append(idx)
            keep = np.sort(np.concatenate(sel))

        X = counts[keep]
        X = np.asarray(X.todense()) if hasattr(X, "todense") else np.asarray(X)
        # de-duplicate gene columns (Rectangle needs unique gene index)
        gser = pd.Index(genes)
        if gser.duplicated().any():
            df = pd.DataFrame(X, columns=genes)
            df = df.T.groupby(level=0).sum().T
            X = df.to_numpy()
            genes = df.columns.to_numpy()
        sc_adata = ad.AnnData(X.astype("float32"))
        sc_adata.var_names = pd.Index(genes)
        sc_adata.obs[args.cell_type_col] = pd.Categorical(labels[keep])

        # --- bulk (TPM, samples x genes) ---------------------------------
        bulk = pd.read_csv(args.bulk_tpm, sep="\t", index_col=0)
        shared = [g for g in bulk.columns if g in set(genes)]
        n_shared = len(shared)
        if n_shared < 50:
            _write_meta(out_dir, "failed", rectanglepy_version=rp_version,
                        n_shared_genes=n_shared,
                        error_message=f"only {n_shared} shared genes ref/bulk")
            print(f"Rectangle: failed - only {n_shared} shared genes")
            return 0

        adv = rp.RectangleAdvancedParameters(number_of_bootstraps=args.bootstraps)
        with warnings.catch_warnings():
            warnings.simplefilter("ignore")
            sig = rp.pp.build_rectangle_signatures(
                sc_adata, cell_type_col=args.cell_type_col, bulks=bulk,
                optimize_cutoffs=bool(args.optimize_cutoffs),
                advanced_parameters=adv)
            variants = []
            for corr in (True, False):
                est = rp.tl.deconvolution(sig, bulk, correct_mrna_bias=corr)
                est_df = est[0] if isinstance(est, tuple) else est
                tag = "mrnacorrected" if corr else "raw"
                predf = out_dir / f"predictions_{tag}.tsv"
                est_df.round(6).to_csv(predf, sep="\t")
                variants.append(tag)

        _write_meta(out_dir, "executed", executed=True,
                    rectanglepy_version=rp_version,
                    runtime_seconds=time.perf_counter() - t0,
                    n_cell_types=int(len(sig.pseudobulk_sig_cpm.columns)),
                    n_bulk_samples=int(bulk.shape[0]),
                    n_shared_genes=n_shared,
                    correct_mrna_bias_variants=variants,
                    optimize_cutoffs=bool(args.optimize_cutoffs),
                    command_run=(f"rectanglepy build_rectangle_signatures + "
                                 f"deconvolution (bootstraps={args.bootstraps}, "
                                 f"optimize_cutoffs={bool(args.optimize_cutoffs)})"),
                    warnings=(f"counts_source={counts_source}; Rectangle emits a "
                              f"scalar Unknown column (preserved)"))
        print(f"Rectangle: executed -> {out_dir} (variants: {variants})")
        return 0
    except Exception as exc:  # noqa: BLE001 — never abort the outer benchmark
        import traceback
        _write_meta(out_dir, "failed", runtime_seconds=time.perf_counter() - t0,
                    error_message=f"{type(exc).__name__}: {exc}",
                    warnings=traceback.format_exc(limit=4))
        print(f"Rectangle: failed - {exc}")
        return 0


if __name__ == "__main__":
    raise SystemExit(main())
