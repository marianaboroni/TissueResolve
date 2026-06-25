#!/usr/bin/env python
"""
09 — Hierarchical (broad → fine) deconvolution on the breast-cancer harness.

Runs the new ``--resolution-mode hierarchical`` workflow on the *existing*
prepared reference, pseudobulk mixtures, and Visium section.  It does **not**
download or modify any data and does not alter the core algorithms.

What it does
------------
1. Loads the saved TissueResolve reference (fine ``cell_type`` labels).
2. Loads the documented fine→broad mapping
   (``config/breast_cancer_cell_type_hierarchy.tsv``) — the reference only has
   fine labels, so the broad families are supplied explicitly.
3. Saves the hierarchical reference bundle (broad + fine references, hierarchy
   tables, metadata).
4. Runs hierarchical **bulk** deconvolution on the pseudobulk mixtures and
   compares flat vs hierarchical against the known ground truth.
5. Runs hierarchical **spatial** deconvolution on the Visium section
   (optionally subsampling spots for tractability — the cap is logged).
6. Writes a reproducible family-aware colour map and a summary.

Outputs are written under ``outputs/hierarchical/``.

Usage
-----
    python scripts/09_run_hierarchical_deconvolution.py
    python scripts/09_run_hierarchical_deconvolution.py --max-spatial-spots 1200
    python scripts/09_run_hierarchical_deconvolution.py --skip-spatial
"""
from __future__ import annotations

import argparse
import json
import sys
import warnings
from pathlib import Path

import numpy as np
import pandas as pd

import _harness as H

HIER_DIR = H.OUTPUTS_DIR / "hierarchical"
HIER_REF_DIR = HIER_DIR / "reference"
HIERARCHY_TSV = H.HARNESS_DIR / "config" / "breast_cancer_cell_type_hierarchy.tsv"


def _families_of(mapping: dict[str, str], props: pd.DataFrame) -> pd.DataFrame:
    from tissueresolve.reference.hierarchy import aggregate_predictions_by_family
    return aggregate_predictions_by_family(props, mapping)


def main(argv: list[str] | None = None) -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--n-bootstrap", type=int, default=0)
    parser.add_argument("--max-spatial-spots", type=int, default=1500,
                        help="Subsample at most this many spots for the spatial "
                             "run (the cap is logged; use 0 for all spots).")
    parser.add_argument("--skip-spatial", action="store_true")
    parser.add_argument("--seed", type=int, default=0)
    args = parser.parse_args(argv)

    if not H.SAVED_REFERENCE_DIR.exists():
        print(f"ERROR: reference not found at {H.SAVED_REFERENCE_DIR}.  "
              "Run 01_prepare_reference.py first.", file=sys.stderr)
        return 1
    if not HIERARCHY_TSV.exists():
        print(f"ERROR: hierarchy mapping not found: {HIERARCHY_TSV}.", file=sys.stderr)
        return 1

    import tissueresolve as tr
    from tissueresolve.results import ReferenceSignature
    from tissueresolve.reference.hierarchy import (
        build_cell_type_hierarchy, load_hierarchy_mapping,
    )
    from tissueresolve.reference.hierarchical_build import (
        save_hierarchical_reference_bundle,
    )
    from tissueresolve.plotting import palette as pal

    HIER_DIR.mkdir(parents=True, exist_ok=True)

    ref = ReferenceSignature.load(H.SAVED_REFERENCE_DIR)
    raw_map = load_hierarchy_mapping(HIERARCHY_TSV)
    # restrict / validate against the reference's actual fine labels
    with warnings.catch_warnings():
        warnings.simplefilter("ignore")
        mapping = build_cell_type_hierarchy(list(ref.cell_types), raw_map)
    families = sorted(set(mapping.values()))
    print(f"Hierarchy: {len(ref.cell_types)} fine types → {len(families)} "
          f"families: {families}")

    # 1. save the hierarchical reference bundle
    with warnings.catch_warnings():
        warnings.simplefilter("ignore")
        save_hierarchical_reference_bundle(ref, mapping, HIER_REF_DIR)
    # also copy the resolved hierarchy to the top of the bundle
    from tissueresolve.reference.hierarchy import save_hierarchy_mapping
    save_hierarchy_mapping(mapping, HIER_DIR / "cell_type_hierarchy.tsv")

    summary: dict = {"families": families, "n_fine": len(ref.cell_types)}

    # ------------------------------------------------------------------
    # 2. Hierarchical BULK
    # ------------------------------------------------------------------
    bulk_metrics = {}
    if H.PSEUDOBULK_COUNTS.exists() and H.PSEUDOBULK_TRUE_PROPS.exists():
        from tissueresolve.io.bulk import read_bulk_counts
        from tissueresolve.bulk.hierarchical import save_hierarchical_bulk_outputs

        bulk_raw = read_bulk_counts(H.PSEUDOBULK_COUNTS)
        true_props = pd.read_csv(H.PSEUDOBULK_TRUE_PROPS, sep="\t",
                                 comment="#", index_col=0)
        bulk, orient = H.orient_bulk_genes_by_samples(bulk_raw, ref.gene_names)
        ref_genes = set(map(str, ref.gene_names))
        bulk = bulk.loc[[g for g in bulk.index if g in ref_genes]]
        print(f"\n[bulk] {bulk.shape[1]} samples × {bulk.shape[0]} shared genes")

        with warnings.catch_warnings():
            warnings.simplefilter("ignore")
            flat = tr.deconv_bulk(bulk, ref, resolution_mode="none",
                                  n_bootstrap=args.n_bootstrap)
            hres = tr.deconv_bulk(bulk, ref, resolution_mode="hierarchical",
                                  hierarchy_mapping=mapping,
                                  n_bootstrap=args.n_bootstrap)
        save_hierarchical_bulk_outputs(hres, HIER_DIR)

        # flat vs hierarchical at the fine level (shared columns only)
        m_flat = H.compute_bulk_metrics(true_props, flat.deconv.proportions)
        m_hier = H.compute_bulk_metrics(true_props, hres.estimates.combined_fine)
        # family-level comparison (aggregate true props to families)
        true_fam = _families_of(mapping, true_props)
        m_fam = H.compute_bulk_metrics(true_fam, hres.estimates.family_proportions)
        bulk_metrics = {"fine_flat": m_flat, "fine_hierarchical": m_hier,
                        "family_hierarchical": m_fam}
        H.write_json(bulk_metrics, HIER_DIR / "bulk_flat_vs_hierarchical_metrics.json")
        H.write_tsv(hres.estimates.qc.set_index("broad_family"),
                    HIER_DIR / "hierarchical_qc.tsv")
        print(f"[bulk] fine flat   : pearson={m_flat['pearson']:.3f} "
              f"rmse={m_flat['rmse']:.4f}")
        print(f"[bulk] fine hier   : pearson={m_hier['pearson']:.3f} "
              f"rmse={m_hier['rmse']:.4f}")
        print(f"[bulk] family hier : pearson={m_fam['pearson']:.3f} "
              f"rmse={m_fam['rmse']:.4f}")
        summary["bulk_metrics"] = bulk_metrics
        summary["bulk_unresolved_families"] = hres.estimates.unresolved_families
    else:
        print("[bulk] pseudobulk not found; skipping bulk hierarchical run.")

    # ------------------------------------------------------------------
    # 3. Hierarchical SPATIAL
    # ------------------------------------------------------------------
    if not args.skip_spatial and H.SPATIAL_H5AD.exists():
        import anndata as ad
        import scipy.sparse as sp
        from tissueresolve.spatial.hierarchical import save_hierarchical_spatial_outputs

        print(f"\n[spatial] loading {H.SPATIAL_H5AD}")
        adata = ad.read_h5ad(H.SPATIAL_H5AD)
        obs = adata.obs
        if not {"array_row", "array_col"} <= set(obs.columns):
            print("ERROR: array_row/array_col missing in Visium obs.",
                  file=sys.stderr)
            return 1

        # subsample spots for tractability (logged, never silent)
        n_total = adata.n_obs
        cap = args.max_spatial_spots
        if cap and n_total > cap:
            rng = np.random.default_rng(args.seed)
            keep = np.sort(rng.choice(n_total, size=cap, replace=False))
            adata = adata[keep].copy()
            print(f"[spatial] subsampled {cap}/{n_total} spots "
                  f"(seed={args.seed}); recorded in metadata.")
        else:
            cap = n_total

        mat, _ = H.resolve_counts_matrix(adata)
        Y = mat if sp.issparse(mat) else np.asarray(mat)
        array_row = adata.obs["array_row"].to_numpy()
        array_col = adata.obs["array_col"].to_numpy()
        lib_sizes = (np.asarray(Y.sum(axis=1)).ravel()
                     if sp.issparse(Y) else Y.sum(axis=1)).astype("float32")
        gene_names, _, _ = H.select_spatial_gene_identifiers(adata)
        spot_ids = list(adata.obs_names)

        overlap = H.gene_overlap(gene_names, ref.gene_names)
        print(f"[spatial] gene overlap: {overlap}")

        with warnings.catch_warnings():
            warnings.simplefilter("ignore")
            sres = tr.deconv_spatial(
                Y, ref, array_row, array_col, lib_sizes, gene_names, spot_ids,
                resolution_mode="hierarchical", hierarchy_mapping=mapping,
                run_neighbourhood=False)
        save_hierarchical_spatial_outputs(sres, HIER_DIR)
        H.write_json(
            {"n_spots_used": int(cap), "n_spots_total": int(n_total),
             "lambda_spatial_family": sres.run_metadata.get("lambda_spatial_family"),
             "lambda_spatial_fine": sres.run_metadata.get("lambda_spatial_fine"),
             "unresolved_families": sres.estimates.unresolved_families},
            HIER_DIR / "spatial_run_metadata.json")
        print(f"[spatial] done; unresolved families: "
              f"{sres.estimates.unresolved_families}")
        summary["spatial_unresolved_families"] = sres.estimates.unresolved_families
        summary["spatial_n_spots_used"] = int(cap)

        # colour map from the spatial fine columns (family-aware)
        cm = pal.build_hierarchical_color_map(
            list(sres.estimates.combined_fine.columns), mapping)
        pal.save_hierarchical_color_map(cm, HIER_DIR)
    else:
        print("[spatial] skipped.")
        # still emit a colour map from the reference's fine types
        cm = pal.build_hierarchical_color_map(list(ref.cell_types), mapping)
        pal.save_hierarchical_color_map(cm, HIER_DIR)

    # ------------------------------------------------------------------
    # 4. Summary
    # ------------------------------------------------------------------
    H.write_json(summary, HIER_DIR / "hierarchical_summary.json")
    _write_summary_md(summary, HIER_DIR / "hierarchical_summary.md")
    print(f"\nWrote hierarchical outputs to {HIER_DIR}/")
    return 0


def _write_summary_md(summary: dict, path: Path) -> None:
    lines = ["# Hierarchical broad-to-fine deconvolution — summary", ""]
    lines.append(f"- Fine cell types: **{summary.get('n_fine', '—')}**")
    lines.append(f"- Broad families: **{len(summary.get('families', []))}** "
                 f"({', '.join(summary.get('families', []))})")
    bm = summary.get("bulk_metrics")
    if bm:
        lines += ["", "## Bulk: flat vs hierarchical (vs ground truth)", "",
                  "| level | pearson | rmse | mae |",
                  "|---|---|---|---|"]
        for key, label in (("fine_flat", "fine (flat)"),
                           ("fine_hierarchical", "fine (hierarchical)"),
                           ("family_hierarchical", "family (hierarchical)")):
            m = bm.get(key, {})
            lines.append(f"| {label} | {m.get('pearson', float('nan')):.3f} | "
                         f"{m.get('rmse', float('nan')):.4f} | "
                         f"{m.get('mae', float('nan')):.4f} |")
        unres = summary.get("bulk_unresolved_families", [])
        lines += ["", f"- Unresolved families (bulk): "
                  f"{', '.join(unres) if unres else 'none'}"]
    if "spatial_unresolved_families" in summary:
        lines += ["", "## Spatial",
                  f"- Spots used: {summary.get('spatial_n_spots_used', '—')}",
                  f"- Unresolved families (spatial): "
                  f"{', '.join(summary['spatial_unresolved_families']) or 'none'}"]
    lines += ["", "_Unresolved families are reported at the broad level only; "
              "their subtype splits are not claimed._", ""]
    path.write_text("\n".join(lines), encoding="utf-8")


if __name__ == "__main__":
    raise SystemExit(main())
