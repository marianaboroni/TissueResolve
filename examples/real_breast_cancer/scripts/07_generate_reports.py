#!/usr/bin/env python
"""
07 — Generate publication HTML reports + figures from existing harness outputs.

Reads the validation outputs already on disk (no re-download, no re-fit):

    outputs/reference/   reference summary, cell-type counts, gene-id source
    outputs/bulk/        bulk estimates + QC + warnings
    outputs/spatial/     spot proportions + QC + Moran's I + warnings
    outputs/resolution/  separability / spillover (optional)

and produces:

    outputs/bulk/figures/*        + outputs/bulk/report.html
    outputs/spatial/figures/*     + outputs/spatial/report.html
    outputs/validation_summary/report.html   (combined)

Every figure also writes its source data (.data.tsv).  Static PDF/SVG/PNG are
written only if kaleido is installed; otherwise interactive HTML + data are
saved and a warning is recorded.

Usage
-----
    python scripts/07_generate_reports.py
"""
from __future__ import annotations

import argparse
import json
import shutil
import warnings as _warnings
from pathlib import Path
from types import SimpleNamespace

import pandas as pd

import _harness as H


def _read_tsv(p: Path):
    try:
        return pd.read_csv(p, sep="\t", comment="#", index_col=0)
    except Exception:
        return None


def _copy_into_tables(results_dir: Path, files: list[Path]) -> None:
    tables = results_dir / "tables"
    tables.mkdir(parents=True, exist_ok=True)
    for f in files:
        if f and Path(f).exists():
            shutil.copy(Path(f), tables / Path(f).name)


def _load_reference():
    from tissueresolve.results import ReferenceSignature

    if H.SAVED_REFERENCE_DIR.exists():
        return ReferenceSignature.load(H.SAVED_REFERENCE_DIR)
    return None


def _write_report_bundle_artifacts(
    bulk_warns: list[str],
    spatial_warns: list[str],
    bulk_meta: dict[str, object],
    spatial_meta: dict[str, object],
) -> None:
    from tissueresolve.report import methods_text
    from tissueresolve.results import BulkDeconvResult, QCReport, SpatialDeconvResult

    H.OUT_SUMMARY_DIR.mkdir(parents=True, exist_ok=True)
    H.write_json({"bulk": bulk_warns, "spatial": spatial_warns},
                 H.OUT_SUMMARY_DIR / "warnings.json")
    run_metadata = {"generated_by": "07_generate_reports.py"}
    if bulk_meta:
        run_metadata["bulk"] = bulk_meta
    if spatial_meta:
        run_metadata["spatial"] = spatial_meta
    H.write_json(run_metadata, H.OUT_SUMMARY_DIR / "run_metadata.json")

    methods: list[str] = []
    # Bulk methods text is inferred from the saved outputs.
    try:
        bulk_df = _read_tsv(H.OUT_BULK_DIR / "bulk_estimated_proportions.tsv")
        bulk_qc = _read_tsv(H.OUT_BULK_DIR / "bulk_qc.tsv")
        if bulk_df is not None:
            recon_r2 = None
            profile_corr = None
            mismatch_flag = None
            if bulk_qc is not None:
                if "recon_r2" in bulk_qc.columns:
                    recon_r2 = bulk_qc["recon_r2"]
                if "profile_corr" in bulk_qc.columns:
                    profile_corr = bulk_qc["profile_corr"]
                if "mismatch_flag" in bulk_qc.columns:
                    mismatch_flag = bulk_qc["mismatch_flag"]
            qc = QCReport(
                modality="bulk",
                recommendations=[],
                metadata={},
                recon_r2=recon_r2,
                profile_corr=profile_corr,
                mismatch_flag=mismatch_flag,
            )
            coverage_r2 = recon_r2 if recon_r2 is not None else pd.Series(dtype=float)
            deconv = BulkDeconvResult(
                proportions=bulk_df,
                coverage_r2=coverage_r2.reindex(bulk_df.index, fill_value=float("nan")),
                gene_panel=list(bulk_df.columns),
                gene_weights=None,
                lower_ci=None,
                upper_ci=None,
                cell_fractions=None,
                run_metadata={"generated_from": "saved TSV outputs"},
            )
            bulk_result = SimpleNamespace(
                deconv=deconv,
                qc=qc,
                protocol_risk=None,
                run_metadata={"generated_from": "07_generate_reports.py"},
            )
            methods.append(methods_text.compose_bulk_methods(bulk_result))
    except Exception:
        pass

    try:
        spatial_df = _read_tsv(H.OUT_SPATIAL_DIR / "spatial_spot_proportions.tsv")
        spatial_qc = _read_tsv(H.OUT_SPATIAL_DIR / "spatial_qc.tsv")
        morans_df = _read_tsv(H.OUT_SPATIAL_DIR / "morans_i.tsv")
        morans_i = morans_df["morans_i"] if morans_df is not None else None
        spatial_run_meta = {}
        spatial_meta_path = H.OUT_SPATIAL_DIR / "spatial_run_metadata.json"
        if spatial_meta_path.exists():
            spatial_run_meta = json.loads(spatial_meta_path.read_text())
        if spatial_df is not None:
            qc = QCReport(
                modality="spatial",
                recommendations=[],
                metadata={},
                spot_qc=spatial_qc,
                morans_i=morans_i,
            )
            deconv = SpatialDeconvResult(
                proportions=spatial_df,
                cell_types=list(spatial_df.columns),
                marker_genes=[],
                n_iter=int(spatial_run_meta.get("n_iter", 0)),
                converged=bool(spatial_run_meta.get("converged", False)),
                convergence_trace=[],
                lambda_spatial=float(spatial_run_meta.get("lambda_spatial", 0.0)),
                mismatch_factors=None,
                lower_ci=None,
                upper_ci=None,
                bootstrap_coverage_note=None,
                n_smooth=spatial_run_meta.get("n_smooth"),
                run_metadata=spatial_run_meta,
            )
            spatial_result = SimpleNamespace(
                deconv=deconv,
                qc=qc,
                morans_i=morans_i,
                spot_qc=spatial_qc,
                run_metadata=spatial_run_meta,
            )
            methods.append(methods_text.compose_spatial_methods(spatial_result))
    except Exception:
        pass

    if methods:
        (H.OUT_SUMMARY_DIR / "methods.txt").write_text(
            "\n\n".join(methods), encoding="utf-8"
        )


def _separability_spillover_figs(ref, figdir: Path, prefix: str, warns: list):
    """Separability + spillover heatmaps from the reference (expression proxy)."""
    if ref is None:
        return
    try:
        from tissueresolve.benchmark.spillover import expression_spillover_proxy
        from tissueresolve.plotting.separability_plots import plot_separability_heatmap
        from tissueresolve.plotting.spillover_plots import (
            plot_spillover_heatmap, plot_spillover_network)
        from tissueresolve.reference.separability import (
            compute_separability, separability_heatmap_data)

        with _warnings.catch_warnings():
            _warnings.simplefilter("ignore")
            sep = compute_separability(ref, raise_on_critical=False)
        M = separability_heatmap_data(sep, list(ref.cell_types))
        plot_separability_heatmap(M, figdir, cell_types=list(ref.cell_types),
                                  name=f"{prefix}_separability_heatmap")
        spill = expression_spillover_proxy(ref)
        plot_spillover_heatmap(spill, figdir, name=f"{prefix}_spillover_heatmap")
        plot_spillover_network(spill, figdir, name="spillover_network")
    except Exception as exc:  # pragma: no cover - defensive
        warns.append(f"separability/spillover figures skipped: {exc}")


def generate_bulk_outputs(ref) -> tuple[Path, list[str], dict[str, object]]:
    from tissueresolve.plotting import bulk_plots
    from tissueresolve.report import generate_report

    rdir = H.OUT_BULK_DIR
    figdir = rdir / "figures"
    figdir.mkdir(parents=True, exist_ok=True)
    warns: list[str] = []

    props = _read_tsv(rdir / "bulk_estimated_proportions.tsv")
    qc = _read_tsv(rdir / "bulk_qc.tsv")
    if props is not None:
        try:
            bulk_plots.plot_bulk_composition_clustered_barplot(props, figdir)
            bulk_plots.plot_bulk_composition_heatmap(props, figdir)
            bulk_plots.plot_bulk_uncertainty(props, None, None, figdir)
        except Exception as exc:
            warns.append(f"bulk composition figures skipped: {exc}")
    if qc is not None:
        try:
            bulk_plots.plot_bulk_qc_summary(qc, figdir)
        except Exception as exc:
            warns.append(f"bulk QC figure skipped: {exc}")
    _separability_spillover_figs(ref, figdir, "bulk", warns)

    _copy_into_tables(rdir, [
        H.OUT_REFERENCE_DIR / "reference_summary.tsv",
        H.OUT_REFERENCE_DIR / "cell_type_counts.tsv",
        H.OUT_REFERENCE_DIR / "selected_gene_identifier_column.txt",
        H.OUT_RESOLUTION_DIR / "pairwise_resolvability.tsv",
        H.OUT_RESOLUTION_DIR / "spillover_risk_by_celltype.tsv",
        H.OUT_RESOLUTION_DIR / "recommended_merges.tsv",
        H.OUT_RESOLUTION_DIR / "bulk_family_proportions.tsv",
    ])
    bulk_meta = {"modality": "bulk", "n_samples":
                 int(props.shape[0]) if props is not None else None}
    H.write_json(bulk_meta, rdir / "run_metadata.json")
    out = generate_report("bulk", rdir, rdir / "report.html", warnings=warns)
    return out, warns, bulk_meta


def generate_spatial_outputs(ref) -> tuple[Path, list[str], dict[str, object]]:
    from tissueresolve.plotting import spatial_plots
    from tissueresolve.report import generate_report

    rdir = H.OUT_SPATIAL_DIR
    figdir = rdir / "figures"
    figdir.mkdir(parents=True, exist_ok=True)
    warns: list[str] = []

    props = _read_tsv(rdir / "spatial_spot_proportions.tsv")
    morans = _read_tsv(rdir / "morans_i.tsv")
    if props is not None:
        try:
            spatial_plots.plot_spatial_mean_composition_barplot(props, figdir)
        except Exception as exc:
            warns.append(f"spatial mean composition figure skipped: {exc}")
    if morans is not None:
        try:
            spatial_plots.plot_spatial_morans_i_barplot(morans["morans_i"], figdir)
        except Exception as exc:
            warns.append(f"Moran's I figure skipped: {exc}")

    # Coordinate-dependent maps need array coords from the Visium h5ad.
    coords = _spatial_coords(props)
    if coords is not None and props is not None:
        ar, ac = coords
        for fn in (lambda: spatial_plots.plot_spatial_abundance_maps(props, ar, ac, figdir),
                   lambda: spatial_plots.plot_spatial_dominant_cell_type_map(props, ar, ac, figdir),
                   lambda: spatial_plots.plot_spatial_spot_pie_charts(props, ar, ac, figdir)):
            try:
                fn()
            except Exception as exc:
                warns.append(f"spatial map skipped: {exc}")
    else:
        warns.append("spatial array coordinates unavailable; coordinate maps "
                     "(abundance/dominant/pie) were skipped.")

    _separability_spillover_figs(ref, figdir, "spatial", warns)
    _copy_into_tables(rdir, [
        H.OUT_REFERENCE_DIR / "reference_summary.tsv",
        H.OUT_REFERENCE_DIR / "cell_type_counts.tsv",
        H.OUT_REFERENCE_DIR / "selected_gene_identifier_column.txt",
        rdir / "spatial_run_metadata.json",
        H.OUT_RESOLUTION_DIR / "pairwise_resolvability.tsv",
        H.OUT_RESOLUTION_DIR / "spillover_risk_by_celltype.tsv",
        H.OUT_RESOLUTION_DIR / "recommended_merges.tsv",
        H.OUT_RESOLUTION_DIR / "spatial_family_proportions.tsv",
    ])
    meta: dict[str, object] = {}
    mp = rdir / "spatial_run_metadata.json"
    if mp.exists():
        try:
            meta = json.loads(mp.read_text())
        except Exception:
            meta = {}
    H.write_json(meta, rdir / "run_metadata.json")
    out = generate_report("spatial", rdir, rdir / "report.html",
                           run_metadata=meta, warnings=warns)
    return out, warns, meta


def _spatial_coords(props):
    """Try to read array_row/array_col from the local Visium h5ad (no download)."""
    if props is None or not H.SPATIAL_H5AD.exists():
        return None
    try:
        import anndata as ad

        adata = ad.read_h5ad(H.SPATIAL_H5AD)
        if "array_row" in adata.obs and "array_col" in adata.obs:
            sub = adata.obs.loc[[i for i in props.index if i in adata.obs.index]]
            if len(sub) == len(props):
                return (sub["array_row"].to_numpy(), sub["array_col"].to_numpy())
    except Exception:
        return None
    return None


def generate_combined_report() -> Path:
    from tissueresolve.report import templates as T

    rdir = H.OUT_SUMMARY_DIR
    rdir.mkdir(parents=True, exist_ok=True)
    summary = ""
    mp = rdir / "validation_summary.md"
    if mp.exists():
        summary = "<pre>" + T.escape(mp.read_text()) + "</pre>"
    body = T.section(0, "Validation summary", summary or "<p>not available</p>")
    body += T.section(1, "Reports", T.file_list([
        "../bulk/report.html", "../spatial/report.html"]))
    body += T.section(2, "Estimate types", T.estimate_box(
        "Bulk: mRNA-derived proportions (not cell fractions). "
        "Spatial: spot-level RNA-derived composition (not cell counts)."))
    out = rdir / "report.html"
    out.write_text(T.page("TissueResolve — Breast cancer validation report", body),
                   encoding="utf-8")
    return out


def main(argv: list[str] | None = None) -> int:
    argparse.ArgumentParser(description=__doc__).parse_args(argv)
    H.ensure_dirs()
    ref = _load_reference()
    generated = []
    bulk_warns: list[str] = []
    spatial_warns: list[str] = []
    bulk_meta: dict[str, object] = {}
    spatial_meta: dict[str, object] = {}

    if H.OUT_BULK_DIR.exists():
        bulk_path, bulk_warns, bulk_meta = generate_bulk_outputs(ref)
        generated.append(bulk_path)
    if H.OUT_SPATIAL_DIR.exists():
        spatial_path, spatial_warns, spatial_meta = generate_spatial_outputs(ref)
        generated.append(spatial_path)
    generated.append(generate_combined_report())
    _write_report_bundle_artifacts(bulk_warns, spatial_warns, bulk_meta, spatial_meta)

    print("Generated reports:")
    for p in generated:
        print(f"  {p}")
    print(f"Bundle metadata: {H.OUT_SUMMARY_DIR / 'warnings.json'}, {H.OUT_SUMMARY_DIR / 'run_metadata.json'}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
