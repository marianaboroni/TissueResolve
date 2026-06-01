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


def _reference_family_counts(ref):
    from tissueresolve.plotting.palette import infer_cell_type_family

    if ref is None or not getattr(ref, "n_cells_per_type", None):
        return None
    fam: dict[str, int] = {}
    for ct, n in ref.n_cells_per_type.items():
        f = infer_cell_type_family(ct)
        fam[f] = fam.get(f, 0) + int(n)
    return pd.Series(fam) if fam else None


def _read_pairs(path: Path):
    try:
        df = pd.read_csv(path, sep="\t", comment="#")
        if {"type_a", "type_b"}.issubset(df.columns):
            return df
    except Exception:
        pass
    return None


def _color_map_and_save(props, figdir: Path):
    from tissueresolve.plotting.palette import assign_family_palette, save_color_map

    cols = list(props.columns) + ["Other"] if props is not None else ["Other"]
    cmap = assign_family_palette(cols)
    save_color_map(cmap, figdir)
    return cmap


def generate_bulk_outputs(ref) -> tuple[Path, list[str], dict[str, object]]:
    from tissueresolve.plotting import bulk_plots
    from tissueresolve.plotting.summary_figures import bulk_main_summary_figure
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
    # Main publication summary figure + stable family palette.
    if props is not None:
        try:
            cmap = _color_map_and_save(props, figdir)
            bulk_main_summary_figure(
                props, figdir, reference_family_counts=_reference_family_counts(ref),
                top_pairs=_read_pairs(H.OUT_RESOLUTION_DIR / "pairwise_resolvability.tsv"),
                color_map=cmap)
        except Exception as exc:
            warns.append(f"bulk main summary figure skipped: {exc}")
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
    from tissueresolve.plotting.export import kaleido_available
    bulk_meta = {"modality": "bulk",
                 "n_samples": int(props.shape[0]) if props is not None else None,
                 "static_export": bool(kaleido_available())}
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

    # Main publication summary figure + family palette.
    has_he = False
    if props is not None:
        try:
            cmap = _color_map_and_save(props, figdir)
            from tissueresolve.plotting.summary_figures import spatial_main_summary_figure
            coords_df = (pd.DataFrame({"array_row": coords[0], "array_col": coords[1]},
                                      index=props.index) if coords is not None else None)
            spatial_main_summary_figure(
                props, figdir, coords=coords_df,
                morans=morans["morans_i"] if morans is not None else None,
                reference_family_counts=_reference_family_counts(ref), color_map=cmap)
        except Exception as exc:
            warns.append(f"spatial main summary figure skipped: {exc}")
        # H&E overlays (from the local Visium .h5ad if present).
        try:
            has_he = _he_overlays(props, figdir, cmap, warns)
        except Exception as exc:
            warns.append(f"H&E overlays skipped: {exc}")

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
    from tissueresolve.plotting.export import kaleido_available
    meta["has_he_image"] = bool(has_he)
    meta["static_export"] = bool(kaleido_available())
    H.write_json(meta, rdir / "run_metadata.json")
    out = generate_report("spatial", rdir, rdir / "report.html",
                           run_metadata=meta, warnings=warns)
    return out, warns, meta


def _he_overlays(props, figdir: Path, cmap: dict, warns: list) -> bool:
    """Generate H&E overlays from the local Visium .h5ad; returns has_image."""
    if not H.SPATIAL_H5AD.exists():
        warns.append("Visium .h5ad not present; H&E overlays skipped "
                     "(coordinate maps used instead).")
        return False
    import anndata as ad

    from tissueresolve.plotting.histology import (
        extract_visium_coordinates, load_visium_histology_image,
        plot_abundance_on_he, plot_dominant_cell_type_on_he, plot_he_with_spots)

    adata = ad.read_h5ad(H.SPATIAL_H5AD)
    image, sf, w = load_visium_histology_image(adata)
    warns.extend(w)
    coords = extract_visium_coordinates(adata, scalefactors=sf)
    coords = coords.set_index("spot")
    common = [s for s in props.index if s in coords.index]
    if not common:
        warns.append("No spot-id overlap between Visium image coords and "
                     "predictions; H&E overlays skipped.")
        return False
    coords_he = coords.loc[common].reset_index()
    props_he = props.loc[common]
    plot_he_with_spots(image, coords_he, figdir, name="he_spots_check")
    plot_dominant_cell_type_on_he(image, coords_he, props_he, cmap, figdir,
                                  name="he_dominant_cell_type")
    for ct in props_he.mean(0).sort_values(ascending=False).index[:3]:
        plot_abundance_on_he(image, coords_he, props_he[ct], figdir,
                             cell_type=str(ct))
    return image is not None


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


def _hierarchical_report_html() -> str:
    """Return an HTML block summarising hierarchical outputs, if present."""
    from tissueresolve.report import templates as T

    hdir = H.OUTPUTS_DIR / "hierarchical"
    if not hdir.exists():
        return ""
    blocks = []
    md = hdir / "hierarchical_summary.md"
    if md.exists():
        blocks.append("<pre>" + T.escape(md.read_text()) + "</pre>")
    qc = hdir / "hierarchical_qc.tsv"
    if qc.exists():
        try:
            df = pd.read_csv(qc, sep="\t", comment="#")
            blocks.append("<h3>Within-family resolvability</h3>"
                          + df.to_html(index=False, border=0))
        except Exception:
            pass
    tsvs = sorted(p.name for p in hdir.glob("*.tsv"))
    if tsvs:
        blocks.append("<h3>Tables</h3>"
                      + T.file_list([f"../hierarchical/{n}" for n in tsvs]))
    blocks.append(T.estimate_box(
        "Hierarchical mode first estimates broad cell-type families, then fine "
        "subpopulations within each family.  Families whose subtypes are not "
        "separable are reported at the broad level as "
        "<code>unresolved_&lt;family&gt;</code> — subtype splits there are not "
        "claimed.  Colours are consistent across figures; see "
        "<code>../hierarchical/cell_type_color_map.tsv</code>."))
    return "".join(blocks)


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
    idx = 2
    hier_html = _hierarchical_report_html()
    if hier_html:
        body += T.section(idx, "Broad-to-fine hierarchical deconvolution", hier_html)
        idx += 1
    body += T.section(idx, "Estimate types", T.estimate_box(
        "Bulk: mRNA-derived proportions (not cell fractions). "
        "Spatial: spot-level RNA-derived composition (not cell counts)."))
    out = rdir / "report.html"
    out.write_text(T.page("TissueResolve — Breast cancer validation report", body),
                   encoding="utf-8")
    return out


def _read_text(path) -> str:
    from pathlib import Path as _P
    p = _P(path)
    return p.read_text(encoding="utf-8") if p.exists() else ""


def _df_collapsible(title, path, max_rows=12, comment="#"):
    """Top rows inline + full table as a collapsible + source-data link."""
    from tissueresolve.report import components as C
    from pathlib import Path as _P
    p = _P(path)
    if not p.exists():
        return ""
    try:
        df = pd.read_csv(p, sep="\t", comment=comment)
    except Exception:
        return ""
    head = df.head(max_rows).to_html(index=False, border=0)
    note = f"(showing top {min(max_rows, len(df))} of {len(df)} rows)" if len(df) > max_rows else ""
    return C.collapsible_table(title, head, note=note)


def _figure_cards(fig_dir, out_dir, manifest, section, captions):
    """Build figure cards for every figure HTML in *fig_dir*, register in manifest."""
    from tissueresolve.report import components as C
    from tissueresolve.report.figures import FigureRecord
    import os
    from pathlib import Path as _P
    fig_dir = _P(fig_dir)
    if not fig_dir.exists():
        return ""
    rel = lambda p: os.path.relpath(p, out_dir)  # noqa: E731
    cards = []
    for fh in sorted(fig_dir.glob("*.html")):
        stem = fh.stem
        cap = captions.get(stem, captions.get("_default"))
        data = fh.with_suffix(".data.tsv")
        links = [("interactive figure", rel(fh))]
        if data.exists():
            links.append(("source data (.tsv)", rel(data)))
        for ext in ("png", "svg", "pdf"):
            ex = fh.with_suffix("." + ext)
            if ex.exists():
                links.append((f"{ext.upper()}", rel(ex)))
        title = stem.replace("_", " ").strip().capitalize()
        cards.append(C.figure_card(
            title=title, subtitle=cap["subtitle"], body_html="",
            caption=cap["caption"], how_to_read=cap["how_to_read"],
            methodology=cap.get("methodology", ""),
            source_links=links))
        manifest.add(FigureRecord(
            figure_id=stem, section=section, title=title,
            html_path=rel(fh), source_data_path=rel(data) if data.exists() else "",
            methodology=cap.get("methodology", ""), status="ok"))
    return "".join(cards)


def generate_unified_report() -> Path:
    """One report.html with sidebar nav, metric cards, methodology boxes,
    glossary, figure cards and collapsible tables (the main entry point)."""
    import os, warnings as _w
    from tissueresolve.report import components as C
    from tissueresolve.report import glossary as G
    from tissueresolve.report.figures import FigureManifest
    from tissueresolve.report.unified import Section, build_unified_report

    out_dir = H.OUTPUTS_DIR
    rel = lambda p: os.path.relpath(p, out_dir)  # noqa: E731
    manifest = FigureManifest()
    sections = []

    # ---- gather context (defensive) ----
    ctx = {"modality": "bulk + spatial", "samples": "—", "spots": "—",
           "ref_cells": "—", "genes": "—", "broad": "—", "fine": "—",
           "overlap": "—", "solver": "auto", "resolution_mode": "hierarchical",
           "high_risk_pairs": "—", "suit": None}
    ref0 = None
    mp = None
    try:
        from tissueresolve.results import ReferenceSignature
        from tissueresolve.reference.hierarchy import (
            load_hierarchy_mapping, build_cell_type_hierarchy)
        if H.SAVED_REFERENCE_DIR.exists():
            ref0 = ReferenceSignature.load(H.SAVED_REFERENCE_DIR)
            ctx["fine"] = ref0.n_cell_types
            ctx["genes"] = ref0.n_genes
            ctx["ref_cells"] = sum(ref0.n_cells_per_type.values()) or "—"
        hmap_p = H.HARNESS_DIR / "config" / "breast_cancer_cell_type_hierarchy.tsv"
        if ref0 is not None and hmap_p.exists():
            mp = build_cell_type_hierarchy(list(ref0.cell_types), load_hierarchy_mapping(hmap_p))
            ctx["broad"] = len(set(mp.values()))
    except Exception:
        pass
    try:
        if H.PSEUDOBULK_TRUE_PROPS.exists():
            tp = pd.read_csv(H.PSEUDOBULK_TRUE_PROPS, sep="\t", comment="#", index_col=0)
            ctx["samples"] = tp.shape[0]
    except Exception:
        pass
    try:
        sp = H.OUT_SPATIAL_DIR / "spatial_spot_proportions.tsv"
        if sp.exists():
            ctx["spots"] = pd.read_csv(sp, sep="\t", comment="#", index_col=0).shape[0]
    except Exception:
        pass
    try:
        if ref0 is not None and H.PSEUDOBULK_COUNTS.exists():
            qg = set(pd.read_csv(H.PSEUDOBULK_COUNTS, sep="\t", comment="#", index_col=0).index.map(str))
            ctx["overlap"] = len(qg & set(map(str, ref0.gene_names)))
    except Exception:
        pass
    bench = H.HARNESS_DIR.parent.parent / "benchmarks" / "outputs" / "benchmark_summary_report.html"
    bench_status = "available" if bench.exists() else "not run"

    # ---- 1. Executive summary ----
    cards = {"Modality": ctx["modality"], "Samples (bulk)": ctx["samples"],
             "Spots (spatial)": ctx["spots"], "Reference cells": ctx["ref_cells"],
             "Genes": ctx["genes"], "Broad families": ctx["broad"],
             "Fine subpopulations": ctx["fine"], "Gene overlap": ctx["overlap"],
             "Solver": ctx["solver"], "Resolution mode": ctx["resolution_mode"],
             "Benchmark": bench_status, "Report": "generated"}
    summ = (C.metric_grid(cards)
            + C.estimate_note(
                "This report summarizes reference quality, input-data checks, "
                "deconvolution predictions, model diagnostics, resolution limits, "
                "and benchmark results. The values shown are <b>RNA-derived "
                "estimates</b> and should be interpreted together with QC, "
                "separability, and uncertainty metrics. The report does not infer "
                "biological meaning beyond what the data support.")
            + C.interpretation_guide(
                "Use the sidebar to jump to a section. Start here for the headline "
                "numbers, then check reference quality and warnings before trusting "
                "fine-grained predictions."))
    md = H.OUT_SUMMARY_DIR / "validation_summary.md"
    if md.exists():
        summ += C.collapsible_table("Text validation summary",
                                    "<pre>" + C.esc(_read_text(md)) + "</pre>")
    sections.append(Section("summary", "1. Executive summary", summ))

    # ---- 2. Reference quality ----
    refq = C.methodology_summary([
        "Reference built from a single-cell/nucleus h5ad by aggregating per-cell "
        "profiles into per-cell-type signatures (CPM).",
        "Broad and fine labels taken from the documented hierarchy mapping; gene "
        "identifiers harmonized to symbols.",
        "Minimum cells per type enforced; cross-donor variability recorded when "
        "multiple donors are present.",
    ])
    try:
        from tissueresolve.reference.suitability import (
            compute_reference_suitability_score, save_reference_suitability)
        if ref0 is not None:
            qgenes = (list(pd.read_csv(H.PSEUDOBULK_COUNTS, sep="\t", index_col=0,
                                       comment="#").index.map(str))
                      if H.PSEUDOBULK_COUNTS.exists() else None)
            with _w.catch_warnings():
                _w.simplefilter("ignore")
                suit = compute_reference_suitability_score(ref0, query_genes=qgenes, mapping=mp)
                save_reference_suitability(suit, H.OUT_REFERENCE_DIR)
            refq += (f"<p>Overall suitability: {C.status_badge(suit.classification)} "
                     f"(score {suit.overall_score})</p>"
                     + C.collapsible_table("Suitability components",
                                           suit.components_frame().to_html(border=0)))
    except Exception:
        pass
    refq += _figure_cards(H.OUT_REFERENCE_DIR / "figures", out_dir, manifest,
                          "reference", _CAPTIONS)
    refq += C.variable_dictionary(G.subset(["separability", "spillover", "gene overlap"]))
    refq += C.interpretation_guide(
        "Use this section to decide whether your reference is balanced, compatible "
        "with your query, and sufficiently annotated. CAUTION/WARNING/FAIL flags "
        "highlight components to check before trusting fine predictions.")
    sections.append(Section("reference", "2. Reference quality", refq,
                            links=[("reference tables", rel(H.OUT_REFERENCE_DIR))]))

    # ---- 3. Input data ----
    inp = C.methodology_summary([
        "Bulk: pseudobulk count mixtures with known ground-truth proportions.",
        "Spatial: a 10x Visium section (counts + array coordinates; H&E when present).",
        "Gene overlap with the reference is computed and reported; nothing is "
        "silently re-normalized.",
    ]) + C.metric_grid({"Bulk samples": ctx["samples"], "Spatial spots": ctx["spots"],
                        "Gene overlap": ctx["overlap"]})
    inp += C.interpretation_guide(
        "Check that gene overlap is high and that the input type/normalization "
        "matches expectations before interpreting predictions.")
    sections.append(Section("input", "3. Input data quality", inp))

    # ---- 4. Bulk ----
    bulk = C.methodology_summary([
        "Solver backbone selected automatically (solver=auto) by gene-masking "
        "cross-validation; the chosen backbone is recorded.",
        "Estimates are mRNA-derived proportions (rows sum to 1), NOT absolute "
        "cell fractions.",
    ])
    bulk += _figure_cards(H.OUT_BULK_DIR / "figures", out_dir, manifest, "bulk", _CAPTIONS)
    bulk += _df_collapsible("Estimated proportions (source data)",
                            H.OUT_BULK_DIR / "bulk_estimated_proportions.tsv")
    bulk += C.variable_dictionary(G.subset(["mRNA-derived proportion", "coverage R²", "unresolved mass"]))
    bulk += C.interpretation_guide(
        "Use this section to compare composition across samples. Check whether "
        "clustering matches expected sample groups and whether high-risk or "
        "high-spillover populations dominate a sample.")
    blinks = [("detailed bulk report", rel(H.OUT_BULK_DIR / "report.html"))] \
        if (H.OUT_BULK_DIR / "report.html").exists() else []
    sections.append(Section("bulk", "4. Bulk deconvolution", bulk, links=blinks))

    # ---- 5. Spatial ----
    he_present = any((H.OUT_SPATIAL_DIR / "figures").glob("*he*")) if (H.OUT_SPATIAL_DIR / "figures").exists() else False
    spat = C.methodology_summary([
        "Counts deconvolved per spot with an NB-CAR model; spatial smoothing "
        "(lambda_spatial) is recorded and never hidden.",
        "Array coordinates define the spot graph; H&E overlays are shown when an "
        "image is available.",
        "Real Visium has NO ground truth — results are concordance/structure, not "
        "accuracy.",
    ])
    spat += _figure_cards(H.OUT_SPATIAL_DIR / "figures", out_dir, manifest, "spatial", _CAPTIONS)
    spat += C.variable_dictionary(G.subset(["entropy", "dominant fraction", "near-zero fraction", "Moran's I"]))
    spat += C.interpretation_guide(
        "Use this section to inspect where predicted RNA-derived compositions "
        "localize in tissue. Compare overlays with H&E morphology, but do not "
        "treat predictions as direct cell counts.")
    if he_present:
        spat = "<p class='muted'>H&E overlay figures are included below.</p>" + spat
    slinks = [("detailed spatial report", rel(H.OUT_SPATIAL_DIR / "report.html"))] \
        if (H.OUT_SPATIAL_DIR / "report.html").exists() else []
    sections.append(Section("spatial", "5. Spatial deconvolution", spat, links=slinks))

    # ---- 6. Hierarchical ----
    hier_dir = H.OUTPUTS_DIR / "hierarchical"
    hier = C.methodology_summary([
        "Broad-to-fine strategy: estimate broad families first, then fine "
        "subpopulations within each family.",
        "Partial resolution: confident subtypes receive mass; ambiguous remainder "
        "is reported as unresolved_<family>.",
        "Fine subtype estimates in low-separability families should be treated "
        "cautiously (interpret at the family level).",
    ])
    hmd = hier_dir / "hierarchical_summary.md"
    if hmd.exists():
        hier += C.collapsible_table("Hierarchical summary",
                                    "<pre>" + C.esc(_read_text(hmd)) + "</pre>", open=True)
    hier += _df_collapsible("Within-family resolvability QC",
                            hier_dir / "hierarchical_qc.tsv")
    hier += C.variable_dictionary(G.subset(["unresolved mass", "separability", "spillover"]))
    hier += C.interpretation_guide(
        "Use this section to see which families could be split into reliable "
        "subtypes and which are reported at the family level (unresolved mass).")
    hlinks = [("hierarchical tables", rel(hier_dir))] if hier_dir.exists() else []
    sections.append(Section("hierarchical", "6. Hierarchical broad→fine deconvolution",
                            hier, links=hlinks))

    # ---- 7. Resolution / separability / spillover ----
    res_dir = H.OUT_RESOLUTION_DIR
    resb = C.methodology_summary([
        "Pairwise separability computed from the reference expression profiles "
        "(1 − Bhattacharyya coefficient).",
        "High-risk pairs (BC > 0.90) and per-family resolution summaries highlight "
        "where fine labels are unreliable.",
    ])
    resb += _df_collapsible("Top non-separable pairs",
                            res_dir / "pairwise_separability.tsv")
    resb += C.variable_dictionary(G.subset(["separability", "spillover", "unresolved mass"]))
    resb += C.interpretation_guide(
        "Use this to judge what resolution you can trust: families with many "
        "non-separable subtypes should be interpreted at the broad level.")
    rlinks = [("resolution outputs", rel(res_dir))] if res_dir.exists() else []
    sections.append(Section("resolution", "7. Resolution, separability & spillover",
                            resb, links=rlinks))

    # ---- 8. Benchmark ----
    benchb = C.methodology_summary([
        "Methods compared on identical harmonized inputs.",
        "Only EXECUTED or IMPORTED methods are ranked; EXPORTED-only and SKIPPED "
        "tools are listed but not scored.",
        "Bulk pseudobulk has ground truth (accuracy valid); real Visium does not "
        "(concordance/structure only).",
    ])
    benchb += _df_collapsible("Best-method summary",
                              bench.parent / "bulk" / "best_method_summary.tsv")
    benchb += C.variable_dictionary(G.subset(
        ["Pearson correlation", "RMSE", "concordance", "composite score"]))
    benchb += C.interpretation_guide(
        "Use this section to compare methods. Only executed or imported methods "
        "are ranked; exported-only tools are listed but not scored.")
    blinks2 = [("benchmark summary report", rel(bench))] if bench.exists() else []
    sections.append(Section("benchmark", "8. Benchmark comparison", benchb, links=blinks2))

    # ---- 9. Warnings & limitations ----
    warns = H.OUT_SUMMARY_DIR / "warnings.json"
    warn_items = []
    if warns.exists():
        try:
            w = json.loads(_read_text(warns))
            warn_items = w if isinstance(w, list) else w.get("warnings", [])
        except Exception:
            warn_items = []
    wbody = C.warning_box(warn_items) if warn_items else \
        "<p class='muted'>No warnings recorded for this run.</p>"
    wbody += C.limitation_box([
        "Bulk estimates are RNA-derived proportions, not absolute cell fractions.",
        "Spatial estimates are spot-level RNA-derived composition, not cell counts; "
        "real Visium has no ground truth.",
        "Fine subtype estimates in non-separable families are reported as "
        "unresolved mass; do not over-interpret them.",
    ])
    sections.append(Section("warnings", "9. Warnings & limitations", wbody))

    # ---- 10. Methods ----
    methods = H.OUT_SUMMARY_DIR / "methods.txt"
    mbody = ("<pre>" + C.esc(_read_text(methods)) + "</pre>" if methods.exists()
             else "<p class='muted'>Methods text not available.</p>")
    sections.append(Section("methods", "10. Methods", mbody))

    # ---- 11. Output files & source data ----
    man_path = manifest.save(out_dir / "figures")
    out_files = sorted(p for p in out_dir.rglob("*.tsv"))[:100]
    obody = (C.estimate_note(f"Figure manifest: {rel(man_path)} "
                             f"({len(manifest.records)} figures registered).")
             + C.collapsible_table("All source-data tables (.tsv)",
                                   "<ul>" + "".join(f"<li>{C.esc(rel(p))}</li>"
                                                    for p in out_files) + "</ul>"))
    sections.append(Section("outputs", "11. Output files & source data", obody))

    path = build_unified_report(
        out_dir / "report.html", sections,
        title="TissueResolve — breast-cancer analysis report",
        subtitle="RNA-derived estimates with QC, resolution limits and benchmarks")
    return path


# Per-figure captions: title is derived from filename; these add the
# subtitle / full caption / how-to-read.  '_default' covers any unlisted figure.
_CAPTIONS = {
    "_default": {
        "subtitle": "TissueResolve figure",
        "caption": "Visual summary of a TissueResolve output. Axes, colors and "
                   "values are described in the linked source data; values are "
                   "RNA-derived estimates, not absolute cell counts.",
        "how_to_read": "Open the interactive figure and its source data (.tsv). "
                       "Check whether patterns match expected groupings and whether "
                       "low-confidence or high-spillover populations dominate.",
        "methodology": "Generated by TissueResolve from the harmonized reference "
                       "and query; see the section methodology box for details.",
    },
}


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
    unified = generate_unified_report()
    generated.append(unified)
    _write_output_index()

    print("Generated reports:")
    for p in generated:
        print(f"  {p}")
    print(f"\n>>> MAIN REPORT: {unified}")
    print(f"Bundle metadata: {H.OUT_SUMMARY_DIR / 'warnings.json'}, {H.OUT_SUMMARY_DIR / 'run_metadata.json'}")
    return 0


def _write_output_index() -> None:
    """Write a short index.html + README.md at the outputs root ('where are my results')."""
    idx = H.OUTPUTS_DIR / "index.html"
    idx.write_text(
        "<!doctype html><meta charset='utf-8'><title>TissueResolve outputs</title>"
        "<h1>TissueResolve outputs</h1><ul>"
        "<li><a href='report.html'><b>report.html</b> — main unified report (start here)</a></li>"
        "<li><a href='validation_summary/report.html'>validation_summary/report.html</a></li>"
        "<li><a href='bulk/report.html'>bulk/report.html</a></li>"
        "<li><a href='spatial/report.html'>spatial/report.html</a></li>"
        "<li>hierarchical/ — broad→fine tables</li>"
        "<li>reference/, resolution/ — reference and diagnostics</li>"
        "</ul>", encoding="utf-8")
    (H.OUTPUTS_DIR / "README.md").write_text(
        "# TissueResolve outputs\n\n"
        "**Start here:** `report.html` (unified report with section navigation).\n\n"
        "- `report.html` — main report (reference, bulk, spatial, hierarchical, "
        "benchmark, warnings, figures, methods, outputs)\n"
        "- `bulk/`, `spatial/` — detailed per-modality reports + tables/figures\n"
        "- `hierarchical/` — broad→fine proportions, unresolved mass, QC\n"
        "- `reference/`, `resolution/` — reference summary and diagnostics\n"
        "- `validation_summary/` — warnings.json, run_metadata.json, methods.txt\n",
        encoding="utf-8")


if __name__ == "__main__":
    raise SystemExit(main())
