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


def generate_unified_report() -> Path:
    """One report.html with nav linking every section (the main entry point)."""
    from tissueresolve.report import templates as T
    from tissueresolve.report.unified import Section, build_unified_report

    out_dir = H.OUTPUTS_DIR
    rel = lambda p: __import__("os").path.relpath(p, out_dir)  # noqa: E731

    sections = []

    md = H.OUT_SUMMARY_DIR / "validation_summary.md"
    sections.append(Section("summary", "1. Executive summary",
                            "<pre>" + T.escape(_read_text(md)) + "</pre>"
                            if md.exists() else ""))

    ref_sum = H.OUT_REFERENCE_DIR / "reference_summary.tsv"
    refq = ""
    if ref_sum.exists():
        try:
            refq = pd.read_csv(ref_sum, sep="\t", comment="#").head(40).to_html(index=False, border=0)
        except Exception:
            refq = ""
    sections.append(Section("reference", "2. Reference quality", refq,
                            links=[("reference tables", rel(H.OUT_REFERENCE_DIR))]))

    sections.append(Section("input", "3. Input data",
                            "<p>Bulk pseudobulk mixtures and a 10x Visium section "
                            "(see the linked detailed reports).</p>"))

    bulk_links = [("detailed bulk report", rel(H.OUT_BULK_DIR / "report.html"))] \
        if (H.OUT_BULK_DIR / "report.html").exists() else []
    sections.append(Section("bulk", "4. Bulk results",
                            "<p>RNA-derived mRNA proportions (not cell fractions).</p>",
                            links=bulk_links))

    spatial_links = [("detailed spatial report", rel(H.OUT_SPATIAL_DIR / "report.html"))] \
        if (H.OUT_SPATIAL_DIR / "report.html").exists() else []
    sections.append(Section("spatial", "5. Spatial results",
                            "<p>Spot-level RNA-derived composition (not cell counts).</p>",
                            links=spatial_links))

    hier_dir = H.OUTPUTS_DIR / "hierarchical"
    hmd = hier_dir / "hierarchical_summary.md"
    hbody = "<pre>" + T.escape(_read_text(hmd)) + "</pre>" if hmd.exists() else ""
    hlinks = [("hierarchical tables", rel(hier_dir))] if hier_dir.exists() else []
    sections.append(Section("hierarchical", "6. Hierarchical broad→fine results",
                            hbody, links=hlinks))

    res_dir = H.OUT_RESOLUTION_DIR
    rlinks = [("resolution / separability / spillover", rel(res_dir))] if res_dir.exists() else []
    sections.append(Section("resolution", "7. Resolution, separability & spillover",
                            "<p>Pairwise separability and spillover diagnostics.</p>",
                            links=rlinks))

    bench = H.HARNESS_DIR.parent.parent / "benchmarks" / "outputs" / "benchmark_summary_report.html"
    blinks = [("benchmark summary report", rel(bench))] if bench.exists() else []
    sections.append(Section("benchmark", "8. Benchmark comparison",
                            "<p>Comparison of TissueResolve (flat + hierarchical) "
                            "against baselines and external tools. Accuracy is "
                            "reported only where ground truth exists.</p>",
                            links=blinks))

    warns = H.OUT_SUMMARY_DIR / "warnings.json"
    wbody = ""
    if warns.exists():
        try:
            w = json.loads(_read_text(warns))
            items = w if isinstance(w, list) else w.get("warnings", [])
            wbody = ("<div class='warn'><ul>" + "".join(
                f"<li>{T.escape(str(x))}</li>" for x in items[:50]) + "</ul></div>") \
                if items else "<p>No warnings recorded.</p>"
        except Exception:
            wbody = ""
    sections.append(Section("warnings", "9. Warnings & recommended actions", wbody))

    figs = []
    for d in (H.OUT_BULK_DIR, H.OUT_SPATIAL_DIR):
        figs += sorted((d / "figures").glob("*.html")) if (d / "figures").exists() else []
    fbody = ("<ul>" + "".join(f"<li><a href='{rel(f)}'>{f.name}</a></li>"
                              for f in figs[:60]) + "</ul>") if figs else ""
    sections.append(Section("figures", "10. Publication figures", fbody))

    methods = H.OUT_SUMMARY_DIR / "methods.txt"
    sections.append(Section("methods", "11. Methods",
                            "<pre>" + T.escape(_read_text(methods)) + "</pre>"
                            if methods.exists() else ""))

    out_files = sorted(p for p in out_dir.rglob("*.tsv"))[:80]
    obody = ("<ul>" + "".join(f"<li>{rel(p)}</li>" for p in out_files) + "</ul>") \
        if out_files else ""
    sections.append(Section("outputs", "12. Output files", obody))

    path = build_unified_report(
        out_dir / "report.html", sections,
        title="TissueResolve — breast-cancer analysis report",
        subtitle="One page; click a section above. Detailed sub-reports are linked.")
    return path


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
