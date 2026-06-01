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


def _collect_warnings() -> list[str]:
    """Aggregate the run's *real* warnings from every diagnostic source.

    The previous report only read ``validation_summary/warnings.json`` and
    looked for a non-existent ``"warnings"`` key, so it always showed
    "No warnings recorded" even though the QC, spillover and separability
    outputs contain genuine cautions.  This gathers them honestly:

    * figure/processing warnings — ``validation_summary/warnings.json``
      (``{"bulk": [...], "spatial": [...]}``),
    * bulk QC recommendations — ``bulk/bulk_warnings.json`` (mismatch flags,
      high-spillover/collinear cell types),
    * spatial QC — ``spatial/spatial_warnings.json`` (non-convergence, low gene
      overlap, recommendations),
    * separability — ``resolution/pairwise_separability.tsv`` (poorly separable
      pairs), and
    * unresolved families — ``resolution/unresolved_families.tsv``.

    Returns a de-duplicated, order-preserving list of warning strings.
    """
    items: list[str] = []

    def _add(msg: str) -> None:
        msg = (msg or "").strip()
        if msg and msg not in items:
            items.append(msg)

    # 1. figure / processing warnings recorded during report generation
    vs = H.OUT_SUMMARY_DIR / "warnings.json"
    if vs.exists():
        try:
            w = json.loads(_read_text(vs))
            if isinstance(w, dict):
                for mod in ("bulk", "spatial"):
                    for m in w.get(mod, []) or []:
                        _add(f"{mod.capitalize()} processing: {m}")
            elif isinstance(w, list):
                for m in w:
                    _add(str(m))
        except Exception:
            pass

    # 2. bulk QC / spillover recommendations
    bw = H.OUT_BULK_DIR / "bulk_warnings.json"
    if bw.exists():
        try:
            w = json.loads(_read_text(bw))
            for rec in (w.get("qc_recommendations") or []):
                _add(f"Bulk QC: {rec}")
            go = w.get("gene_overlap") or {}
            if go.get("n_shared") is not None and go.get("n_reference"):
                frac = go["n_shared"] / max(go["n_reference"], 1)
                if frac < 0.5:
                    _add(f"Bulk gene overlap is low: {go['n_shared']} of "
                         f"{go['n_reference']} reference genes shared "
                         f"({frac:.0%}).")
        except Exception:
            pass

    # 3. spatial QC
    sw = H.OUT_SPATIAL_DIR / "spatial_warnings.json"
    if sw.exists():
        try:
            w = json.loads(_read_text(sw))
            if w.get("converged") is False:
                _add("Spatial model did not converge; spot estimates are "
                     "unreliable.")
            for rec in (w.get("qc_recommendations") or []):
                _add(f"Spatial QC: {rec}")
            go = w.get("gene_overlap") or {}
            if go.get("n_shared") is not None and go.get("n_reference"):
                frac = go["n_shared"] / max(go["n_reference"], 1)
                if frac < 0.5:
                    _add(f"Spatial gene overlap is low: {go['n_shared']} of "
                         f"{go['n_reference']} reference genes shared "
                         f"({frac:.0%}).")
        except Exception:
            pass

    # 4. separability — poorly separable cell-type pairs
    sep = H.OUTPUTS_DIR / "resolution" / "pairwise_separability.tsv"
    if sep.exists():
        try:
            df = pd.read_csv(sep, sep="\t", comment="#")
            if "resolvability" in df.columns:
                n_bad = int((df["resolvability"].astype(str)
                             .str.lower() == "unresolved").sum())
                n_tot = len(df)
                if n_bad:
                    _add(f"{n_bad} of {n_tot} cell-type pairs are not "
                         "separable from the reference; their fine subtype "
                         "estimates are unreliable and are reported as "
                         "unresolved family mass.")
        except Exception:
            pass

    # 5. unresolved families
    uf = H.OUTPUTS_DIR / "resolution" / "unresolved_families.tsv"
    if uf.exists():
        try:
            df = pd.read_csv(uf, sep="\t", comment="#")
            n = len(df)
            if n:
                _add(f"{n} cell-type family/families could not be resolved into "
                     "subtypes at this reference's resolution; subtype mass is "
                     "pooled as unresolved.")
        except Exception:
            pass

    return items


def _component_status(suit, name: str, default: str = "CAUTION") -> str:
    """Status (PASS/CAUTION/WARNING/FAIL/UNKNOWN) of one suitability component.

    ``components_frame`` indexes by component name, so reset the index to look it
    up robustly whether ``component`` is the index or a column.
    """
    try:
        df = suit.components_frame().reset_index()
        col = "component" if "component" in df.columns else df.columns[0]
        row = df[df[col].astype(str) == name]
        if not row.empty:
            return str(row.iloc[0]["status"])
    except Exception:
        pass
    return default


def _hierarchy_summary() -> dict:
    """Read hierarchical_summary.json (families, n_fine, unresolved families)."""
    p = H.OUTPUTS_DIR / "hierarchical" / "hierarchical_summary.json"
    if not p.exists():
        return {}
    try:
        return json.loads(_read_text(p))
    except Exception:
        return {}


def _decision_statuses(ctx: dict, suit, mp) -> dict:
    """Derive the executive-decision status cards from real diagnostics.

    Returns a dict with PASS/CAUTION/WARNING/FAIL (or UNKNOWN) per dimension plus
    availability flags.  Nothing biological is inferred — these are QC verdicts.
    """
    bulk_ok = (H.OUT_BULK_DIR / "bulk_estimated_proportions.tsv").exists()
    spatial_ok = (H.OUT_SPATIAL_DIR / "spatial_spot_proportions.tsv").exists()
    bench = H.HARNESS_DIR.parent.parent / "benchmarks" / "outputs" / \
        "benchmark_summary_report.html"
    ref_status = (suit.classification if suit is not None else "UNKNOWN")
    return {
        "reference": ref_status,
        "hierarchy": (_component_status(suit, "hierarchy_quality", "UNKNOWN")
                      if suit is not None else
                      ("PASS" if mp is not None else "UNKNOWN")),
        "separability": (_component_status(suit, "fine_label_separability", "UNKNOWN")
                         if suit is not None else "UNKNOWN"),
        "input": (_component_status(suit, "gene_overlap", "UNKNOWN")
                  if suit is not None else "UNKNOWN"),
        "bulk_available": bulk_ok,
        "spatial_available": spatial_ok,
        "benchmark_available": bench.exists(),
    }


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


def _caption_for(stem: str, section: str, captions: dict) -> dict:
    """Resolve a *specific* caption for a figure.

    Order: exact stem match -> prefix match -> a non-generic caption built from
    the humanized title.  Never returns the old generic "Visual summary of a
    TissueResolve output" text.
    """
    if stem in captions:
        return captions[stem]
    for key, cap in captions.items():
        if key not in ("_default",) and stem.startswith(key):
            return cap
    title = stem.replace("_", " ").strip()
    title = title[:1].upper() + title[1:] if title else "Figure"
    return {
        "subtitle": f"{section} diagnostic",
        "caption": (f"{title}. Generated from the harmonized reference and "
                    f"{section} query data; values are RNA-derived estimates "
                    "(not absolute cell counts). Axes, colors and values are "
                    "given in the linked source data (.tsv). Colors follow the "
                    "deterministic hierarchical cell-type palette used "
                    "throughout this report."),
        "how_to_read": ("Compare the shown pattern with the QC verdicts in "
                        "sections 2–4; treat fine subtypes cautiously where "
                        "separability is low."),
        "methodology": ("Produced by TissueResolve plotting from saved outputs; "
                        "see the section methodology box and the figure manifest "
                        "for parameters."),
    }


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
        cap = _caption_for(stem, section, captions)
        data = fh.with_suffix(".data.tsv")
        links = [("interactive figure", rel(fh))]
        if data.exists():
            links.append(("source data (.tsv)", rel(data)))
        for ext in ("png", "svg", "pdf"):
            ex = fh.with_suffix("." + ext)
            if ex.exists():
                links.append((f"{ext.upper()}", rel(ex)))
        title = stem.replace("_", " ").strip().capitalize()
        # Embed the figure inline: prefer a static PNG (<img>, always renders);
        # otherwise embed the interactive Plotly HTML via <iframe> so the figure
        # body is never empty.
        png = fh.with_suffix(".png")
        if png.exists():
            body_html = (f"<img src='{rel(png)}' alt='{title}' loading='lazy' "
                         "style='max-width:100%;height:auto;border:1px solid "
                         "#e3e8ef;border-radius:6px'/>")
        else:
            body_html = (f"<iframe src='{rel(fh)}' title='{title}' loading='lazy' "
                         "style='width:100%;height:480px;border:1px solid "
                         "#e3e8ef;border-radius:6px'></iframe>")
        cards.append(C.figure_card(
            title=title, subtitle=cap["subtitle"], body_html=body_html,
            caption=cap["caption"], how_to_read=cap["how_to_read"],
            methodology=cap.get("methodology", ""),
            source_links=links))
        png_rel = rel(png) if png.exists() else ""
        manifest.add(FigureRecord(
            figure_id=stem, section=section, title=title,
            html_path=rel(fh), png_path=png_rel,
            source_data_path=rel(data) if data.exists() else "",
            caption=cap.get("caption", ""),
            methodology=cap.get("methodology", ""), status="generated"))
    return "".join(cards)


def _build_mapping(ref):
    """Build the fine→broad mapping from the documented hierarchy (or {})."""
    try:
        from tissueresolve.reference.hierarchy import (
            build_cell_type_hierarchy, load_hierarchy_mapping)
        hmap_p = H.HARNESS_DIR / "config" / "breast_cancer_cell_type_hierarchy.tsv"
        if ref is not None and hmap_p.exists():
            return build_cell_type_hierarchy(list(ref.cell_types),
                                             load_hierarchy_mapping(hmap_p))
    except Exception:
        pass
    return {}


def generate_diagnostic_figures(ref) -> list[dict]:
    """Generate the reference-QC, signature-QC and resolution figures.

    Writes figures (HTML + PNG + .data.tsv) into the section ``figures/`` dirs so
    the unified report embeds them automatically.  Returns a list of
    ``{figure_id, section, reason}`` for figures whose inputs were missing
    (recorded as ``missing_data`` in the manifest) — never crashes.
    """
    missing: list[dict] = []

    def _safe(fig_id, section, fn):
        try:
            fn()
        except Exception as exc:  # missing data / optional input
            missing.append({"figure_id": fig_id, "section": section,
                            "reason": str(exc)[:300]})

    mp = _build_mapping(ref)
    # cell counts per type (prefer the saved table; fall back to the ref object)
    counts = {}
    ctc = _read_tsv(H.OUT_REFERENCE_DIR / "cell_type_counts.tsv")
    if ctc is not None and "n_cells" in ctc.columns:
        counts = ctc["n_cells"].to_dict()
    elif ref is not None and getattr(ref, "n_cells_per_type", None):
        counts = dict(ref.n_cells_per_type)

    sep = _read_tsv_plain(H.OUTPUTS_DIR / "resolution" / "pairwise_separability.tsv")
    hsum = _hierarchy_summary()
    unresolved = (hsum.get("bulk_unresolved_families")
                  or hsum.get("spatial_unresolved_families") or [])

    # ---- Reference QC figures (section 2) ----
    ref_fig = H.OUT_REFERENCE_DIR / "figures"
    from tissueresolve.plotting import reference_qc_plots as RQ
    _safe("reference_broad_family_composition", "reference",
          lambda: RQ.plot_reference_broad_family_composition(counts, mp, ref_fig))
    _safe("reference_fine_subpopulation_support", "reference",
          lambda: RQ.plot_reference_fine_subpopulation_support(counts, mp, ref_fig))
    _safe("reference_celltype_imbalance", "reference",
          lambda: RQ.plot_reference_celltype_imbalance(counts, ref_fig))
    _safe("gene_overlap_by_modality", "reference",
          lambda: RQ.plot_gene_overlap_by_modality(_gene_overlap_by_modality(), ref_fig))
    comp = _read_tsv_plain(H.OUT_REFERENCE_DIR / "reference_suitability_components.tsv")
    _safe("reference_suitability_components", "reference",
          lambda: RQ.plot_reference_suitability_components(comp, ref_fig))

    # ---- Signature quality & hierarchy figures (section 3) ----
    sig_fig = H.OUTPUTS_DIR / "signature" / "figures"
    from tissueresolve.plotting import signature_qc_plots as SQ
    if ref is not None and getattr(ref, "R_log", None) is not None:
        _safe("signature_matrix_heatmap", "signature",
              lambda: SQ.plot_signature_matrix_heatmap(
                  ref.R_log, ref.gene_names, ref.cell_types, mp, sig_fig))
    else:
        missing.append({"figure_id": "signature_matrix_heatmap",
                        "section": "signature",
                        "reason": "reference signature matrix (R_log) unavailable"})
    _safe("top_confusable_pairs", "signature",
          lambda: SQ.plot_top_confusable_pairs(sep, sig_fig))
    _safe("within_vs_between_family_separability", "signature",
          lambda: SQ.plot_within_vs_between_family_separability(sep, mp, sig_fig))
    _safe("hierarchy_map", "signature",
          lambda: SQ.plot_hierarchy_map(counts, mp, sig_fig,
                                        unresolved_families=unresolved))
    _safe("marker_support_by_family", "signature",
          lambda: SQ.plot_marker_support_by_family(sep, mp, sig_fig))

    # ---- Resolution / spillover / unresolved figures (section 8) ----
    res_fig = H.OUT_RESOLUTION_DIR / "figures"
    from tissueresolve.plotting import resolution_plots as RES
    _safe("separability_distribution", "resolution",
          lambda: RES.plot_separability_distribution(sep, res_fig))
    um = _read_tsv_plain(H.OUTPUTS_DIR / "hierarchical" / "bulk_unresolved_family_mass.tsv")
    _safe("unresolved_mass_by_family", "resolution",
          lambda: RES.plot_unresolved_mass_by_family(um, res_fig))
    _safe("trusted_resolution_summary", "resolution",
          lambda: RES.plot_trusted_resolution_summary(
              mp, sep, res_fig, unresolved_families=unresolved))
    return missing


def generate_benchmark_figures() -> list[dict]:
    """Generate benchmark figures from the benchmark harness TSVs (if present)."""
    missing: list[dict] = []
    bench_out = H.HARNESS_DIR.parent.parent / "benchmarks" / "outputs"
    fig_dir = H.OUTPUTS_DIR / "benchmark" / "figures"
    status_p = bench_out / "real_external_method_status.tsv"
    comp_p = bench_out / "composite_scores.tsv"
    if not status_p.exists() and not comp_p.exists():
        return [{"figure_id": "benchmark_figures", "section": "benchmark",
                 "reason": "no benchmark outputs found (run benchmarks/ first)"}]
    from tissueresolve.plotting import benchmark_report_plots as BP

    def _safe(fig_id, fn):
        try:
            fn()
        except Exception as exc:
            missing.append({"figure_id": fig_id, "section": "benchmark",
                            "reason": str(exc)[:300]})

    status = _read_tsv_plain(status_p)
    comp = _read_tsv_plain(comp_p)
    if status is not None:
        _safe("benchmark_method_status_summary",
              lambda: BP.plot_benchmark_method_status(status, fig_dir))
        _safe("bulk_accuracy_leaderboard",
              lambda: BP.plot_bulk_accuracy_leaderboard(status, fig_dir))
        _safe("runtime_comparison",
              lambda: BP.plot_runtime_comparison(status, fig_dir))
    if comp is not None:
        _safe("composite_scorecard",
              lambda: BP.plot_composite_scorecard(comp, fig_dir))
    return missing


def _read_tsv_plain(p):
    """Read a TSV with a comment header but WITHOUT forcing an index column."""
    try:
        return pd.read_csv(p, sep="\t", comment="#")
    except Exception:
        return None


def _gene_overlap_by_modality() -> dict:
    """Collect {modality: {n_reference, n_query, n_shared}} from QC warnings."""
    out = {}
    for modality, path in (("bulk", H.OUT_BULK_DIR / "bulk_warnings.json"),
                           ("spatial", H.OUT_SPATIAL_DIR / "spatial_warnings.json")):
        if path.exists():
            try:
                go = (json.loads(_read_text(path)) or {}).get("gene_overlap") or {}
                if go:
                    out[modality] = {"n_reference": go.get("n_reference", 0),
                                     "n_query": go.get("n_query", 0),
                                     "n_shared": go.get("n_shared", 0)}
            except Exception:
                pass
    return out


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
    # Compute reference suitability ONCE (reused by the executive decision
    # summary AND the reference-quality section) so headline statuses reflect
    # this run, not a stale file.
    suit = None
    try:
        from tissueresolve.reference.suitability import (
            compute_reference_suitability_score, save_reference_suitability)
        if ref0 is not None:
            qgenes = (list(pd.read_csv(H.PSEUDOBULK_COUNTS, sep="\t", index_col=0,
                                       comment="#").index.map(str))
                      if H.PSEUDOBULK_COUNTS.exists() else None)
            with _w.catch_warnings():
                _w.simplefilter("ignore")
                suit = compute_reference_suitability_score(
                    ref0, query_genes=qgenes, mapping=mp)
                save_reference_suitability(suit, H.OUT_REFERENCE_DIR)
    except Exception:
        suit = None
    ctx["suit"] = suit
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

    # ---- 1. Executive decision summary ----
    st = _decision_statuses(ctx, suit, mp)
    hs = _hierarchy_summary()
    n_unres = len(hs.get("bulk_unresolved_families")
                  or hs.get("spatial_unresolved_families") or [])
    yn = lambda b: "yes" if b else "no"  # noqa: E731
    status_cards = C.decision_status_grid([
        ("Reference suitability", st["reference"],
         "is the reference good enough?"),
        ("Hierarchy", st["hierarchy"], "broad/fine labels usable?"),
        ("Signature separability", st["separability"],
         "can cell types be told apart?"),
        ("Input compatibility", st["input"], "query vs reference overlap"),
        ("Bulk results", yn(st["bulk_available"]), "available?"),
        ("Spatial results", yn(st["spatial_available"]), "available?"),
        ("Benchmark", yn(st["benchmark_available"]), "available?"),
    ])
    checklist = C.checklist([
        ("Gene overlap between query and reference checked",
         ctx["overlap"] not in ("—", None)),
        ("Broad/fine cell-type labels present", mp is not None),
        ("Minimum cells per cell type checked", suit is not None),
        ("Cell-type imbalance checked", suit is not None),
        ("Signature separability checked",
         (H.OUTPUTS_DIR / "resolution" / "pairwise_separability.tsv").exists()),
        ("Spillover risk checked",
         (H.OUTPUTS_DIR / "resolution" / "spillover_risk_by_celltype.tsv").exists()),
        ("Solver selected", bool(ctx["solver"])),
        ("Spatial 'no ground truth' caveat shown",
         st["spatial_available"]),
    ])
    cards = {"Modality": ctx["modality"], "Samples (bulk)": ctx["samples"],
             "Spots (spatial)": ctx["spots"], "Reference cells": ctx["ref_cells"],
             "Genes": ctx["genes"], "Broad families": ctx["broad"],
             "Fine subpopulations": ctx["fine"], "Gene overlap": ctx["overlap"],
             "Solver": ctx["solver"], "Resolution mode": ctx["resolution_mode"],
             "Unresolved families": n_unres or "—", "Benchmark": bench_status}
    guide_bits = [
        "Check the QC panels (sections 2–4) <b>before</b> interpreting any "
        "deconvolution result.",
    ]
    if str(st["reference"]).upper() in ("WARNING", "FAIL"):
        guide_bits.append(
            f"Reference suitability is <b>{C.esc(st['reference'])}</b> — interpret "
            "predictions cautiously.")
    if str(st["separability"]).upper() in ("CAUTION", "WARNING", "FAIL"):
        guide_bits.append(
            "Fine-level separability is limited — prefer family-level results "
            "where subtypes are not separable.")
    summ = (C.estimate_note(
                "<b>Analysis status.</b> The cards below are quality-control "
                "verdicts, not biological conclusions.")
            + status_cards
            + "<h3 style='margin:14px 0 4px;font-size:14px'>Before interpreting "
              "results</h3>" + checklist
            + C.metric_grid(cards)
            + C.estimate_note(
                "This report provides QC, estimates, diagnostics, and visualization "
                "tools. It does <b>not</b> determine biological causality or "
                "validate cell identities without external evidence. Values shown "
                "are <b>RNA-derived estimates</b>.")
            + C.interpretation_guide(" ".join(guide_bits)))
    md = H.OUT_SUMMARY_DIR / "validation_summary.md"
    if md.exists():
        summ += C.collapsible_table("Text validation summary",
                                    "<pre>" + C.esc(_read_text(md)) + "</pre>")
    sections.append(Section("summary", "1. Executive decision summary", summ))

    # ---- 2. Reference quality ----
    refq = C.methodology_summary([
        "Reference built from a single-cell/nucleus h5ad by aggregating per-cell "
        "profiles into per-cell-type signatures (CPM).",
        "Broad and fine labels taken from the documented hierarchy mapping; gene "
        "identifiers harmonized to symbols.",
        "Minimum cells per type enforced; cross-donor variability recorded when "
        "multiple donors are present.",
    ])
    if suit is not None:
        try:
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

    # ---- 3. Signature quality & hierarchy (QC, BEFORE any result) ----
    sig = C.methodology_summary([
        "Pairwise separability is computed from the reference expression profiles "
        "(1 − Bhattacharyya coefficient); higher = more distinguishable.",
        "Hierarchy usability checks whether broad families can be split into "
        "reliably separable fine subtypes.",
        "This section answers 'are the signatures good enough to tell cell types "
        "apart?' and must be read BEFORE the deconvolution results.",
    ])
    if suit is not None:
        sig += ("<p>Signature separability: "
                f"{C.status_badge(_component_status(suit, 'fine_label_separability', 'UNKNOWN'))}"
                " &nbsp; Hierarchy quality: "
                f"{C.status_badge(_component_status(suit, 'hierarchy_quality', 'UNKNOWN'))}"
                "</p>")
    # Signature-QC figures (heatmap, confusable pairs, within/between, hierarchy,
    # marker support) — visual summaries before the raw table.
    sig += _figure_cards(H.OUTPUTS_DIR / "signature" / "figures", out_dir,
                         manifest, "signature", _CAPTIONS)
    # Top confusable pairs (compact view; full table is source data)
    sep_p = H.OUTPUTS_DIR / "resolution" / "pairwise_separability.tsv"
    if sep_p.exists():
        try:
            sdf = pd.read_csv(sep_p, sep="\t", comment="#")
            keep = [c for c in ("type_a", "type_b", "separability_score",
                                "bhattacharyya", "resolvability") if c in sdf.columns]
            sort_col = "separability_score" if "separability_score" in sdf.columns \
                else ("bhattacharyya" if "bhattacharyya" in sdf.columns else keep[0])
            asc = sort_col == "separability_score"
            top = sdf.sort_values(sort_col, ascending=asc)[keep].head(15)
            sig += ("<p class='muted'>The 15 hardest-to-distinguish cell-type "
                    f"pairs (of {len(sdf)} total). Low separability means fine "
                    "labels for these types are unreliable.</p>")
            sig += C.collapsible_table(
                "15 most confusable cell-type pairs",
                top.to_html(index=False, border=0), open=True)
        except Exception:
            pass
    # Hierarchy usability summary
    if hs:
        fams = hs.get("families") or []
        unres = (hs.get("bulk_unresolved_families")
                 or hs.get("spatial_unresolved_families") or [])
        sig += C.metric_grid({
            "Broad families": len(fams) or "—",
            "Fine subpopulations": hs.get("n_fine", ctx["fine"]),
            "Families resolved to subtypes": (len(fams) - len(unres)) if fams else "—",
            "Families kept at family level": len(unres) or "—",
        })
        if unres:
            sig += C.warning_box([
                "These families cannot be split into separable subtypes at this "
                "reference's resolution and are reported at the family level "
                "(unresolved mass): " + ", ".join(map(str, unres)) + "."])
    sig += C.variable_dictionary(G.subset(["separability", "spillover", "unresolved mass"]))
    sig += C.interpretation_guide(
        "If signature separability is CAUTION/WARNING/FAIL, or a family is listed "
        "as unresolved, interpret those cell types at the family level rather than "
        "as confident subtypes. Fine subtype labels should be treated cautiously "
        "when separability is low.")
    siglinks = [("resolution outputs", rel(H.OUTPUTS_DIR / "resolution"))] \
        if (H.OUTPUTS_DIR / "resolution").exists() else []
    sections.append(Section("signature", "3. Signature quality & hierarchy",
                            sig, links=siglinks))

    # ---- 4. Input data ----
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
    sections.append(Section("input", "4. Input data quality", inp))

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
    sections.append(Section("bulk", "5. Bulk deconvolution", bulk, links=blinks))

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
    sections.append(Section("spatial", "6. Spatial deconvolution", spat, links=slinks))

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
    sections.append(Section("hierarchical", "7. Hierarchical broad→fine deconvolution",
                            hier, links=hlinks))

    # ---- 7. Resolution / separability / spillover ----
    res_dir = H.OUT_RESOLUTION_DIR
    resb = C.methodology_summary([
        "Pairwise separability computed from the reference expression profiles "
        "(1 − Bhattacharyya coefficient).",
        "High-risk pairs (BC > 0.90) and per-family resolution summaries highlight "
        "where fine labels are unreliable.",
    ])
    resb += _figure_cards(res_dir / "figures", out_dir, manifest,
                          "resolution", _CAPTIONS)
    resb += _df_collapsible("Top non-separable pairs",
                            res_dir / "pairwise_separability.tsv")
    resb += C.variable_dictionary(G.subset(["separability", "spillover", "unresolved mass"]))
    resb += C.interpretation_guide(
        "Use this to judge what resolution you can trust: families with many "
        "non-separable subtypes should be interpreted at the broad level.")
    rlinks = [("resolution outputs", rel(res_dir))] if res_dir.exists() else []
    sections.append(Section("resolution", "8. Resolution, separability & spillover",
                            resb, links=rlinks))

    # ---- 8. Benchmark ----
    benchb = C.methodology_summary([
        "Methods compared on identical harmonized inputs.",
        "Only EXECUTED or IMPORTED methods are ranked; EXPORTED-only and SKIPPED "
        "tools are listed but not scored.",
        "Bulk pseudobulk has ground truth (accuracy valid); real Visium does not "
        "(concordance/structure only).",
        "Measured metrics (status, accuracy, runtime) come first; the weighted "
        "composite is a scorecard, not objective accuracy.",
    ])
    benchb += _figure_cards(H.OUTPUTS_DIR / "benchmark" / "figures", out_dir,
                            manifest, "benchmark", _CAPTIONS)
    benchb += _df_collapsible("Best-method summary",
                              bench.parent / "bulk" / "best_method_summary.tsv")
    benchb += C.variable_dictionary(G.subset(
        ["Pearson correlation", "RMSE", "concordance", "composite score"]))
    benchb += C.interpretation_guide(
        "Use this section to compare methods. Only executed or imported methods "
        "are ranked; exported-only tools are listed but not scored.")
    blinks2 = [("benchmark summary report", rel(bench))] if bench.exists() else []
    sections.append(Section("benchmark", "9. Benchmark comparison", benchb, links=blinks2))

    # ---- 9. Warnings & limitations ----
    # Aggregate the run's real warnings from QC, spillover and separability
    # outputs (not just the figure-generation log) so the section is honest.
    warn_items = _collect_warnings()
    wbody = C.warning_box(warn_items) if warn_items else \
        "<p class='muted'>No warnings recorded for this run.</p>"
    wbody += C.limitation_box([
        "Bulk estimates are RNA-derived proportions, not absolute cell fractions.",
        "Spatial estimates are spot-level RNA-derived composition, not cell counts; "
        "real Visium has no ground truth.",
        "Fine subtype estimates in non-separable families are reported as "
        "unresolved mass; do not over-interpret them.",
    ])
    sections.append(Section("warnings", "10. Warnings & limitations", wbody))

    # ---- 10. Methods ----
    methods = H.OUT_SUMMARY_DIR / "methods.txt"
    mbody = ("<pre>" + C.esc(_read_text(methods)) + "</pre>" if methods.exists()
             else "<p class='muted'>Methods text not available.</p>")
    sections.append(Section("methods", "11. Methods", mbody))

    # ---- 11. Output files & source data ----
    # Record figures that could not be generated (missing input data) so the
    # manifest is honest rather than silently omitting them.
    from tissueresolve.report.figures import FigureRecord, STATUS_MISSING_DATA
    mf = H.OUT_SUMMARY_DIR / "missing_figures.json"
    n_missing = 0
    if mf.exists():
        try:
            for rec in (json.loads(_read_text(mf)) or []):
                manifest.add(FigureRecord(
                    figure_id=str(rec.get("figure_id", "")),
                    section=str(rec.get("section", "")),
                    title=str(rec.get("figure_id", "")).replace("_", " ").capitalize(),
                    status=STATUS_MISSING_DATA,
                    reason_if_missing=str(rec.get("reason", ""))))
                n_missing += 1
        except Exception:
            pass
    man_path = manifest.save(out_dir / "figures")
    n_gen = len(manifest.records) - n_missing
    out_files = sorted(p for p in out_dir.rglob("*.tsv"))[:100]
    obody = (C.estimate_note(
                 f"Figure manifest: {rel(man_path)} — {n_gen} figures generated, "
                 f"{n_missing} recorded as missing_data (with reasons).")
             + C.metric_grid({"Figures generated": n_gen,
                              "Figures missing data": n_missing,
                              "Source-data tables": len(out_files)})
             + C.collapsible_table("All source-data tables (.tsv)",
                                   "<ul>" + "".join(f"<li>{C.esc(rel(p))}</li>"
                                                    for p in out_files) + "</ul>"))
    sections.append(Section("outputs", "12. Output files & source data", obody))

    path = build_unified_report(
        out_dir / "report.html", sections,
        title="TissueResolve — breast-cancer analysis report",
        subtitle="RNA-derived estimates with QC, resolution limits and benchmarks")
    return path


# Per-figure captions: title is derived from filename; these add the
# subtitle / full caption / how-to-read / methodology.  Unlisted figures get a
# specific (non-generic) caption built from their humanized title via
# `_caption_for`, so no figure carries a generic "TissueResolve output" caption.
_CAPTIONS = {
    # ---- Reference QC (section 2) ----
    "reference_broad_family_composition": {
        "subtitle": "reference composition by broad family",
        "caption": "Reference composition by broad cell family. Horizontal bars "
                   "show the number of reference single-cell/nucleus profiles per "
                   "broad family (x-axis = cells; y-axis = family), sorted by "
                   "abundance; colours are the deterministic broad-family palette "
                   "reused throughout the report. Source: cell_type_counts.tsv + "
                   "hierarchy mapping. A reference dominated by one family, or "
                   "lacking support for rare families, can yield less stable "
                   "signatures.",
        "how_to_read": "Check whether one family dominates or a needed family is "
                       "barely represented.",
        "methodology": "Per-cell-type counts aggregated to broad families via the "
                       "documented hierarchy mapping.",
    },
    "reference_fine_subpopulation_support": {
        "subtitle": "fine-label support within families",
        "caption": "Reference support for fine subpopulations. Bars show reference "
                   "cells per fine label (x-axis), grouped by broad family and "
                   "coloured with the family's subtone; the dashed line marks the "
                   "minimum-cell support threshold. Source: cell_type_counts.tsv + "
                   "mapping. Fine labels below the threshold have weak signatures.",
        "how_to_read": "Treat fine labels below the threshold line cautiously.",
        "methodology": "Per-fine-label counts; family subtones from the "
                       "hierarchical palette.",
    },
    "reference_celltype_imbalance": {
        "subtitle": "abundance imbalance",
        "caption": "Reference imbalance across cell types. Cell types are ranked "
                   "by abundance (x-axis) with cell counts (y-axis); the subtitle "
                   "reports a Gini imbalance score (0 = even, 1 = one type "
                   "dominates). Source: cell_type_counts.tsv. Strong imbalance can "
                   "destabilise rare-type signatures.",
        "how_to_read": "A high Gini or one towering bar means rare types are "
                       "under-supported.",
        "methodology": "Counts ranked; Gini computed from the count distribution.",
    },
    "gene_overlap_by_modality": {
        "subtitle": "reference vs query gene overlap",
        "caption": "Gene overlap between reference and query data. Grouped bars "
                   "show reference, query and shared gene counts per input "
                   "modality (bulk / spatial). Source: per-modality gene_overlap "
                   "in the QC warnings. Only shared genes inform deconvolution; "
                   "low overlap weakens estimates.",
        "how_to_read": "Check the shared-gene bar is a large fraction of the "
                       "reference genes.",
        "methodology": "Gene-set intersection computed during QC; nothing is "
                       "silently re-normalized.",
    },
    "reference_suitability_components": {
        "subtitle": "suitability traffic light",
        "caption": "Reference suitability components. Horizontal bars per component "
                   "(gene overlap, protocol, library, batch, marker stability, "
                   "cell-type balance, separability, hierarchy) coloured "
                   "PASS/CAUTION/WARNING/FAIL/UNKNOWN; bar length = component score "
                   "(full bar when N/A). Source: reference_suitability_components."
                   "tsv. WARNING/FAIL components should be addressed before "
                   "trusting fine predictions.",
        "how_to_read": "Focus on WARNING/FAIL components; UNKNOWN means metadata "
                       "was not provided.",
        "methodology": "Reference suitability scoring (see suitability report).",
    },
    # ---- Signature quality & hierarchy (section 3) ----
    "signature_matrix_heatmap": {
        "subtitle": "top signature genes × cell types",
        "caption": "Reference signature matrix. Heatmap of the most variable "
                   "signature genes (rows) across cell types (columns, grouped by "
                   "broad family); colour is the row z-scored signature value "
                   "(display only). Source: reference signature matrix (full "
                   "matrix is the saved reference). Distinct column blocks mean "
                   "cell types have separable signatures.",
        "how_to_read": "Look for genes that are high in one cell type and low in "
                       "others (good markers); flat rows are uninformative.",
        "methodology": "Top-variance genes selected; rows z-scored for display.",
    },
    "top_confusable_pairs": {
        "subtitle": "least-separable pairs",
        "caption": "Most difficult-to-distinguish cell-type pairs. Lollipops show "
                   "the least-separable pairs (x-axis = 1 − Bhattacharyya; y-axis "
                   "= pair), coloured by severity (red = worst). Source: "
                   "pairwise_separability.tsv. Low-separability pairs have "
                   "unreliable fine labels — interpret at the family level.",
        "how_to_read": "Pairs near 0 (red) cannot be reliably told apart.",
        "methodology": "Pairs ranked by ascending separability (1 − Bhattacharyya).",
    },
    "within_vs_between_family_separability": {
        "subtitle": "within vs between family",
        "caption": "Separability within and between broad families. Box + points "
                   "of separability (y-axis) for within-family vs between-family "
                   "pairs. Source: pairwise_separability.tsv + hierarchy mapping. "
                   "Lower within-family separability means subtypes of a family "
                   "are hard to resolve.",
        "how_to_read": "If the within-family box sits much lower, prefer "
                       "family-level results.",
        "methodology": "Pairs grouped by whether the two types share a broad "
                       "family.",
    },
    "hierarchy_map": {
        "subtitle": "broad→fine annotation tree",
        "caption": "Broad-to-fine annotation hierarchy. Treemap where broad "
                   "families contain their fine subpopulations and box size = "
                   "number of reference cells; colours are the hierarchical "
                   "palette and ⚠ marks families reported only at the family level "
                   "(unresolved). Source: cell_type_counts.tsv + mapping + "
                   "unresolved-family list.",
        "how_to_read": "Large boxes are well-supported; ⚠ families should be read "
                       "at the family level.",
        "methodology": "Counts arranged as a broad→fine treemap.",
    },
    "marker_support_by_family": {
        "subtitle": "subtype marker support",
        "caption": "Subtype marker support within each broad family. Bars show the "
                   "mean number of discriminating markers for within-family "
                   "subtype pairs (x-axis) per family; the dashed line marks weak "
                   "support. Source: pairwise_separability.tsv "
                   "(n_discriminating_genes). Families below the line have weak "
                   "subtype marker support.",
        "how_to_read": "Families below the threshold cannot reliably separate "
                       "their subtypes.",
        "methodology": "Mean discriminating-gene count over within-family pairs.",
    },
    # ---- Resolution / spillover / unresolved (section 8) ----
    "separability_distribution": {
        "subtitle": "global separability histogram",
        "caption": "Distribution of pairwise cell-type separability. Histogram of "
                   "separability (x-axis = 1 − Bhattacharyya; y-axis = pair count) "
                   "with a dashed high-risk threshold. Source: "
                   "pairwise_separability.tsv. Mass below the threshold indicates "
                   "many confusable pairs.",
        "how_to_read": "A large left tail means many cell types are hard to "
                       "separate.",
        "methodology": "Histogram of all pairwise separability scores.",
    },
    "unresolved_mass_by_family": {
        "subtitle": "mean unresolved mass",
        "caption": "Unresolved mass by broad family. Bars show the mean "
                   "RNA-derived mass not split into subtypes (x-axis) per family. "
                   "Source: bulk_unresolved_family_mass.tsv. Higher values mean "
                   "more mass stayed at the family level because subtypes were not "
                   "separable.",
        "how_to_read": "High-mass families should be interpreted at the family "
                       "level.",
        "methodology": "Per-sample unresolved mass averaged across samples.",
    },
    "trusted_resolution_summary": {
        "subtitle": "recommended interpretation level",
        "caption": "Recommended interpretation level by family. Bars rank each "
                   "family's recommended resolution (1 = broad/family-level … "
                   "4 = fine), coloured by level, combining within-family "
                   "separability and unresolved-family status. Source: mapping + "
                   "separability + unresolved list. Families marked "
                   "broad/family-level should not be read as confident subtypes.",
        "how_to_read": "Read each family at (or below) its recommended level.",
        "methodology": "Worst within-family separability + unresolved status → "
                       "recommended level.",
    },
    # ---- Benchmark (section 9) ----
    "benchmark_method_status_summary": {
        "subtitle": "tools by status",
        "caption": "Benchmark method status. Bars count methods per status "
                   "(executed / imported / exported-only / skipped / failed). "
                   "Source: real_external_method_status.tsv. Only executed/"
                   "imported tools carry measured metrics; the rest are listed but "
                   "not scored.",
        "how_to_read": "Check how many tools actually ran vs were skipped/failed.",
        "methodology": "Counts of the benchmark status column.",
    },
    "bulk_accuracy_leaderboard": {
        "subtitle": "measured bulk accuracy",
        "caption": "Bulk benchmark accuracy. Bars show measured accuracy (x-axis) "
                   "for executed/imported bulk methods only, against pseudobulk "
                   "ground truth. Source: benchmark status/metrics. Spatial "
                   "methods and skipped/failed tools are excluded (no comparable "
                   "ground truth).",
        "how_to_read": "Compare measured accuracy across methods; higher is "
                       "better.",
        "methodology": "Accuracy vs known pseudobulk proportions; executed/"
                       "imported only.",
    },
    "runtime_comparison": {
        "subtitle": "runtime (seconds)",
        "caption": "Benchmark runtime. Bars show wall-clock runtime in seconds "
                   "(x-axis) for executed/imported methods. Source: benchmark "
                   "status runtime_seconds. Runtimes depend on this environment "
                   "and settings (e.g. fast mode).",
        "how_to_read": "Compare cost; very fast or very slow tools stand out.",
        "methodology": "Recorded wall-clock time per executed/imported method.",
    },
    "composite_scorecard": {
        "subtitle": "weighted scorecard (NOT accuracy)",
        "caption": "Composite benchmark scorecard. Stacked bars show each scored "
                   "method's weighted contribution per dimension (accuracy, "
                   "robustness, usability, interpretability, resolution-awareness, "
                   "runtime). Source: composite_scores.tsv. This is a weighted "
                   "scorecard, NOT an objective accuracy measure; bulk and spatial "
                   "are scored separately.",
        "how_to_read": "Use as a multi-criteria summary, not a single accuracy "
                       "ranking.",
        "methodology": "Weighted sum of per-dimension scores (see composite-score "
                       "weights).",
    },
    "bulk_composition_clustered_barplot": {
        "subtitle": "bulk composition by sample",
        "caption": "Bulk RNA-derived composition by sample. Each stacked bar is "
                   "one pseudobulk sample; segment height is the estimated mRNA "
                   "proportion of a cell type/family (bar sums to 1). Colors are "
                   "the hierarchical cell-type palette; samples are ordered by "
                   "hierarchical clustering of their composition. RNA-derived "
                   "proportions, NOT absolute cell fractions.",
        "how_to_read": "Look for sample groups with similar composition and for "
                       "families that dominate a sample. Interpret subtypes from "
                       "low-separability families at the family level only.",
        "methodology": "Weighted-NNLS family of solvers (solver=auto, chosen by "
                       "gene-masking cross-validation); proportions normalized per "
                       "sample.",
    },
    "bulk_composition_heatmap": {
        "subtitle": "composition heatmap",
        "caption": "Bulk composition heatmap. Rows are cell types/families, "
                   "columns are samples; cell color encodes the RNA-derived "
                   "proportion (legend on the figure). Rows and columns ordered by "
                   "clustering. RNA-derived proportions, not absolute cell counts.",
        "how_to_read": "Scan for blocks of co-varying cell types across samples; "
                       "compare with the clustered barplot.",
        "methodology": "Same per-sample proportions as the composition barplot, "
                       "shown as a matrix.",
    },
    "bulk_qc_summary": {
        "subtitle": "per-sample reconstruction QC",
        "caption": "Bulk per-sample QC. Bars/points show reconstruction quality "
                   "(coverage R² / profile correlation) and mismatch flags per "
                   "sample. Higher reconstruction = the reference explains the "
                   "sample better.",
        "how_to_read": "Samples with low reconstruction or a mismatch flag should "
                       "be interpreted cautiously.",
        "methodology": "QC metrics recomputed from the saved fit; see bulk_qc.tsv.",
    },
    "bulk_uncertainty_plot": {
        "subtitle": "estimate uncertainty",
        "caption": "Bulk estimate uncertainty. Shows per-cell-type confidence "
                   "intervals where bootstrap/CV uncertainty was computed; wider "
                   "intervals mean less certain estimates.",
        "how_to_read": "Prefer cell types with narrow intervals; treat wide-"
                       "interval estimates as indicative only.",
        "methodology": "Bootstrap/gene-masking CV when enabled; otherwise the "
                       "figure notes that uncertainty was not computed.",
    },
    "bulk_separability_heatmap": {
        "subtitle": "reference separability",
        "caption": "Reference separability heatmap (bulk view). Rows/columns are "
                   "cell types; color is pairwise separability (1 − Bhattacharyya "
                   "coefficient). Darker/low values mark pairs that are hard to "
                   "tell apart in the reference.",
        "how_to_read": "Clusters of low separability indicate cell types whose "
                       "fine estimates are unreliable — interpret at the family "
                       "level.",
        "methodology": "Bhattacharyya coefficient between per-cell-type expression "
                       "profiles.",
    },
    "spatial_separability_heatmap": {
        "subtitle": "reference separability",
        "caption": "Reference separability heatmap (spatial view). Rows/columns "
                   "are cell types; color is pairwise separability "
                   "(1 − Bhattacharyya). Low values mark confusable pairs.",
        "how_to_read": "Low-separability pairs should be interpreted at the family "
                       "level on the spatial maps.",
        "methodology": "Bhattacharyya coefficient between per-cell-type expression "
                       "profiles.",
    },
    "bulk_spillover_heatmap": {
        "subtitle": "spillover risk",
        "caption": "Spillover-risk heatmap (bulk). Color encodes the correlation "
                   "between cell-type signatures; high values mean signal can leak "
                   "between those types during deconvolution.",
        "how_to_read": "High off-diagonal values flag cell types whose estimates "
                       "may be confounded with collinear partners.",
        "methodology": "Expression-signature correlation proxy for spillover.",
    },
    "spatial_spillover_heatmap": {
        "subtitle": "spillover risk",
        "caption": "Spillover-risk heatmap (spatial). Color encodes signature "
                   "correlation; high values mean signal can leak between types.",
        "how_to_read": "High off-diagonal values flag confounded cell types.",
        "methodology": "Expression-signature correlation proxy for spillover.",
    },
    "spillover_network": {
        "subtitle": "high-risk spillover network",
        "caption": "High-risk spillover network. Nodes are cell types (colored by "
                   "broad family); edges connect pairs with high signature "
                   "correlation (spillover risk), edge weight = risk. Only "
                   "high-risk edges are drawn.",
        "how_to_read": "Tightly connected groups share signal; estimates within "
                       "them are less independent.",
        "methodology": "Edges thresholded on the signature-correlation spillover "
                       "proxy.",
    },
    "spatial_mean_composition_barplot": {
        "subtitle": "tissue-average composition",
        "caption": "Mean RNA-derived composition across the tissue section. Bars "
                   "show the section-averaged proportion per cell type/family "
                   "(colors = hierarchical palette). Spot-level RNA-derived "
                   "composition, not cell counts; real Visium has no ground truth.",
        "how_to_read": "Read as the average tissue makeup; inspect the spatial "
                       "maps for where each population localizes.",
        "methodology": "Mean over per-spot NB-CAR proportions.",
    },
    "spatial_dominant_cell_type_map": {
        "subtitle": "dominant family per spot",
        "caption": "Dominant predicted population per spot. Each spot is colored "
                   "by its highest-proportion cell type/family (hierarchical "
                   "palette) at its array coordinates. Spot-level RNA-derived "
                   "composition, not direct cell identity.",
        "how_to_read": "Compare regions with H&E morphology; do not treat the "
                       "dominant label as a single-cell call.",
        "methodology": "Argmax of per-spot NB-CAR proportions.",
    },
    "spatial_abundance_maps": {
        "subtitle": "abundance maps",
        "caption": "Spatial abundance maps for top broad families. Each panel "
                   "shows one family's RNA-derived proportion across spots "
                   "(continuous color scale). Spot-level composition, not counts.",
        "how_to_read": "Look for spatially coherent regions of high abundance; "
                       "compare with H&E.",
        "methodology": "Per-spot NB-CAR proportions for the top families by mean "
                       "abundance.",
    },
    "spatial_morans_i_barplot": {
        "subtitle": "spatial autocorrelation",
        "caption": "Spatial autocorrelation (Moran's I) of predicted populations. "
                   "Bars rank cell types/families by Moran's I; higher = more "
                   "spatially clustered, near 0 = spatially random.",
        "how_to_read": "High Moran's I indicates structured localization; it is a "
                       "structure statistic, not an accuracy measure.",
        "methodology": "Moran's I computed on the spot graph from array "
                       "coordinates.",
    },
    "spatial_spot_pie_charts": {
        "subtitle": "exploratory — per-spot pies",
        "caption": "EXPLORATORY: per-spot composition pie charts. Each pie shows a "
                   "spot's RNA-derived composition. Dense and hard to read at "
                   "scale; provided for exploration only, not as a primary view.",
        "how_to_read": "Use the dominant-family and abundance maps as the primary "
                       "spatial views; treat pies as exploratory.",
        "methodology": "Per-spot NB-CAR proportions drawn as pies.",
    },
    "he_spots_check": {
        "subtitle": "H&E + spot alignment",
        "caption": "H&E and Visium spot alignment. Gray spots are overlaid on the "
                   "H&E image at their array coordinates to verify registration "
                   "before interpreting any abundance map.",
        "how_to_read": "Confirm spots fall on tissue and align with morphology "
                       "before trusting spatial overlays.",
        "methodology": "Spot coordinates overlaid on the section H&E.",
    },
    "he_dominant_cell_type": {
        "subtitle": "dominant family over H&E",
        "caption": "Dominant predicted broad family over H&E. Spots colored by "
                   "their top family (hierarchical palette) on the H&E image. "
                   "Spot-level RNA-derived composition, not cell identity.",
        "how_to_read": "Compare colored regions with visible tissue morphology; "
                       "do not treat as single-cell identity.",
        "methodology": "Argmax of per-spot proportions over the H&E.",
    },
    "he_abundance": {  # prefix match for he_abundance_<celltype>
        "subtitle": "abundance over H&E",
        "caption": "Abundance of one broad family over H&E. Spot color encodes "
                   "that family's RNA-derived proportion on the H&E image "
                   "(continuous scale). Spot-level composition, not counts.",
        "how_to_read": "Look for high-abundance regions and compare with "
                       "morphology; this is RNA-derived, not a direct count.",
        "methodology": "Per-spot NB-CAR proportion for the named family over H&E.",
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
    # New informative QC / signature / resolution / benchmark figures.
    missing_figs = generate_diagnostic_figures(ref)
    missing_figs += generate_benchmark_figures()
    H.write_json(missing_figs, H.OUT_SUMMARY_DIR / "missing_figures.json")
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
