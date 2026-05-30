"""
Thin top-level convenience API for TissueResolve.

These functions are **thin wrappers** over the existing, discovered pipeline
APIs — they do not introduce new behaviour, only a stable entry point:

================  =========================================================
deconv_bulk       → :meth:`tissueresolve.bulk.pipeline.BulkPipeline.run`
deconv_spatial    → :meth:`tissueresolve.spatial.pipeline.SpatialPipeline.run`
build_reference   → :class:`tissueresolve.reference.build.ReferenceBuilder`
generate_report   → :mod:`tissueresolve.report.html`
plot_results      → :mod:`tissueresolve.plotting`
================  =========================================================

All heavy / optional dependencies (scanpy, matplotlib) are imported lazily
inside the functions, so ``import tissueresolve`` stays cheap and works
without the optional extras installed.
"""
from __future__ import annotations

from pathlib import Path
from typing import Any, Optional, Sequence, Union

__all__ = [
    "deconv_bulk",
    "deconv_spatial",
    "build_reference",
    "generate_report",
    "plot_results",
]


def deconv_bulk(bulk, ref, *, config=None, **kwargs):
    """Run the bulk deconvolution pipeline.  Returns ``BulkPipelineResult``.

    Wrapper over :meth:`BulkPipeline.run`; extra keyword arguments
    (``gene_panel``, ``n_bootstrap``, ``protocol_meta``, ``mrna_corrector`` …)
    are forwarded verbatim.
    """
    from tissueresolve.bulk.pipeline import BulkPipeline

    return BulkPipeline(config).run(bulk, ref, **kwargs)


def deconv_spatial(
    Y, ref, array_row, array_col, lib_sizes, gene_names,
    spot_ids=None, *, config=None, **kwargs,
):
    """Run the spatial deconvolution pipeline.  Returns ``SpatialPipelineResult``.

    Wrapper over :meth:`SpatialPipeline.run`; extra keyword arguments
    (``marker_genes``, ``run_neighbourhood`` …) are forwarded verbatim.
    """
    from tissueresolve.spatial.pipeline import SpatialPipeline

    return SpatialPipeline(config).run(
        Y, ref, array_row, array_col, lib_sizes, gene_names, spot_ids, **kwargs
    )


def build_reference(
    source,
    *,
    metadata=None,
    config=None,
    cell_type_col: Optional[str] = None,
    estimate_overdispersion: bool = False,
):
    """Build a :class:`ReferenceSignature` from a DataFrame, AnnData, or .h5ad.

    Dispatches on *source* type:

    * ``pd.DataFrame`` → :meth:`ReferenceBuilder.build_from_df` (``metadata``
      required: per-cell cell-type labels);
    * AnnData (duck-typed ``.obs``/``.var``/``.X``) →
      :meth:`ReferenceBuilder.build_from_adata`;
    * ``str``/``Path`` → :meth:`ReferenceBuilder.build_from_h5ad`.
    """
    import pandas as pd

    from tissueresolve.config import TissueResolveConfig
    from tissueresolve.reference.build import ReferenceBuilder

    cfg = config or TissueResolveConfig()
    if cell_type_col is not None:
        cfg.reference.celltype_col = cell_type_col
    builder = ReferenceBuilder(cfg.reference)

    if isinstance(source, pd.DataFrame):
        if metadata is None:
            raise ValueError(
                "build_reference from a DataFrame requires `metadata` "
                "(per-cell cell-type labels)."
            )
        return builder.build_from_df(
            source, metadata, estimate_overdispersion=estimate_overdispersion
        )
    if isinstance(source, (str, Path)):
        return builder.build_from_h5ad(
            str(source), estimate_overdispersion=estimate_overdispersion
        )
    if hasattr(source, "obs") and hasattr(source, "var") and hasattr(source, "X"):
        return builder.build_from_adata(
            source, estimate_overdispersion=estimate_overdispersion
        )
    raise TypeError(
        f"build_reference: unsupported source type {type(source)!r}.  "
        "Expected DataFrame, AnnData, or a path to an .h5ad file."
    )


def generate_report(
    result: Any,
    output_path: Union[str, Path],
    **kwargs,
) -> Path:
    """Write an HTML report, dispatching on bulk vs spatial result type.

    Accepts a ``BulkPipelineResult`` or ``SpatialPipelineResult`` (detected by
    the embedded ``deconv`` estimate type).  Extra kwargs (``separability``,
    ``figures``, ``output_files`` …) are forwarded to the report generator.
    """
    from tissueresolve.report import html
    from tissueresolve.results import BulkDeconvResult, SpatialDeconvResult

    deconv = getattr(result, "deconv", None)
    if isinstance(deconv, BulkDeconvResult):
        return html.generate_bulk_report(result, output_path, **kwargs)
    if isinstance(deconv, SpatialDeconvResult):
        return html.generate_spatial_report(result, output_path, **kwargs)
    raise TypeError(
        "generate_report: result must be a BulkPipelineResult or "
        "SpatialPipelineResult (with a .deconv result object)."
    )


def plot_results(
    result: Any,
    output_dir: Union[str, Path],
    *,
    array_row: Optional[Sequence[float]] = None,
    array_col: Optional[Sequence[float]] = None,
    max_cell_types: int = 6,
) -> list:
    """Produce a standard set of figures for a pipeline result.

    Returns a list of :class:`~tissueresolve.plotting.style.PlotResult`.
    For spatial results, ``array_row`` and ``array_col`` are required (the
    coordinates are not stored on the result object).
    """
    from tissueresolve.results import BulkDeconvResult, SpatialDeconvResult

    deconv = getattr(result, "deconv", None)
    out = Path(output_dir)
    figs: list = []

    if isinstance(deconv, BulkDeconvResult):
        from tissueresolve.plotting import bulk_plots

        figs.append(bulk_plots.composition_barplot(deconv, out))
        figs.append(bulk_plots.proportion_heatmap(deconv, out))
        qc = getattr(result, "qc", None)
        if qc is not None and (qc.recon_r2 is not None or deconv.coverage_r2 is not None):
            figs.append(bulk_plots.qc_barplot(qc, out, deconv=deconv))
        if qc is not None and qc.marker_recall is not None:
            figs.append(bulk_plots.marker_recall_plot(qc, out))
        if deconv.lower_ci is not None:
            figs.append(bulk_plots.bootstrap_ci_plot(deconv, out))
        return figs

    if isinstance(deconv, SpatialDeconvResult):
        if array_row is None or array_col is None:
            raise ValueError(
                "plot_results for spatial results needs array_row and array_col "
                "(spot coordinates are not stored on the result)."
            )
        from tissueresolve.plotting import spatial_plots

        cts = list(deconv.proportions.columns)[:max_cell_types]
        figs.append(spatial_plots.abundance_map_multi(
            deconv, array_row, array_col, cts, out))
        figs.append(spatial_plots.dominant_type_map(
            deconv, array_row, array_col, out))
        morans = getattr(result, "morans_i", None)
        if morans is not None:
            figs.append(spatial_plots.morans_i_barplot(morans, out))
        return figs

    raise TypeError(
        "plot_results: result must be a BulkPipelineResult or "
        "SpatialPipelineResult."
    )
