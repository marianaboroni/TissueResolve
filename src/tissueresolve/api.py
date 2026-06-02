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


def _solver_bulk_result(bulk, ref, solver: str, cfg):
    """Run a chosen solver backbone (``nnls``/``ridge_nnls``/``auto``/…) and wrap
    the flat estimate in a ``BulkPipelineResult``-compatible object.

    This is the ``--solver`` path; it is additive and only used when a solver is
    explicitly requested.  The default ``deconv_bulk`` path (protocol-aware
    weighted pipeline) is unchanged."""
    import numpy as np
    from tissueresolve.solver import get_solver
    from tissueresolve.bulk.pipeline import BulkPipelineResult
    from tissueresolve.results import BulkDeconvResult, QCReport

    res = get_solver(solver, n_splits=2).solve(bulk, ref) if solver in ("auto", "ensemble_nnls") \
        else get_solver(solver).solve(bulk, ref)
    props = res.proportions
    sub = ref.subset_genes(res.genes_used)
    R = sub.as_R_cpm()                       # K × g
    cols = [c for c in props.columns if c in set(sub.cell_types)]
    idx = [list(sub.cell_types).index(c) for c in cols]
    recon = props[cols].to_numpy(float) @ R[idx]   # samples × g
    bsub = bulk.copy(); bsub.index = bsub.index.map(str)
    obs = bsub.loc[list(sub.gene_names)].to_numpy(float).T  # samples × g
    r2 = []
    for i in range(obs.shape[0]):
        o, p = obs[i], recon[i]
        ss_res = float(np.sum((o - p) ** 2)); ss_tot = float(np.sum((o - o.mean()) ** 2))
        r2.append(1.0 - ss_res / ss_tot if ss_tot > 0 else float("nan"))
    import pandas as pd
    deconv = BulkDeconvResult(
        proportions=props, coverage_r2=pd.Series(r2, index=props.index),
        gene_panel=list(res.genes_used),
        run_metadata={"solver": solver, **res.diagnostics, "resolution_mode": "flat"})
    qc = QCReport(modality="bulk",
                  recommendations=[f"solver={solver}: "
                                   f"{res.diagnostics.get('selection_reason', '')}".strip()])
    return BulkPipelineResult(deconv=deconv, qc=qc, gene_selection=None,
                              protocol_risk=None,
                              run_metadata={"solver": solver, **res.diagnostics})


def deconv_bulk(
    bulk, ref, *, config=None, resolution_mode: Optional[str] = None,
    hierarchy_mapping: Optional[dict] = None, solver: Optional[str] = None,
    state_aware: bool = False, reference_adata=None,
    state_to_celltype: Optional[dict] = None,
    broad_col: Optional[str] = None, cell_type_col: Optional[str] = None,
    state_col: Optional[str] = None,
    **kwargs,
):
    """Run the bulk deconvolution pipeline.  Returns ``BulkPipelineResult``.

    Wrapper over :meth:`BulkPipeline.run`; extra keyword arguments
    (``gene_panel``, ``n_bootstrap``, ``protocol_meta``, ``mrna_corrector`` …)
    are forwarded verbatim.

    Parameters
    ----------
    resolution_mode:
        ``"none"`` or ``"suggest"`` for flat results, ``"auto"`` to merge
        non-separable types before deconvolution, or ``"hierarchical"`` to run
        broad-family deconvolution first then refine via within-family subtype
        estimates (returning a
        :class:`~tissueresolve.bulk.hierarchical.HierarchicalBulkResult`).
    hierarchy_mapping:
        For ``"hierarchical"`` mode, an explicit ``{fine_cell_type:
        broad_family}`` mapping.  When ``None``, families are inferred from the
        cell-type labels (a warning is emitted) — providing explicit broad/fine
        annotations or a mapping file is strongly preferred.
    """
    from tissueresolve.config import TissueResolveConfig
    from tissueresolve.bulk.pipeline import BulkPipeline

    import warnings as _warnings

    cfg = config or TissueResolveConfig()
    # Explicit solver backbone (additive): only when requested.  Hierarchical +
    # solver combination is not yet wired, so a solver implies a flat estimate.
    if solver is not None and resolution_mode in (None, "flat", "none", "auto"):
        return _solver_bulk_result(bulk, ref, solver, cfg)
    if hasattr(cfg, "resolution") and resolution_mode is None:
        resolution_mode = getattr(cfg.resolution, "resolution_mode", None)
    resolution_mode = resolution_mode or "auto"
    if resolution_mode == "flat":
        resolution_mode = "none"
    if resolution_mode not in ("none", "suggest", "auto", "hierarchical"):
        raise ValueError(
            "resolution_mode must be one of: auto, hierarchical, flat, none, suggest"
        )

    # auto: prefer hierarchical broad→fine when a hierarchy mapping is available,
    # otherwise fall back to flat (fine-only) with a caution.
    if resolution_mode == "auto":
        if hierarchy_mapping is not None:
            resolution_mode = "hierarchical"
        else:
            _warnings.warn(
                "deconv_bulk(resolution_mode='auto'): no hierarchy_mapping was "
                "provided, so flat (fine-only) deconvolution is used.  Provide a "
                "fine→broad mapping to enable the recommended hierarchical "
                "broad→fine workflow.", stacklevel=2)
            result = BulkPipeline(cfg).run(bulk, ref, **kwargs)
            result.deconv.run_metadata["resolution_mode"] = "none"
            result.deconv.run_metadata["resolution_mode_requested"] = "auto"
            return result

    if resolution_mode in ("none", "suggest"):
        result = BulkPipeline(cfg).run(bulk, ref, **kwargs)
        result.deconv.run_metadata["resolution_mode"] = resolution_mode
        return result

    # ---- hierarchical mode ----
    from tissueresolve.bulk.hierarchical import run_hierarchical_bulk
    from tissueresolve.reference.hierarchy import build_cell_type_hierarchy

    hcfg = getattr(cfg, "hierarchical", None)
    gate = dict(
        allow_unresolved=getattr(hcfg, "allow_unresolved", True),
        unresolved_threshold=getattr(hcfg, "unresolved_threshold", 0.10),
        min_discriminating_genes=getattr(hcfg, "min_discriminating_genes", 10),
        within_family_spillover_threshold=getattr(
            hcfg, "within_family_spillover_threshold", 0.30),
        allow_partial_resolution=getattr(hcfg, "allow_partial_resolution", True),
        subtype_confidence_threshold=getattr(hcfg, "subtype_confidence_threshold", 0.10),
    )
    # Explicit gating kwargs override the cfg-derived defaults (avoids a
    # duplicate-keyword collision when callers pass e.g. min_discriminating_genes).
    for _k in list(kwargs):
        if _k in gate:
            gate[_k] = kwargs.pop(_k)

    # ---- experimental state-aware (broad → cell type → state) mode ----
    # Routed BEFORE building the fine→broad mapping: with state labels the
    # reference's cell_types are *states*, so the cell_type→broad mapping comes
    # from the raw hierarchy_mapping, not from build_cell_type_hierarchy(states).
    if state_aware:
        return _run_state_aware_bulk(
            bulk, ref, hierarchy_mapping, cfg, gate,
            reference_adata=reference_adata, state_to_celltype=state_to_celltype,
            broad_col=broad_col, cell_type_col=cell_type_col, state_col=state_col,
            **kwargs)

    mapping = build_cell_type_hierarchy(list(ref.cell_types), hierarchy_mapping)
    return run_hierarchical_bulk(bulk, ref, mapping, config=cfg, **gate, **kwargs)


def _run_state_aware_bulk(bulk, ref, hierarchy_mapping, cfg, gate, *,
                          reference_adata=None, state_to_celltype=None,
                          broad_col=None, cell_type_col=None, state_col=None,
                          **kwargs):
    """Route to the experimental state-aware 3-level solver (default off).

    When state labels exist, the ``cell_type → broad`` mapping is the raw
    *hierarchy_mapping* (the reference's ``cell_types`` are then *states*, so the
    broad mapping cannot be derived from them).  When no state labels exist, runs
    a two-level fallback (broad → cell type) built from the reference cell types
    and records ``fallback_reason``.  Granular gene panels are built only when a
    ``reference_adata`` (per-cell) is provided.
    """
    import warnings as _warnings
    from tissueresolve.bulk.state_aware_hierarchical import (
        run_state_aware_hierarchical_bulk)
    from tissueresolve.reference.hierarchy import build_cell_type_hierarchy
    from tissueresolve.reference.three_level_hierarchy import (
        build_three_level_hierarchy, build_three_level_from_two_level)

    fallback_reason = None
    if state_to_celltype:
        if not hierarchy_mapping:
            raise ValueError(
                "state_aware=True with state_to_celltype requires a "
                "hierarchy_mapping (cell_type→broad).")
        hierarchy = build_three_level_hierarchy(state_to_celltype,
                                                dict(hierarchy_mapping))
    else:
        mapping = build_cell_type_hierarchy(list(ref.cell_types), hierarchy_mapping)
        hierarchy = build_three_level_from_two_level(mapping)
        fallback_reason = ("no state labels available; running broad→cell_type "
                           "two-level fallback (no third-level states)")
        _warnings.warn("deconv_bulk(state_aware=True): " + fallback_reason,
                       stacklevel=2)

    celltype_panels = state_panels = None
    panel_source = "global_genes"
    if reference_adata is not None and broad_col and cell_type_col:
        from tissueresolve.reference.granular_signatures import (
            build_multigranularity_panels)
        with _warnings.catch_warnings():
            _warnings.simplefilter("ignore")
            panels = build_multigranularity_panels(
                reference_adata, broad_col, cell_type_col, state_col)
        celltype_panels = panels["celltype_panels"].family_panels
        state_panels = (panels["state_panels"].family_panels
                        if panels["state_panels"] is not None else None)
        panel_source = "granular_within_family"

    result = run_state_aware_hierarchical_bulk(
        bulk, ref, hierarchy, config=cfg,
        celltype_panels=celltype_panels, state_panels=state_panels, **gate, **kwargs)
    result.metadata.update({
        "hierarchy_mode": "state_aware",
        "state_aware_enabled": True,
        "feature_status": "experimental",
        "broad_col": broad_col, "cell_type_col": cell_type_col,
        "state_col": state_col,
        "has_states": hierarchy.has_states,
        "panel_source": panel_source,
        "fallback_reason": fallback_reason,
    })
    return result


def deconv_spatial(
    Y, ref, array_row, array_col, lib_sizes, gene_names,
    spot_ids=None, *, config=None, resolution_mode: Optional[str] = None,
    hierarchy_mapping: Optional[dict] = None, **kwargs,
):
    """Run the spatial deconvolution pipeline.  Returns ``SpatialPipelineResult``.

    Wrapper over :meth:`SpatialPipeline.run`; extra keyword arguments
    (``marker_genes``, ``run_neighbourhood`` …) are forwarded verbatim.

    Parameters
    ----------
    resolution_mode:
        ``"none"`` or ``"suggest"`` for flat results, ``"auto"`` to merge
        non-separable types before deconvolution, or ``"hierarchical"`` to run
        broad-family deconvolution first then refine via within-family subtype
        estimates (returning a
        :class:`~tissueresolve.spatial.hierarchical.HierarchicalSpatialResult`).
    hierarchy_mapping:
        For ``"hierarchical"`` mode, an explicit ``{fine_cell_type:
        broad_family}`` mapping.  When ``None``, families are inferred from the
        labels (a warning is emitted).
    """
    from tissueresolve.config import TissueResolveConfig
    from tissueresolve.spatial.pipeline import SpatialPipeline

    import warnings as _warnings

    cfg = config or TissueResolveConfig()
    if hasattr(cfg, "resolution") and resolution_mode is None:
        resolution_mode = getattr(cfg.resolution, "resolution_mode", None)
    resolution_mode = resolution_mode or "auto"
    if resolution_mode == "flat":
        resolution_mode = "none"
    if resolution_mode not in ("none", "suggest", "auto", "hierarchical"):
        raise ValueError(
            "resolution_mode must be one of: auto, hierarchical, flat, none, suggest"
        )

    if resolution_mode == "auto":
        if hierarchy_mapping is not None:
            resolution_mode = "hierarchical"
        else:
            _warnings.warn(
                "deconv_spatial(resolution_mode='auto'): no hierarchy_mapping was "
                "provided, so flat (fine-only) deconvolution is used.  Provide a "
                "fine→broad mapping to enable the recommended hierarchical "
                "broad→fine workflow.", stacklevel=2)
            result = SpatialPipeline(cfg).run(
                Y, ref, array_row, array_col, lib_sizes, gene_names, spot_ids,
                **kwargs)
            result.deconv.run_metadata["resolution_mode"] = "none"
            result.deconv.run_metadata["resolution_mode_requested"] = "auto"
            return result

    if resolution_mode in ("none", "suggest"):
        result = SpatialPipeline(cfg).run(
            Y, ref, array_row, array_col, lib_sizes, gene_names, spot_ids, **kwargs
        )
        result.deconv.run_metadata["resolution_mode"] = resolution_mode
        return result

    # ---- hierarchical mode ----
    from tissueresolve.spatial.hierarchical import run_hierarchical_spatial
    from tissueresolve.reference.hierarchy import build_cell_type_hierarchy

    mapping = build_cell_type_hierarchy(list(ref.cell_types), hierarchy_mapping)
    hcfg = getattr(cfg, "hierarchical", None)
    gate = dict(
        allow_unresolved=getattr(hcfg, "allow_unresolved", True),
        unresolved_threshold=getattr(hcfg, "unresolved_threshold", 0.10),
        min_discriminating_genes=getattr(hcfg, "min_discriminating_genes", 10),
        within_family_spillover_threshold=getattr(
            hcfg, "within_family_spillover_threshold", 0.30),
        allow_partial_resolution=getattr(hcfg, "allow_partial_resolution", True),
        subtype_confidence_threshold=getattr(hcfg, "subtype_confidence_threshold", 0.10),
    )
    return run_hierarchical_spatial(
        Y, ref, array_row, array_col, lib_sizes, gene_names, mapping,
        spot_ids, config=cfg, **gate, **kwargs)


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
    from tissueresolve.report import orchestration
    from tissueresolve.results import BulkDeconvResult, SpatialDeconvResult

    deconv = getattr(result, "deconv", None)
    if isinstance(deconv, BulkDeconvResult):
        return orchestration.generate_bulk_report(result, output_path, **kwargs)
    if isinstance(deconv, SpatialDeconvResult):
        return orchestration.generate_spatial_report(result, output_path, **kwargs)
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
