"""
Spatial deconvolution pipeline for TissueResolve (SpatCAR algorithm).

Orchestration order
-------------------
1.  Input validation and gene name harmonisation
2.  Marker gene selection
3.  Separability diagnostics
4.  Spatial graph construction
5.  NNLS warm start → SpatCAR model fit
6.  (Mismatch correction is internal to the model)
7.  Spot-level QC
8.  Moran's I
9.  Optional neighbourhood statistics
10. SpatialDeconvResult creation

Output type
-----------
``SpatialDeconvResult.proportions`` always records
``ESTIMATE_TYPE = "spot_rna_composition"``.  Spatial smoothing parameters
(``lambda_spatial``, ``n_smooth``) are always recorded and never hidden.
"""
from __future__ import annotations

import logging
import warnings
from dataclasses import dataclass, field
from typing import Any, Optional

import numpy as np
import pandas as pd
import scipy.sparse as sp

from tissueresolve.config import TissueResolveConfig
from tissueresolve.io.validation import check_gene_overlap
from tissueresolve.results import QCReport, ReferenceSignature, SpatialDeconvResult

__all__ = ["SpatialPipeline", "SpatialPipelineResult"]

logger = logging.getLogger("tissueresolve.spatial.pipeline")


# ---------------------------------------------------------------------------
# Result bundle
# ---------------------------------------------------------------------------


@dataclass
class SpatialPipelineResult:
    """All outputs of a completed spatial pipeline run.

    Attributes
    ----------
    deconv:
        :class:`~tissueresolve.results.SpatialDeconvResult`.
    qc:
        :class:`~tissueresolve.results.QCReport` for the spatial run.
    spot_qc:
        Per-spot QC DataFrame (entropy, nb_loglik, spatial_resid, …).
    morans_i:
        Moran's I per cell type.
    neighbourhood:
        :class:`~tissueresolve.spatial.neighbourhood.NeighbourhoodStats`
        or ``None`` when not requested.
    model:
        Fitted :class:`~tissueresolve.spatial.model.SpatCARModel`.
    run_metadata:
        Dict of pipeline parameters.
    """

    deconv: SpatialDeconvResult
    qc: QCReport
    spot_qc: pd.DataFrame
    morans_i: pd.Series
    neighbourhood: Optional[Any] = None
    model: Optional[Any] = None
    run_metadata: dict[str, Any] = field(default_factory=dict)


# ---------------------------------------------------------------------------
# Pipeline
# ---------------------------------------------------------------------------


class SpatialPipeline:
    """Orchestrate the full spatial deconvolution workflow.

    Parameters
    ----------
    config:
        Full :class:`~tissueresolve.config.TissueResolveConfig`.
        When ``None``, defaults are used.
    """

    def __init__(self, config: Optional[TissueResolveConfig] = None) -> None:
        self.cfg = config or TissueResolveConfig()

    def run(
        self,
        Y: np.ndarray,
        ref: ReferenceSignature,
        array_row: np.ndarray,
        array_col: np.ndarray,
        lib_sizes: np.ndarray,
        visium_gene_names: list[str],
        spot_ids: Optional[list[str]] = None,
        *,
        marker_genes: Optional[list[str]] = None,
        run_neighbourhood: bool = False,
        n_neighbourhood_perm: int = 99,
        family_map: Optional[dict] = None,
        rare_protection: Optional[list] = None,
    ) -> SpatialPipelineResult:
        """Run the full spatial deconvolution pipeline.

        Parameters
        ----------
        Y:
            Full gene count matrix for Visium, shape ``(N, G_visium)``
            (dense or sparse).  Rows are spots; columns are genes.
        ref:
            :class:`~tissueresolve.results.ReferenceSignature` (pre-built).
        array_row / array_col:
            Integer Visium array coordinates ``(N,)``.
        lib_sizes:
            Per-spot library sizes ``(N,)`` float32.
        visium_gene_names:
            Gene names corresponding to columns of *Y*.
        spot_ids:
            Optional barcode labels for spots.
        marker_genes:
            Pre-selected marker gene names.  When ``None``, marker selection
            is run automatically via :class:`~tissueresolve.reference.markers.GeneSelector`.
        run_neighbourhood:
            Whether to run neighbourhood co-occurrence statistics.
        n_neighbourhood_perm:
            Permutations for the neighbourhood null (default 99 for speed;
            use 999 for publication).

        Returns
        -------
        SpatialPipelineResult
        """
        from tissueresolve.spatial.graph import build_hex_graph_from_arrays
        from tissueresolve.spatial.model import SpatCARModel
        from tissueresolve.spatial.qc import compute_spot_qc, compute_morans_i, compute_model_qc
        from tissueresolve.reference.separability import compute_separability

        N = len(array_row)
        array_row = np.asarray(array_row, dtype=np.int32)
        array_col = np.asarray(array_col, dtype=np.int32)
        lib_sizes = np.asarray(lib_sizes, dtype=np.float32)

        spot_id_arr = (
            np.array(spot_ids) if spot_ids is not None
            else np.array([f"spot_{i}" for i in range(N)])
        )

        # 1. Gene overlap
        shared_genes = check_gene_overlap(visium_gene_names, ref.gene_names)
        if not shared_genes:
            raise ValueError(
                "No genes overlap between Visium and reference.  "
                "Check gene naming conventions."
            )

        # 2. Marker selection
        if marker_genes is None:
            from tissueresolve.reference.markers import GeneSelector
            selector = GeneSelector(config=self.cfg.genes)
            sel_result = selector.select(ref, query_genes=shared_genes)
            marker_genes = sel_result.selected_genes
            logger.info("Marker selection: %d genes.", len(marker_genes))
        else:
            marker_genes = [g for g in marker_genes if g in set(shared_genes)]
            if not marker_genes:
                raise ValueError("No provided marker genes overlap with shared genes.")

        # 3. Separability diagnostics
        try:
            sep_report = compute_separability(
                ref, warn_threshold=0.90, raise_on_critical=False
            )
            if sep_report.has_problems:
                logger.warning(
                    "Reference has %d HIGH/CRITICAL separability pairs — "
                    "estimates for those types will be unreliable.",
                    sep_report.n_critical + sep_report.n_high,
                )
        except Exception as e:
            logger.debug("Separability check skipped: %s", e)
            sep_report = None

        # 4. Subset Y to marker genes (needed before the graph for edge-aware mode)
        visium_gene_to_col = {g: i for i, g in enumerate(visium_gene_names)}
        marker_col_idx = np.array(
            [visium_gene_to_col[g] for g in marker_genes if g in visium_gene_to_col],
            dtype=np.intp,
        )
        if sp.issparse(Y):
            Y_marker = np.asarray(Y[:, marker_col_idx].todense(), dtype=np.float32)
        else:
            Y_marker = np.asarray(Y[:, marker_col_idx], dtype=np.float32)

        ref_marker = ref.subset_genes(marker_genes)

        # 5. Spatial graph (default: coordinate-only hex graph). Experimental
        # opt-in: an edge-weighted graph so the existing CAR penalty smooths less
        # across likely boundaries. Default behaviour is unchanged.
        cfg_sp = self.cfg.spatial_solver
        if getattr(cfg_sp, "edge_aware", False):
            from tissueresolve.spatial.graph import build_edge_aware_spatial_graph
            graph = build_edge_aware_spatial_graph(
                array_row, array_col, Y_marker, spot_ids=spot_id_arr,
                k_neighbors=cfg_sp.edge_aware_k_neighbors,
                expression_weight=cfg_sp.edge_aware_expression_weight,
                composition_weight=cfg_sp.edge_aware_composition_weight,
                min_edge_weight=cfg_sp.edge_aware_min_weight,
                max_edge_weight=cfg_sp.edge_aware_max_weight,
                batch_size=cfg_sp.n_jobs or 500,
            )
            logger.info("Edge-aware spatial graph: %s", graph.metadata)
        else:
            graph = build_hex_graph_from_arrays(
                array_row, array_col, spot_ids=spot_id_arr,
                batch_size=cfg_sp.n_jobs or 500,
            )

        # 6. Fit model
        cfg_s = self.cfg.spatial_solver
        model = SpatCARModel(
            lambda_spatial=cfg_s.lambda_spatial,
            max_iter=cfg_s.max_iter,
            tol=cfg_s.tol,
            update_mismatch_every=cfg_s.update_mismatch_every,
            n_jobs=cfg_s.n_jobs,
            random_state=cfg_s.random_state,
            verbose=False,
        )
        model.fit(Y_marker, ref_marker, graph, lib_sizes)

        Pi = model.proportions_
        d_g = model._mismatch.d_g if model._mismatch is not None else None
        R_d = (model._R_lin * d_g[np.newaxis, :]).astype(np.float32) if (
            model._R_lin is not None and d_g is not None
        ) else None

        # 6b. Experimental, opt-in state-similarity / sparsity refinement (post-fit,
        # within-family, mass-conserving). Default OFF → Pi unchanged. No solver change.
        state_meta: dict[str, Any] = {}
        cfg_sr = getattr(self.cfg, "state_regularization", None)
        if cfg_sr is not None and getattr(cfg_sr, "enabled", False):
            Pi, state_meta = self._apply_state_regularization(
                Pi, ref_marker, model.cell_types_, cfg_sr, family_map, rare_protection)

        # 7. Spot QC — always pass model arrays so nb_loglik is non-NaN
        spot_qc_df = compute_spot_qc(
            Pi, graph,
            cell_types=model.cell_types_,
            Y_marker=Y_marker,
            R_d=R_d,
            phi_g=model._phi_g,
            lib_sizes=lib_sizes,
        )

        # 8. Moran's I
        morans_i = compute_morans_i(Pi, graph, cell_types=model.cell_types_)

        # 9. Neighbourhood stats (optional)
        neighbourhood = None
        if run_neighbourhood:
            from tissueresolve.spatial.neighbourhood import compute_neighbourhood_stats
            neighbourhood = compute_neighbourhood_stats(
                Pi, graph, model.cell_types_,
                n_permutations=n_neighbourhood_perm,
                seed=cfg_s.random_state,
            )

        # 10. Assemble result
        props_df = pd.DataFrame(
            Pi,
            index=spot_id_arr,
            columns=model.cell_types_,
        )
        props_df.index.name = "spot"

        model_qc = compute_model_qc(
            model.n_iter_, model.converged_, model.convergence_trace_,
            d_g=d_g, marker_genes=model.marker_genes_,
        )

        run_meta: dict[str, Any] = {
            "pipeline": "SpatialPipeline",
            "n_spots": N,
            "n_cell_types": ref.n_cell_types,
            "n_genes_panel": len(marker_genes),
            "lambda_spatial": model.lambda_spatial,
            "alpha": model.alpha,
            "estimate_type": "spot_rna_composition",
            "edge_aware_smoothing_used": bool(getattr(cfg_sp, "edge_aware", False)),
            "state_regularization_used": bool(state_meta.get("enabled", False)),
        }
        if state_meta:
            for k, v in state_meta.items():
                run_meta[f"state_reg_{k}"] = v
        if getattr(cfg_sp, "edge_aware", False):
            gm = getattr(graph, "metadata", {}) or {}
            for key in ("edge_weight_min", "edge_weight_max", "edge_weight_mean",
                        "edge_weight_median", "k_neighbors", "used_initial_proportions",
                        "fallback", "expression_weight", "composition_weight"):
                if key in gm:
                    run_meta[f"edge_{key}" if not key.startswith("edge_") else key] = gm[key]
        run_meta.update(model_qc)

        deconv = SpatialDeconvResult(
            proportions=props_df,
            cell_types=model.cell_types_,
            marker_genes=model.marker_genes_,
            n_iter=model.n_iter_,
            converged=model.converged_,
            convergence_trace=model.convergence_trace_,
            lambda_spatial=model.lambda_spatial,
            mismatch_factors=d_g,
            run_metadata=run_meta,
        )

        qc_report = QCReport(
            modality="spatial",
            recommendations=[],
            spot_qc=spot_qc_df,
            morans_i=morans_i,
            model_qc=model_qc,
        )

        return SpatialPipelineResult(
            deconv=deconv,
            qc=qc_report,
            spot_qc=spot_qc_df,
            morans_i=morans_i,
            neighbourhood=neighbourhood,
            model=model,
            run_metadata=run_meta,
        )

    # ------------------------------------------------------------------
    # Experimental, opt-in post-fit state-similarity refinement
    # ------------------------------------------------------------------

    def _apply_state_regularization(
        self, Pi, ref_marker, cell_types, cfg_sr, family_map, rare_protection,
    ):
        """Apply the experimental within-family state refinement to ``Pi``.

        Default behaviour is never reached (gated on ``cfg_sr.enabled``). Returns
        ``(Pi_refined, metadata)``. On any failure it falls back to the input
        ``Pi`` and records the reason — it must never break a spatial run.
        """
        from tissueresolve.experimental.state_similarity_regularization import (
            compute_state_similarity_graph, apply_state_regularized_refinement,
            recommend_state_groups,
        )

        meta: dict[str, Any] = {
            "enabled": True, "mode": cfg_sr.mode,
            "lambda_state": float(cfg_sr.lambda_state),
            "lambda_sparse": float(cfg_sr.lambda_sparse),
            "within_family_only": bool(cfg_sr.within_family_only),
            "preserve_broad_mass": bool(cfg_sr.preserve_broad_mass),
            "family_map_provided": family_map is not None,
            "feature_status": "experimental",
        }
        if cfg_sr.within_family_only and family_map is None:
            warnings.warn(
                "state_regularization enabled with within_family_only=True but no "
                "family_map was provided; the state graph has no within-family edges "
                "so the refinement is a no-op.  Pass family_map={fine: broad} to "
                "deconv_spatial to enable within-family refinement.", stacklevel=2)

        try:
            profiles = ref_marker.as_R_cpm()  # (K, G_m)
            graph = compute_state_similarity_graph(
                profiles, list(cell_types), family_map=family_map,
                method=cfg_sr.similarity_method, k_states=cfg_sr.k_states,
                min_similarity=cfg_sr.min_similarity,
                within_family_only=cfg_sr.within_family_only)
            meta["state_graph_n_edges"] = graph.n_edges
            meta["state_graph_n_components"] = graph.metadata["n_connected_components"]

            df = pd.DataFrame(Pi, columns=list(cell_types))

            if cfg_sr.mode == "adaptive_resolution":
                rec = recommend_state_groups(
                    graph, family_map=family_map, threshold=cfg_sr.group_threshold,
                    rare_protection=rare_protection)
                meta["adaptive_recommended_groups"] = rec["recommended_groups"]
                meta["adaptive_n_groups"] = rec["n_groups"]
                meta["estimates_modified"] = False
                return Pi, meta  # diagnostic only — estimates unchanged

            lam_state = 0.0 if cfg_sr.mode == "sparsity" else cfg_sr.lambda_state
            refined, rmeta = apply_state_regularized_refinement(
                df, state_graph=graph, family_map=family_map,
                lambda_state=lam_state, lambda_sparse=cfg_sr.lambda_sparse,
                preserve_broad_mass=cfg_sr.preserve_broad_mass,
                rare_protection=rare_protection)
            meta.update({
                "n_families_refined": rmeta["n_families_refined"],
                "mass_shifted_state": rmeta["mass_shifted_state"],
                "mass_shifted_sparse": rmeta["mass_shifted_sparse"],
                "broad_mass_max_deviation": rmeta["broad_mass_max_deviation"],
                "estimates_modified": True,
            })
            return refined.to_numpy(dtype=np.float32), meta
        except Exception as exc:  # noqa: BLE001 — never break a run
            warnings.warn(
                f"state_regularization refinement failed ({exc}); returning "
                "un-refined spatial estimates.", stacklevel=2)
            meta["error"] = str(exc)
            meta["estimates_modified"] = False
            return Pi, meta
