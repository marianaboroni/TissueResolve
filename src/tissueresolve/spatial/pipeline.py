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

        # 4. Spatial graph
        graph = build_hex_graph_from_arrays(
            array_row, array_col, spot_ids=spot_id_arr,
            batch_size=self.cfg.spatial_solver.n_jobs or 500,
        )

        # 5. Subset Y to marker genes
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
        }
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
