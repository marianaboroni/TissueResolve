"""
Bulk deconvolution pipeline (CHIMERA algorithm).

Orchestration order
-------------------
1.  Input validation (gene names, sample names, counts matrix)
2.  Gene overlap check
3.  Protocol risk assessment (optional)
4.  Marker gene selection
5.  Gene weight computation
6.  Weighted NNLS solve
7.  Bootstrap CI (optional)
8.  mRNA content correction (optional)
9.  Bulk QC

The pipeline does not implement HTML reports or full publication plotting.

Output type
-----------
``BulkDeconvResult.proportions`` always records ``ESTIMATE_TYPE = "mRNA_proportion"``.
mRNA content correction, when applied, populates ``BulkDeconvResult.cell_fractions``
as a separate field and never overwrites ``proportions``.
"""
from __future__ import annotations

import logging
import warnings
from dataclasses import dataclass, field
from typing import Any, Optional

import numpy as np
import pandas as pd

from tissueresolve.config import BootstrapConfig, BulkQCConfig, GeneConfig, TissueResolveConfig
from tissueresolve.io.validation import check_gene_overlap, validate_counts_matrix
from tissueresolve.results import BulkDeconvResult, QCReport, ReferenceSignature

__all__ = ["BulkPipeline", "BulkPipelineResult"]

logger = logging.getLogger("tissueresolve.bulk.pipeline")


# ---------------------------------------------------------------------------
# Result bundle
# ---------------------------------------------------------------------------


@dataclass
class BulkPipelineResult:
    """All outputs of a completed bulk pipeline run.

    Attributes
    ----------
    deconv:
        :class:`~tissueresolve.results.BulkDeconvResult` with mRNA proportions
        (and optionally bootstrap CIs and cell fractions).
    qc:
        :class:`~tissueresolve.results.QCReport` for the bulk run.
    gene_selection:
        :class:`~tissueresolve.reference.markers.MarkerSelectionResult`.
        ``None`` if marker selection was skipped.
    protocol_risk:
        :class:`~tissueresolve.protocol.risk.ProtocolRiskReport`.
        ``None`` if protocol assessment was skipped.
    run_metadata:
        Dict of all pipeline parameters and version info.
    """

    deconv: BulkDeconvResult
    qc: QCReport
    gene_selection: Optional[Any] = None
    protocol_risk: Optional[Any] = None
    run_metadata: dict[str, Any] = field(default_factory=dict)
    #: Optional calibrated identifiability certificate (experimental, opt-in via
    #: ``deconv_bulk(..., identifiability=True)``).  A reported diagnostic only — it
    #: does NOT modify ``deconv.proportions`` or any estimate.  ``None`` unless requested.
    identifiability: Optional[Any] = None


# ---------------------------------------------------------------------------
# Pipeline
# ---------------------------------------------------------------------------


class BulkPipeline:
    """Orchestrate the full bulk RNA-seq deconvolution workflow.

    Parameters
    ----------
    config:
        Full :class:`~tissueresolve.config.TissueResolveConfig`.
        When ``None``, defaults are used.

    Examples
    --------
    Basic usage::

        from tissueresolve.bulk.pipeline import BulkPipeline
        result = BulkPipeline().run(bulk_counts, ref)

    With bootstrap::

        result = BulkPipeline().run(bulk_counts, ref, n_bootstrap=200)

    With mRNA correction::

        from tissueresolve.bulk.solver import MRNAContentCorrector
        corrector = MRNAContentCorrector(mrna_content_series)
        result = BulkPipeline().run(bulk_counts, ref, mrna_corrector=corrector)
    """

    def __init__(self, config: Optional[TissueResolveConfig] = None) -> None:
        self.cfg = config or TissueResolveConfig()

    # ------------------------------------------------------------------
    # Public API
    # ------------------------------------------------------------------

    def run(
        self,
        bulk: pd.DataFrame,
        ref: ReferenceSignature,
        *,
        gene_panel: Optional[list[str]] = None,
        gene_weights: Optional[pd.Series] = None,
        protocol_meta: Optional[Any] = None,
        mrna_corrector: Optional[Any] = None,
        n_bootstrap: Optional[int] = None,
        run_qc: bool = True,
    ) -> BulkPipelineResult:
        """Run the full bulk deconvolution pipeline.

        Parameters
        ----------
        bulk:
            Genes × samples count matrix.  Index is gene names; columns are
            sample identifiers.
        ref:
            :class:`~tissueresolve.results.ReferenceSignature` (pre-built).
        gene_panel:
            Pre-selected gene panel.  When ``None``, marker selection is run
            automatically via :class:`~tissueresolve.reference.markers.GeneSelector`.
        gene_weights:
            Pre-computed per-gene weights.  When ``None`` and *gene_panel* is
            also ``None``, weights are computed by :class:`~GeneSelector`.
        protocol_meta:
            :class:`~tissueresolve.protocol.metadata.ProtocolMetadata`.
            When provided, protocol risk assessment is run before marker
            selection.
        mrna_corrector:
            :class:`~tissueresolve.bulk.solver.MRNAContentCorrector`.
            When provided, ``BulkDeconvResult.cell_fractions`` is populated.
            **Does not modify ``proportions``.**
        n_bootstrap:
            Number of bootstrap iterations.  Overrides ``config.bootstrap.n_bootstrap``.
            ``0`` → skip bootstrap.
        run_qc:
            Whether to compute QC metrics (default ``True``).

        Returns
        -------
        BulkPipelineResult
        """
        # 1. Validate inputs
        validate_counts_matrix(bulk)

        # 2. Gene overlap
        shared_genes = check_gene_overlap(list(bulk.index), ref.gene_names)
        if not shared_genes:
            raise ValueError(
                "No genes overlap between bulk matrix and reference.  "
                "Ensure gene identifiers use the same naming convention."
            )

        # 3. Protocol risk assessment
        protocol_risk = None
        protocol_risk_scores: Optional[pd.Series] = None
        n_genes_excluded_proto = 0
        if protocol_meta is not None:
            from tissueresolve.protocol.risk import ProtocolRiskAssessor
            assessor = ProtocolRiskAssessor(genome=ref.genome)
            protocol_risk = assessor.assess(protocol_meta, shared_genes)
            protocol_risk_scores = assessor.gene_risk_scores(protocol_meta, shared_genes)
            n_genes_excluded_proto = protocol_risk.n_genes_excluded
            if protocol_risk.excluded_genes:
                logger.info(
                    "Protocol filter: excluding %d genes (risk_level=%s).",
                    n_genes_excluded_proto, protocol_risk.risk_level,
                )

        # 4. Marker selection (when no pre-selected panel)
        gene_selection = None
        # Opt-in: a reference built with a stored panel (e.g. donor-aware DE selection
        # at build time) uses it as the default panel. Backward-compatible: references
        # without `selected_genes` fall through to the standard marker selection.
        if gene_panel is None and getattr(ref, "selected_genes", None):
            stored = [g for g in ref.selected_genes if g in set(shared_genes)]
            if len(stored) >= 2:
                gene_panel = stored
                logger.info("Using reference-stored gene panel: %d genes "
                            "(of %d stored, in shared).", len(stored),
                            len(ref.selected_genes))
        if gene_panel is None:
            from tissueresolve.reference.markers import GeneSelector
            selector = GeneSelector(config=self.cfg.genes)
            bulk_mean = bulk.reindex(shared_genes).fillna(0.0).mean(axis=1)
            gene_selection = selector.select(
                ref,
                query_genes=shared_genes,
                protocol_risk_scores=protocol_risk_scores,
                bulk_mean=bulk_mean,
            )
            gene_panel = gene_selection.selected_genes
            logger.info(
                "Marker selection: %d genes selected.", len(gene_panel)
            )

        # 5. Gene weights (when not pre-computed)
        if gene_weights is None and gene_selection is not None:
            from tissueresolve.reference.markers import GeneSelector
            selector = GeneSelector(config=self.cfg.genes)
            bulk_mean_w = bulk.reindex(shared_genes).fillna(0.0).mean(axis=1)
            gene_weights = selector.compute_gene_weights(
                ref,
                panel=gene_panel,
                protocol_risk_scores=protocol_risk_scores,
                bulk_mean=bulk_mean_w,
                bulk_full=bulk.reindex(gene_panel).fillna(0.0),
            )

        # 6. Solve. Default = weighted NNLS (unchanged). Experimental, opt-in:
        # count-likelihood (Poisson/NB) GLM solver selected via cfg.bulk_solver.method.
        method = getattr(self.cfg.bulk_solver, "method", "wNNLS")
        if method in ("poisson_glm_experimental", "nb_glm_experimental"):
            from tissueresolve.experimental.nb_bulk_solver import NBGLMBulkSolver
            loss = "poisson" if method == "poisson_glm_experimental" else "nb"
            solver = NBGLMBulkSolver(loss=loss, seed=self.cfg.bulk_solver.seed,
                                     max_iter=self.cfg.bulk_solver.max_iter,
                                     tol=self.cfg.bulk_solver.tol)
            deconv = solver.solve(bulk, ref, gene_panel, gene_weights)
        elif method == "wNNLS":
            from tissueresolve.bulk.solver import WNNLSSolver
            solver = WNNLSSolver(seed=self.cfg.bulk_solver.seed)
            deconv = solver.solve(bulk, ref, gene_panel, gene_weights)
        else:
            raise ValueError(
                f"Unknown bulk_solver.method {method!r}. Use 'wNNLS' (default), "
                "'poisson_glm_experimental', or 'nb_glm_experimental'.")

        # Update metadata
        deconv.run_metadata.update({
            "n_genes_excluded_protocol": n_genes_excluded_proto,
            "protocol_risk_level": protocol_risk.risk_level if protocol_risk else "unknown",
        })

        # 7. Bootstrap CI
        n_boot = n_bootstrap if n_bootstrap is not None else self.cfg.bootstrap.n_bootstrap
        if n_boot > 0:
            from tissueresolve.uncertainty.bootstrap import BulkBootstrapCI
            boot_cfg = BootstrapConfig(
                n_bootstrap=n_boot,
                bootstrap_frac=self.cfg.bootstrap.bootstrap_frac,
                ci_level=self.cfg.bootstrap.ci_level,
                seed=self.cfg.bootstrap.seed,
            )
            bootstrapper = BulkBootstrapCI(config=boot_cfg)
            bootstrapper.attach(deconv, bulk, ref)
            deconv.run_metadata.update(bootstrapper.ci_metadata)
            # Bootstrap resampling uses wNNLS regardless of the point solver; record
            # this so GLM-point + bootstrap-CI provenance is never silently mislabelled.
            if method != "wNNLS":
                deconv.run_metadata["bootstrap_solver"] = "wNNLS"
                deconv.run_metadata["point_solver"] = method

        # 8. mRNA content correction
        if mrna_corrector is not None:
            cell_fracs = mrna_corrector.correct(deconv.proportions)
            deconv.cell_fractions = cell_fracs

        # 9. QC
        qc_report: QCReport
        if run_qc:
            from tissueresolve.bulk.qc import BulkQC
            qc = BulkQC(config=self.cfg.bulk_qc)
            qc_report = qc.compute(deconv, ref, bulk, protocol_risk_report=protocol_risk)
        else:
            qc_report = QCReport(modality="bulk")

        run_metadata: dict[str, Any] = {
            "pipeline": "BulkPipeline",
            "n_samples": deconv.n_samples,
            "n_cell_types": deconv.n_cell_types,
            "n_genes_panel": len(gene_panel),
            "has_bootstrap": deconv.lower_ci is not None,
            "has_mrna_correction": deconv.cell_fractions is not None,
        }
        run_metadata.update(deconv.run_metadata)

        return BulkPipelineResult(
            deconv=deconv,
            qc=qc_report,
            gene_selection=gene_selection,
            protocol_risk=protocol_risk,
            run_metadata=run_metadata,
        )
