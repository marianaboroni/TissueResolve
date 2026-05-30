"""
Bulk RNA-seq deconvolution workflow (CHIMERA algorithm).

This subpackage contains all bulk-specific algorithmic components.
It must not import from ``spatial`` or mix bulk and spatial model assumptions.

Public API
----------
``solver.WNNLSSolver``
    Weighted NNLS solver per bulk sample.  Returns ``BulkDeconvResult``.

``solver.MRNAContentCorrector``
    Converts mRNA proportions to cell fractions.  Requires external per-cell-type
    mRNA content data.  Never run silently with trivial (uniform) content.

``qc.BulkQC``
    Computes per-sample R², RMSE, profile correlation, mismatch flags;
    per-cell-type marker recall, spillover risk, condition number.

``pipeline.BulkPipeline``
    Orchestrates: input validation → gene matching → protocol risk →
    marker selection → NNLS → optional bootstrap → optional mRNA correction
    → QC → ``BulkPipelineResult``.
"""
from tissueresolve.bulk.solver import MRNAContentCorrector, WNNLSSolver
from tissueresolve.bulk.qc import BulkQC
from tissueresolve.bulk.pipeline import BulkPipeline, BulkPipelineResult

__all__ = [
    "WNNLSSolver",
    "MRNAContentCorrector",
    "BulkQC",
    "BulkPipeline",
    "BulkPipelineResult",
]
