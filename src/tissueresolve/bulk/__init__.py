"""
Bulk RNA-seq deconvolution workflow (CHIMERA algorithm).

This subpackage contains all bulk-specific algorithmic components.
It must not import from ``spatial`` or mix bulk and spatial model assumptions.

Submodules (implemented in Stage 3):

``solver``
    WNNLSSolver — weighted NNLS per bulk sample.
    MRNAContentCorrector — converts mRNA proportions to cell fractions.
    Both are part of the public API.

``pipeline``
    ChimeraPipeline — orchestrates: reference construction →
    protocol risk assessment → gene selection → optional discordance
    filter → NNLS solve → bootstrap CI → QC → report.

``qc``
    BulkQCReport — per-sample R², profile correlation, mismatch flag;
    per-cell-type marker recall, spillover risk, condition number.

``benchmark``
    Five benchmark scenarios: SANITY_PSEUDOBULK, TISSUE_MISMATCH_CONTEXT,
    PROTOCOL_MISMATCH_CENTRAL, EXTERNAL_REFERENCE, PAIRED_CALIBRATION.
"""
