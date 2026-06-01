# TissueResolve feature status

Honest classification of each feature. **stable** = validated + default-safe;
**experimental** = usable behind an explicit flag, not validated across datasets;
**partial** = implemented but not fully wired/measured; **planned** = designed,
not implemented; **not implemented** = absent (and must not be claimed).

| Feature | Status | Notes |
|---|---|---|
| Bulk deconvolution (wNNLS / pipeline) | **stable** | default path; broad test coverage |
| Spatial deconvolution (NB-CAR) | **stable** | default path |
| Solver `auto` (gene-masking CV) | **stable** | matches NNLS on real bulk |
| Hierarchical broad→fine | **stable** | default `resolution-mode auto/hierarchical`; unresolved mass at family level |
| Granular signatures (broad / cell-type / state panels) | **experimental** | `reference/granular_signatures.py`, tested; built from per-cell AnnData at build time |
| State-aware hierarchy (broad→cell type→state) | **experimental** | `bulk/state_aware_hierarchical.py`; behind `--state-aware` / `deconv_bulk(state_aware=True)`; **not default**; falls back to 2-level when no state labels; **not yet validated on real data** |
| Reference suitability score | **stable** | PASS/CAUTION/WARNING/FAIL with worst-component override |
| Report generation (QC-first, embedded figures) | **stable** | unified report; figure manifest |
| Spatial benchmark (concordance/structure, no-truth) | **partial** | bulk/spatial split done; multi-metric ranking module added; marker-recovery/synthetic-truth/dashboard not yet wired into the run |
| External-tool benchmark (MuSiC/CARD/cell2location/…) | **partial** | several tools executed; others skipped/failed honestly |
| Reference adaptation (residual-driven gene-weight update) | **not implemented** | deferred; do not enable/claim |
| Cell-type-specific expression reconstruction | **not implemented** | deferred; do not enable/claim |

Deferred features (**reference adaptation**, **expression reconstruction**) are
intentionally absent in this step and must not be presented as available in the
README, CLI help, or report.
