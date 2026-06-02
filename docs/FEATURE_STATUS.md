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
| Hyperparameter tuning | **not implemented** | `tuning/` is a quarantined stub that raises `NotImplementedError`; deferred to v0.2 |
| Full BayesPrism-like Bayesian model (Gibbs) | **not implemented** | deferred |

## Required wordings (use verbatim in README/CLI/report)

- **State-aware:** "Experimental: state-aware deconvolution runs behind
  `--state-aware`. It is not part of the default v0.1 workflow and has not been
  validated across real datasets."
- **State-aware report integration:** "State-aware outputs are experimental; full
  standard-report integration is limited."
- **Expression reconstruction:** "Cell-type-specific expression reconstruction is
  planned/deferred and not implemented in v0.1." Use "cell-type-level
  deconvolution estimates" / "RNA-derived composition estimates" instead.
- **External benchmark:** "Only executed or imported tools are ranked.
  Exported-only or skipped tools are not benchmarked."
- **Spatial benchmark:** "The real Visium benchmark does not measure accuracy
  without ground truth. It evaluates concordance, spatial structure, marker
  agreement, runtime, and output completeness."
- **Estimate type:** "Bulk estimates are RNA-derived proportions, not absolute
  cell fractions. Spatial estimates are spot-level RNA-derived composition, not
  single-cell labels."
- **Publication status:** "TissueResolve is currently alpha / early-access
  research software. It is not yet a fully publication-ready method until
  external benchmarks, multi-dataset validation, and final API stabilization are
  complete."

Deferred features (**reference adaptation**, **expression reconstruction**,
**hyperparameter tuning**, **full Bayesian model**) are intentionally absent and
must not be presented as available in the README, CLI help, or report.
