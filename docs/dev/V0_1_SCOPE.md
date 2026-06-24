# TissueResolve v0.1 scope

> **Status: TissueResolve is currently alpha / early-access research software.
> It is not yet a fully publication-ready method until external benchmarks,
> multi-dataset validation, and final API stabilization are complete.**

What ships, what is experimental, and what is deferred for v0.1. See
[`FEATURE_STATUS.md`](FEATURE_STATUS.md) for per-feature status and
[`PRODUCT_SCOPE_AUDIT.md`](PRODUCT_SCOPE_AUDIT.md) /
[`V0_1_PRUNING_PLAN.md`](V0_1_PRUNING_PLAN.md) for the module-level audit.

## Core (v0.1, default-safe, validated)
- Reference-based **bulk** RNA-seq deconvolution (wNNLS / protocol-aware pipeline).
- Reference-based **spatial** (10x Visium) deconvolution (NB-CAR).
- **Broad → fine hierarchy** with **unresolved mass** kept at the family level.
- **Solver `auto`** (gene-masking cross-validation).
- **Reference suitability** score (PASS/CAUTION/WARNING/FAIL).
- **Separability / spillover** diagnostics.
- **QC-first simplified report** (`report.html`) + **technical appendix**
  (`technical_appendix.html`).
- **Honest basic benchmark** (internal baselines; bulk and spatial kept separate).

Core outputs are **RNA-derived composition estimates**, never absolute cell
counts unless explicit mRNA-content correction was requested.

## Experimental (behind explicit flags; not part of the core claim)
- **State-aware three-level hierarchy** (broad → cell type → state) —
  `--state-aware` / `deconv_bulk(state_aware=True)`.
- **Granular signatures** (broad / cell-type-within-family / state-within-cell-type
  gene panels).
- **Spatial multi-metric ranking** (`benchmarks/shared/spatial_multimetric_ranking.py`).
- **External-tool benchmark runners** (MuSiC, CARD, cell2location, …) — skip
  gracefully when not installed.
- **Composite scorecard** (a weighted multi-criteria summary, not objective
  accuracy).
- **Synthetic state-aware benchmark** (`benchmarks/synthetic/`).

Experimental features are labelled experimental in the report and docs, never
change defaults, and are not used to support accuracy claims.

## Deferred / not implemented (must not be presented as available)
- **Reference adaptation** (residual-driven reference/gene-weight updating).
- **Cell-type-specific expression reconstruction.**
- **Full BayesPrism-like Bayesian model** (Gibbs sampling).

## Reporting in v0.1
The **main report** (`report.html`) is concise and QC-first: executive summary →
reference/signature QC → input compatibility → resolution/uncertainty → final
bulk predictions → final spatial predictions → bulk benchmark → spatial benchmark
→ warnings → methods/source-data links. All **full heatmaps, spillover networks,
spot pies, the large signature heatmap, exploratory figures and the full
source-data listing** live in the **technical appendix** (`technical_appendix.html`),
which the main report links to (and which links back). Every figure also writes
a `.data.tsv`; the figure manifest maps each figure to its data.
