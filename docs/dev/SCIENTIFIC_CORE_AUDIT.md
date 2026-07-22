# TissueResolve — scientific core audit

Component-by-component status of the deconvolution core, from reading the actual
source under `src/tissueresolve/`, `benchmarks/`, and `examples/` (not the README).
Status legend: **implemented** (validated, default-safe) · **partial** · **experimental**
(opt-in, not validated across datasets) · **benchmark-only** · **deferred** · **absent** ·
**insufficient** (implemented but scientifically weak for the stated goal).

## Core components

| Component | Where | Status | Notes / scientific assessment |
|---|---|---|---|
| Reference construction (donor-aware aggregation, CPM/L1, NB φ) | `reference/build.py` | implemented | Donor-mean averaging + cross-donor CV. Stores mean profiles only (no per-cell). Now also optional build-time gene selection (this session). |
| Broad signature | `reference/markers.py` (`GeneSelector`) + `reference/gene_selection.py` (`donor_de`) | partial→improving | Legacy `GeneSelector` = composite log-additive score on mean profiles (specificity×stability×protocol×concordance). **Insufficient for fine/rare** (benchmark: worst conditional RMSE). New `donor_de` (this session) beats it donor-disjoint on breast+lung+cross-platform. |
| Within-family fine signature | `reference/within_family_markers.py`, `pairwise_markers.py`; `gene_selection.donor_aware_de(mode="sibling")` | partial | Sibling DE exists; not packaged as a first-class per-family fine panel with donor-disjoint validation. |
| Rare-population confirmation signature | — | **absent** | No specificity-first small confirmation panel exists. Stage 1 §3.2C gap. |
| Donor-aware marker selection | `gene_selection.py` (donor pseudobulk + one-vs-rest/sibling DE) | implemented (experimental) | Real per-donor fold-change support metric; single-donor genes penalised, not dropped. NOT pydeseq2 (py3.9). |
| Batch/study-aware marker stability | `benchmarks/shared/batch_effects.py`; `gene_selection` (partial) | partial | Batch-effect diagnostics exist; leave-one-**study**-out selection not implemented (Stage 1 §3.3 gap). |
| Protocol-aware gene weighting | `reference/gene_filters.py`, `markers.py`, `protocol/` | implemented | Curated bias lists + composite weighting. Solid; unchanged. |
| Bulk wNNLS | `bulk/solver.py` (`WNNLSSolver`), `bulk/pipeline.py` | implemented (default) | Default bulk backbone. Benchmark: **beaten by plain NNLS and by Poisson** on donor-disjoint pseudobulk. |
| Bulk Poisson GLM | `experimental/nb_bulk_solver.py`; `solver/poisson_glm.py` (this session) | experimental (strong) | Best TR bulk backbone (pseudobulk + real-bulk coherence). Wired as `solver="poisson"`, flat + hierarchical. Not default (real-bulk gate incomplete). |
| NB GLM | `experimental/nb_bulk_solver.py`, `nb_dispersion.py`; `solver` `"nb"` | experimental (negative) | NB ≈ or < Poisson consistently (donor-level dispersion near-inert). Documented negative. |
| Spatial NB-CAR | `spatial/model.py`, `graph.py`, `pipeline.py` | implemented (default) | NB likelihood + CAR prior, block descent, global λ. **Not yet audited for boundary/rare-niche behaviour** (Stage 7). |
| Hierarchical broad→fine | `bulk/hierarchical.py`, `spatial/hierarchical.py`, `reference/hierarchy.py` | implemented (default) | Family-level + fine + conditional split + gating. Poisson now composes (this session). |
| Resolution decision layer | `resolution.py`, `reference/resolution.py` | implemented | Per-family `broad_only`/`selected_fine`/`full_fine` from mean-profile reliability. **Global, not sample/spot-specific** (Stage 4 gap). |
| Soft gating | `experimental/soft_hierarchy/partial_gating.py` + `hierarchy.assemble_hierarchical_estimates` | implemented (default) | Confidence from **mean-profile separability** → pessimistic for collinear-learnable families. Audit (this session): abstention mostly correct (13/15); confidence not independently calibrated (Stage 5 gap). |
| Separability | `reference/separability.py` | implemented | Bhattacharyya/Jeffreys/Pearson; drives gating. |
| Spillover | `benchmark/spillover.py` | implemented | Mixture-based spillover matrix. |
| Bootstrap uncertainty | `uncertainty/bootstrap.py` | partial/insufficient | Resampling CIs; **pools cells / single source** — not decomposed into donor/reference/signature (Stage 8 gap). |
| Rare-cell evaluation | `experimental/rare_detection.py` | experimental/insufficient | Presence gate + PR curve exist, but **no calibrated detection probability, no LOD/LOQ, no absent-vs-not-testable status** as first-class per-sample output (Stage 2 gap). |
| Reference suitability | `reference/suitability.py` | implemented | 8-component PASS/CAUTION/WARNING/FAIL. |
| Query/reference compatibility | `protocol/mismatch.py`, `suitability.py` | implemented | Diagnostic only; **no adaptation** (Stage 6 gap). |
| Unknown / out-of-reference | — | **absent** | Proportions renormalise over known types; a true absent population is force-assigned. Only Rectangle-style scalar exists externally. Stage 3 gap. |
| Spatial smoothing | `spatial/auto_params.py`, `experimental/spatial_adaptive_smoothing.py`, `spatial_presets.py` | partial/experimental | Global λ auto-select (default 0.1); edge-aware/family-specific re-smoothing is post-hoc experimental, not promoted. Stage 7. |
| AutoSolver | `solver/auto.py`; `solver/auto_composition.py` (this session) | implemented / experimental | `auto` (gene-masking CV) **misaligned** with composition (benchmark). New `auto_composition` (composition-calibrated, +Poisson) fixes regret 0.012→0 on breast. |
| Benchmark generation | `benchmarks/`, `rectangle_comparison.py` (this session), `spatial/benchmark.py` | implemented | Donor-disjoint pseudobulk; config-driven Rectangle harness; cross-platform via 10x chemistries. Study-disjoint / multi-tissue beyond breast+lung not yet (Stage 9). |
| External-tool integration | `benchmarks/**`, `rectangle_env_runner.py` | partial | MuSiC/Bisque/CARD/RCTD/cell2location + Rectangle (isolated py3.11) executable; others export/skip honestly. |

## Duplication / dead paths / mismatches found

- **Two gene-selection systems** now coexist: legacy `GeneSelector` (mean-profile composite;
  the pipeline default) and new `gene_selection.donor_de` (donor-aware; opt-in via
  `build_reference(gene_selection="donor_de")`). Intended (opt-in, benchmark-gated) but must
  be consolidated once donor_de clears the real-bulk gate.
- **Two Poisson entry points**: `experimental/nb_bulk_solver.NBGLMBulkSolver` (pipeline path,
  `cfg.bulk_solver.method`) and `solver/poisson_glm.PoissonGLMSolver` (backbone, `solver=`).
  Same `fit_bulk_nb_glm` core; different gene-panel behaviour (pipeline honors `gene_panel`;
  backbone honors `ref.selected_genes` after this session's fix). Documented, not a bug.
- **Resolvability computed in two places**: `reference/resolution.py` and
  `reference/hierarchy.evaluate_within_family_resolvability`. Both mean-profile-based;
  potential drift — consolidate later.
- No doc/code contradiction blocking Stage 1 was found; FEATURE_STATUS.md matches code.

## No bug found requiring a regression fix in this audit.
Every change this session has been additive/opt-in with tests; the full suite is green (1457).

## Gap analysis vs the Stage programme (§3–§11)

- **Stage 1 (signature optimisation):** donor-aware pseudobulk, broad + within-family
  candidate scoring, deterministic selection, donor-disjoint benchmark, and comparison to the
  current path are **already implemented and validated this session** (donor_de beats current
  markers on breast+lung+cross-platform; hybrid/ml = negatives). **Missing:** first-class
  packaging (`ReferenceSignatureOptimizer`) producing the **three signature levels** (broad /
  within-family fine / **rare-confirmation, absent**), **minimal-within-tolerance** selection,
  the **manifest TSVs**, and the **fallback status** vocabulary (PASS / PASS_WITH_RESTRICTIONS
  / BROAD_ONLY / EXPERIMENTAL_FINE / REFERENCE_INADEQUATE). → This session's vertical slice.
- **Stage 2 (rare-cell detection):** calibrated probability, LOD/LOQ, status vocabulary — gap.
- **Stage 3 (unknown component):** absent — gap.
- **Stage 4 (sample/spot resolution):** global only — gap.
- **Stage 5 (gating calibration):** audited this session (mostly-correct abstention;
  learnability-gate falsified) — confidence still not independently calibrated.
- **Stage 6 (reference adaptation):** absent/deferred — gap.
- **Stage 7 (spatial smoothing):** partial/experimental — gap.
- **Stage 8 (uncertainty decomposition):** single-source bootstrap — gap.
- **Stage 9 (leakage-safe benchmark):** donor + cross-platform done; study/multi-tissue gap.
