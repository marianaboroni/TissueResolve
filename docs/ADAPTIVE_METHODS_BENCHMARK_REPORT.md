# Adaptive Methods Benchmark Report

Conservative, interpretable, **experimental / opt-in** improvements targeting the
weaknesses from prior benchmarks (spatial oversmoothing, tissue-dependent
rare-niche behaviour, inflated effective-N, weak within-family conditional
recovery, unresolved-mass calibration, bulk raw accuracy, diagnostic auto mode).

> No Redeconve-like state-graph regularization. No deep model. No solver rewrite.
> No stable defaults changed. Everything below is opt-in and labelled
> experimental unless explicitly promoted (nothing is promoted to default in
> this task).

## Phase 0 — Implementation map, hooks, plan

### Current implementation map (audited)

- **Spatial smoothing:** `SpatialSolverConfig.lambda_spatial` (default **0.1**;
  `src/tissueresolve/config.py`). The NB-CAR solver (`spatial/model.py`) consumes
  a `SpatialGraph` (`spatial/graph.py`) with row-normalized weighted adjacency
  `A` and Laplacian `L = I − A`; the penalty is `lambda_spatial · Θᵀ L Θ`.
  `SpatialPipeline.run` (`spatial/pipeline.py`) builds the graph internally via
  `build_hex_graph_from_arrays` — **no public graph-injection hook**.
- **weak_smoothing preset:** `experimental/spatial_presets.py`
  (`apply_spatial_preset`), CLI `--spatial-preset`. Sets λ only; default is a
  no-op.
- **Soft gating / unresolved mass:** `bulk/hierarchical.py`,
  `experimental/soft_hierarchy/partial_gating.py`, `resolution.py`. Predictions
  carry `unresolved_<family>` columns; `trusted_resolution` per family
  (broad_only / selected_fine / full_fine) recorded in metadata.
- **Resolution Decision Layer:** `resolution.py` (recent addition) — decides
  per-family resolution from separability/support diagnostics.
- **Gene weighting / protocol:** `compute_modality_aware_gene_weights`
  (results.py), donor-CV gene filters (`reference/gene_filters.py`), WNNLS
  weighting in `bulk/solver.py`.
- **Auto mode:** `api.deconv_bulk/deconv_spatial` resolve `resolution_mode`
  (`auto` → currently maps to a fixed hierarchical path), not yet
  diagnostic-driven.
- **Benchmarks:** `benchmarks/diagnostics/gold_truth_performance_benchmark.py`
  (bulk + spatial-like, per-sample predictions, paired stats via
  `gold_truth_paired_stats.py`); `benchmarks/spatial/run_synthetic_spatial.py`
  and `run_weak_smoothing_grid.py` (breast + lung, multi-seed, generalized
  `domain_families`); external runners (`run_card_synthetic.R`, RCTD/c2l).

### Existing hooks to reuse (no rewrite)

- Edge-aware + family-adaptive smoothing: **post-fit experimental wrapper** over
  the proportions from an initial fit + the spot graph (reuses
  `benchmarks.shared.spatial_metrics` neighbour utilities and `SpatialGraph`
  semantics). The spec explicitly allows "an experimental wrapper that reweights
  smoothing after an initial fit" — chosen to avoid touching the solver.
- Unresolved calibration: operates on the existing `unresolved_<family>` +
  per-family `trusted_resolution` outputs.
- Auto mode: consumes existing reference/query QC + resolution decisions +
  (optional) spatial diagnostics.
- Donor-stable weighting: reuses donor-CV machinery conceptually; pure function
  over the reference AnnData.

### Files to create (minimal)

- `src/tissueresolve/experimental/spatial_adaptive_smoothing.py` (Phases 1+2)
- `src/tissueresolve/experimental/unresolved_calibration.py` (Phase 3)
- `src/tissueresolve/experimental/auto_mode.py` (Phase 4)
- `src/tissueresolve/experimental/bulk_donor_stable_weighting.py` (Phase 5)
- `benchmarks/diagnostics/adaptive_method_benchmark.py` (screening harness)
- tests: `test_spatial_adaptive_smoothing.py`, `test_unresolved_calibration.py`,
  `test_bulk_donor_stable_weighting.py`, `test_auto_mode_diagnostics.py`

### Files to modify

- None of the core algorithm/default paths. (Optional, deferred: a
  `--spatial-preset family_adaptive_smoothing` CLI hook once the post-fit
  wrapper is validated — not wired in this task to keep defaults untouched.)

### Risks

- Post-fit smoothing is a different operator than the in-solver CAR penalty;
  results are indicative of the edge-aware *principle*, not identical to an
  in-solver implementation. Documented as such.
- Full benchmark matrix (2 datasets × 8 scenarios × 5 seeds × 7 spatial configs
  + external reruns) is **runtime-infeasible** in one session (hundreds of
  fits). Per the spec's staged approach, components are **screened** on a subset
  and the full matrix is deferred/labelled.

### Benchmarks to run (staged)

1. Unit tests for every module (offline, deterministic) — primary correctness.
2. Screening benchmark: edge-aware + family-adaptive post-fit smoothing vs
   default/weak on the breast synthetic scenario (validation seeds); bulk
   donor-stable weighting screen on breast pseudobulk. Promotion gates computed.
3. Deferred (documented, not silently skipped): full 2-dataset × 8-scenario ×
   5-seed × 7-config spatial matrix, external-tool reruns, lung replication.

---

## Results (screen run 2026-06-23; breast synthetic spatial, 2 seeds)

### 1. Motivation
Address prior weaknesses (oversmoothing, tissue-dependent rare niches, inflated
effective-N, weak conditional recovery, unresolved calibration, bulk accuracy,
diagnostic auto mode) with conservative, opt-in, interpretable components — no
Redeconve-like regularization, no deep model, no solver rewrite, no default change.

### 2. Components tested
A. Edge-aware spatial smoothing (post-fit). B. Family-specific spatial lambda
(post-fit). C. Family-specific unresolved calibration. D. Diagnostic-driven auto
mode. E. Donor-stable bulk gene weighting. (A,B benchmarked against gold truth;
C,D,E validated by unit tests — end-to-end pipeline wiring deferred, recorded in
`method_status.tsv`.)

### 3. Datasets and scenarios
Breast synthetic spatial (sharp border + gradient + rare niche + mixed/collinear
spots), 2 validation seeds. Full 2-dataset x 8-scenario x 5-seed x 7-config
matrix + external (CARD/RCTD/cell2location) reruns DEFERRED for runtime
(`method_status.tsv`), not silently skipped.

### 5/6/7. Spatial results, edge-aware, family-lambda (mean over 2 seeds)

| method | fine Pear | broad Pear | local RMSE | oversmooth | boundary F1 | rare sens | FP-subtype |
|---|---|---|---|---|---|---|---|
| default (λ=0.1) | 0.853 | 0.938 | 0.024 | 1.838 | 0.298 | 0.556 | 0.002 |
| weak_smoothing (λ=0.02) | 0.857 | 0.954 | 0.024 | 1.357 | 0.357 | 0.611 | 0.002 |
| no_smoothing (λ=0.0) | 0.757 | 0.916 | 0.038 | 0.373 | 0.296 | 0.444 | 0.001 |
| edge_aware (post-fit) | 0.853 | 0.938 | 0.024 | 1.935 | 0.303 | 0.667 | 0.002 |
| family_adaptive (post-fit) | 0.853 | 0.938 | 0.024 | 1.911 | 0.301 | 0.556 | 0.002 |

**Key negative result:** the **post-fit** edge-aware and family-adaptive smoothers
do **not** reduce oversmoothing — they slightly *increase* it (1.84 → 1.91–1.94),
because a post-fit operator can only *add* smoothing on top of the solver's
already-smoothed output; it cannot undo the in-solver CAR penalty. Oversmoothing
is best reduced **in-solver** (weak_smoothing λ, −26% here). Edge-aware did show a
side benefit — higher rare-niche sensitivity (0.556 → 0.667) with no accuracy or
spillover cost — but that is not its stated target and is single-dataset/2-seed.

### 8. Unresolved calibration
Unit-tested only: mass-conserving overall and per broad family; broad-only stays
broad-only; supported subtypes retain more mass than unsupported; rare
marker-supported subtypes are not zeroed; disableable (`mode='none'`). End-to-end
benchmark impact deferred (pipeline wiring).

### 9. Diagnostic auto mode
Unit-tested only: conservative under incomplete evidence; high boundary →
edge-aware (experimental), strong rare-niche / effN-inflation → weak_smoothing;
poor fine reliability → broad_only; **never overrides a user-specified
mode/preset** (warns only); reasons recorded. Advisory; not wired to override
defaults.

### 10. Donor-stable bulk weighting
Unit-tested only: finite/non-negative weights, no gene removal; donor-unstable
and non-specific genes downweighted relative to stable type-specific markers;
safe fallback when the donor column is absent. End-to-end bulk accuracy impact
deferred (requires solver wiring; default weighting unchanged).

### 11. External spatial comparison
Reused from prior reports (breast/lung gold-truth spatial): CARD oversmooths far
less than TissueResolve; weak_smoothing partially closes that gap in-solver. No
new external reruns in this screen (runtime).

### 12. Promotion gates (vs default; spatial screen, breast, 2 seeds)

| method | gates passed | oversmoothing gate | verdict |
|---|---|---|---|
| weak_smoothing | 7/7 | PASS (−26%) | experimental (already a preset; lung not fully passing — see weak-smoothing report) |
| edge_aware (post-fit) | 6/7 | **FAIL** (adds smoothing) | rejected for oversmoothing; rare-niche side-benefit → further validation (in-solver) |
| family_adaptive (post-fit) | 6/7 | **FAIL** | rejected for oversmoothing |
| no_smoothing | 3/7 | PASS but accuracy collapses | control only |

No component meets the **default-promotion** bar (needs ≥2 datasets, multiple
scenario classes, paired stats). **Nothing promoted to default.**

### 13. Supported claims
- weak_smoothing (in-solver λ=0.02) reduces oversmoothing (~26% on breast here)
  and preserves/improves broad/fine accuracy — consistent with the dedicated
  multi-dataset weak-smoothing report; it remains an experimental opt-in preset.
- The experimental modules (unresolved calibration, auto mode, donor-stable
  weighting) behave as specified in unit tests (mass conservation, advisory
  non-override, sensible gene ranking).

### 14. Unsupported claims
- **NOT** claimed: post-fit edge-aware/family-adaptive smoothing reduces
  oversmoothing — it does not (negative result).
- **NOT** claimed: any component improves end-to-end bulk accuracy — not yet
  benchmarked end-to-end.
- **NOT** claimed: any component is ready for default; none is promoted.
- **NOT** claimed: real-Visium accuracy; cross-tissue generalization of these
  screens (breast only, 2 seeds).

### 15. Recommendations
- Keep weak_smoothing as the experimental lever for oversmoothing (in-solver).
- Re-implement the edge-aware **principle in-solver** (inject an edge-reweighted
  `SpatialGraph` into `SpatialPipeline` — a pipeline graph-injection hook, not a
  solver rewrite) before judging it; the post-fit form is rejected.
- Carry unresolved calibration, auto mode, and donor-stable weighting to
  end-to-end benchmarks (wire opt-in) as the next staged step.

### 16. Remaining limitations
- Spatial screen is breast-only, 2 seeds, one scenario type; lung + the full
  scenario/seed/config matrix and external reruns are deferred (runtime).
- C/D/E are validated behaviorally (unit tests), not end-to-end.
- Post-fit smoothing is a different operator than the in-solver CAR penalty.
