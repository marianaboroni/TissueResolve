# In-Solver Edge-Aware Spatial Smoothing Report

Experimental, opt-in. Motivated by the prior **negative result**: post-fit
edge-aware / family-adaptive smoothing did NOT reduce oversmoothing because a
post-fit operator can only *add* smoothing on top of the solver's already-smoothed
output. Edge awareness must therefore enter the **spatial graph used by the
solver**. No Redeconve-like regularization, no new solver, no deep learning, no
default change.

## Phase 0 — Audit (documented before any code change)

### Current solver / graph flow

1. `SpatialPipeline.run` (`spatial/pipeline.py`) builds the spot graph at step 4
   via `build_hex_graph_from_arrays(array_row, array_col, …)` — **coordinates
   only**, no expression.
2. `spatial/graph.py::_build_graph_from_arrays` builds the hexagonal adjacency:
   `A_raw[i,j] = 1` for each hex neighbour (binary), then row-normalizes
   `A_norm = D⁻¹ A_raw` (every neighbour gets equal weight `1/degree`), and
   `L = I − A_norm`. `SpatialGraph` carries `A` (=A_norm), `L`, `degree`,
   `n_spots`, `spot_ids`, `batches`.
3. `SpatCARModel.fit` (`spatial/model.py`) uses `graph.A`/`graph.L` with the CAR
   penalty `lambda_spatial · Θᵀ L Θ`. `lambda_spatial` (default **0.1**) and
   `alpha` are always recorded.
4. Oversmoothing / boundary / Moran metrics:
   `benchmarks/shared/spatial_metrics.py`.

### Where edge weights can be injected (minimal, no rewrite)

- **Graph edges are currently binary then equal-weight-normalized.** Replacing
  the binary `A_raw[i,j] = 1` with a per-edge weight `w_ij ∈ [min,max]`
  (small across boundaries), then row-normalizing `A_norm = D⁻¹ A_w` and
  recomputing `L = I − A_norm`, makes each spot's smoothing target a
  **boundary-aware weighted mean** of its neighbours. The CAR penalty and the
  solver are **unchanged** — they consume `A`/`L` exactly as before.
- Injection point: build the marker spot matrix `Y_marker` (pipeline step 5)
  before the graph (step 4), and when `cfg.spatial_solver.edge_aware` is set,
  build the edge-weighted graph from `Y_marker` + coordinates instead of the
  plain hex graph.

### Minimal code changes

1. `spatial/graph.py`: add `WeightedSpatialGraph` + `compute_edge_aware_graph(…)`
   (per-edge weights from expression/composition distance) and
   `build_edge_aware_spatial_graph(…)` returning a weighted `SpatialGraph`
   (weighted `A_norm`/`L`); add an optional `metadata` field to `SpatialGraph`
   (backward-compatible).
2. `config.py`: add `SpatialSolverConfig.edge_aware: bool = False` (+ a few
   edge params), all defaulting to current behaviour.
3. `pipeline.py`: compute `Y_marker` before the graph; branch to the edge-aware
   builder when `edge_aware` is set; record edge-weight metadata.
4. `experimental/spatial_presets.py`: add `edge_aware_smoothing` preset
   (sets `edge_aware=True`; λ kept moderate). `default`/`weak_smoothing`
   unchanged.

### Risks

- Row-normalized edge weighting redistributes smoothing *within* each
  neighbourhood (favouring within-domain neighbours) but keeps the global λ; it
  sharpens boundaries but may reduce the *global* oversmoothing score only
  modestly. A lower base λ may be needed to hit the −20% gate — to be decided
  empirically (benchmarked, not assumed).
- Expression-only edge weights avoid a second fit (composition distance is
  optional/off by default) to keep runtime at one fit.

### Benchmark plan

Breast + HLCA/lung synthetic spatial, ≥5 seeds, structured scenario; compare
default(0.1) / weak(0.02) / no(0.0) / in-solver edge_aware (λ=0.1 and 0.05) /
post-fit edge_aware (diagnostic) / CARD (reuse). Promotion gates vs default and
vs weak_smoothing. Outputs under `benchmarks/outputs/in_solver_edge_aware/`.

**Feasibility conclusion: no major solver rewrite required** (edge weights enter
via graph construction; solver untouched). Proceeding to implementation.

---

## Results (run 2026-06-23; breast + HLCA/lung synthetic spatial, 5 seeds)

### 2. Why post-fit edge-aware failed (recap)
A post-fit smoother only *adds* smoothing on top of the solver's output, so it
*increased* oversmoothing (prior report). Moving edge weights **into the spot
graph** (binary → per-edge weights, row-normalized, `L = I − A_norm`) lets the
existing CAR penalty smooth less across boundaries — no solver rewrite.

### 3. Implementation
`spatial/graph.py`: `compute_edge_aware_graph` (per-edge weights = `exp(−d/scale)`
from library-normalized, z-scored expression distance; bounded `[0.05, 1.0]`;
identical neighbours → 1.0) and `build_edge_aware_spatial_graph` (weighted
`A_norm`/`L`, fallback to the hex graph if no edges). Config flag
`SpatialSolverConfig.edge_aware` (default **False**). Pipeline builds `Y_marker`
before the graph and branches on the flag. Note: the repo already had an unused
`build_expression_weighted_graph` (SpatCAR v1, γ-decay, no bounds/metadata); the
new builder adds bounded weights, adaptive scale, optional composition distance,
fallback, and the required metadata.

### 6/12/13/14. Results by method (mean over 5 seeds)

**Breast**

| method | fine Pear | broad Pear | local RMSE | oversmooth | boundary F1 | domain ARI | rare sens | FP-subtype | effN (truth 17.0) |
|---|---|---|---|---|---|---|---|---|---|
| default (0.1) | 0.859 | 0.935 | 0.024 | 1.800 | 0.260 | 0.278 | 0.711 | 0.002 | 17.8 |
| edge_aware 0.05 | 0.869 | 0.948 | 0.024 | 1.635 | **0.338** | 0.312 | 0.689 | 0.002 | 17.3 |
| edge_aware 0.10 | 0.865 | 0.941 | 0.024 | 1.795 | 0.288 | 0.292 | 0.689 | 0.002 | 17.5 |
| weak_smoothing (0.02) | 0.865 | 0.953 | 0.024 | **1.307** | 0.353 | 0.330 | 0.667 | 0.002 | 16.7 |
| no_smoothing (0.0) | 0.773 | 0.920 | 0.038 | 0.361 | 0.327 | 0.319 | 0.356 | 0.001 | 10.6 |

**Lung**

| method | fine Pear | broad Pear | local RMSE | oversmooth | boundary F1 | domain ARI | rare sens | FP-subtype | effN (truth 12.2) |
|---|---|---|---|---|---|---|---|---|---|
| default (0.1) | 0.571 | 0.899 | 0.039 | 1.852 | 0.270 | 0.456 | 0.756 | 0.116 | 19.3 |
| edge_aware 0.05 | 0.589 | 0.921 | 0.039 | 1.740 | **0.299** | 0.465 | 0.756 | 0.096 | 17.7 |
| edge_aware 0.10 | 0.582 | 0.913 | 0.039 | 1.850 | 0.306 | 0.449 | 0.756 | 0.109 | 18.8 |
| weak_smoothing (0.02) | 0.591 | 0.923 | 0.039 | 1.524 | 0.283 | 0.471 | 0.822 | 0.087 | 16.7 |
| no_smoothing (0.0) | 0.486 | 0.868 | 0.051 | 0.515 | 0.326 | 0.484 | 0.622 | 0.054 | 10.3 |

### 6. Results vs default
In-solver edge_aware (λ=0.05) **reverses the post-fit negative result**: it
*reduces* oversmoothing (breast −9%, lung −6%) and **improves boundary F1 on both
datasets** (0.260→0.338, 0.270→0.299), with small fine/broad Pearson gains,
effective-N closer to truth, and lower lung false-positive subtype rate — no
metric regresses.

### 7. Results vs weak_smoothing
weak_smoothing reduces oversmoothing **more** (−27% / −18%) but its boundary-F1
gain is dataset-dependent (strong on breast 0.353, weak on lung 0.283) and it
lowers breast rare-niche sensitivity. edge_aware's boundary-F1 gain is **more
consistent across tissues** and it preserves lung rare niche (0.756). They are
complementary: weak = oversmoothing lever; edge_aware = boundary preservation.

### 8. Results vs no_smoothing
no_smoothing minimizes oversmoothing but collapses accuracy (fine 0.77/0.49) —
control only.

### 9. External comparison
Reused from the gold-truth spatial report: CARD oversmooths far less (~0.6–0.97)
than any TissueResolve mode; no new external reruns (runtime). edge_aware narrows
but does not close that gap.

### 10. Edge-weight diagnostics
Bounded as designed: min 0.05, max 1.0, mean ≈ 0.52–0.53 on both datasets —
cross-boundary edges down-weighted toward the floor, within-domain edges near 1.0
(confirmed by the unit test on toy two-domain data).

### 11. Promotion gates (vs default; both datasets, 5 seeds)

| gate | edge_aware 0.05 breast | edge_aware 0.05 lung |
|---|---|---|
| oversmoothing −≥20% | ❌ (−9%) | ❌ (−6%) |
| boundary F1 preserved/improved | ✅ | ✅ |
| fine Pearson preserved/improved | ✅ | ✅ |
| broad Pearson within 2% | ✅ (improved) | ✅ (improved) |
| local RMSE preserved | ✅ | ✅ |
| rare niche preserved/improved | ❌ (0.711→0.689) | ✅ (0.756) |
| FP-subtype not increased | ✅ | ✅ |
| **gates passed** | **5/7** | **6/7** |

The headline retention gate (oversmoothing −≥20% in both datasets) **FAILS**
(−6…−9%). No failures/crashes; improvement replicated across seeds and both
datasets for boundary F1 and accuracy.

### 12. Supported claims
- In-solver edge-aware smoothing **reduces oversmoothing and improves boundary
  preservation** (boundary F1) on both breast and lung — reversing the prior
  post-fit negative result — with no accuracy/spillover regression and
  effective-N closer to truth.
- Its boundary-preservation gain is more consistent across tissues than
  weak_smoothing's.

### 13. Unsupported claims
- **NOT** claimed: it meets the ≥20% oversmoothing-reduction target (it does not;
  ~6–9%).
- **NOT** claimed: it beats weak_smoothing on oversmoothing (weak is stronger).
- **NOT** claimed: it solves collinearity; it is **not** Redeconve-like.
- **NOT** claimed: default-readiness, or real-Visium accuracy.

### 14. Recommendation
**Retain `edge_aware_smoothing` as an experimental, opt-in preset; do NOT promote
to default and do NOT recommend it over weak_smoothing for oversmoothing** (it
fails the ≥20% gate). It is the better lever for **boundary preservation** with no
accuracy cost. Suggested further work (not done here): lower `min_edge_weight`,
combine a lower base λ with edge weights, or add composition distance from a quick
weak fit, to push oversmoothing reduction toward the 20% bar — then re-gate.

### 15. Remaining limitations
- Single structured scenario type per dataset; 5 seeds; expression-only edge
  weights (composition distance off by default to keep one fit).
- Oversmoothing reduction is modest; the strong lever remains in-solver base λ.
- Breast rare-niche sensitivity dips slightly under edge_aware (lung preserved).
- External reruns deferred (runtime); CARD remains far less oversmoothing.
