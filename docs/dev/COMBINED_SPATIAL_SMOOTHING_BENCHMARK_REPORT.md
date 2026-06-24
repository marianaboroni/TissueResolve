# Combined Spatial Smoothing Benchmark Report

Experimental, opt-in. Tests whether **combining** the two existing spatial
smoothing levers helps: `weak_smoothing` (low global λ=0.02) + `edge_aware_smoothing`
(in-solver edge-weighted graph). No Redeconve-like / state-similarity
regularization, no new solver, no default change.

## Phase 0 — Audit (documented before code changes)

### 1. Current `weak_smoothing`
`experimental/spatial_presets.py`: sets `lambda_spatial = 0.02`, `edge_aware = False`.
Prior result: oversmoothing −27% (breast) / −18% (lung); boundary-F1 gain
dataset-dependent; breast rare-niche dips. (Default λ=0.1 unchanged.)

### 2. Current `edge_aware_smoothing`
`spatial/graph.py::build_edge_aware_spatial_graph` builds a row-normalized
**weighted** graph (per-edge `w_ij = exp(−d_ij/scale)`, bounded `[min,max]`,
expression distance) consumed by the unchanged NB-CAR solver; preset sets
`lambda_spatial = 0.05`, `edge_aware = True`. Prior result: boundary-F1 improved
on both datasets, accuracy preserved, but oversmoothing only −9% / −6% (failed
the ≥20% gate).

### 3. Where combined behaviour is configured
`SpatialSolverConfig` already carries `lambda_spatial`, `edge_aware`,
`edge_aware_min_weight`, etc. The pipeline already branches to the weighted graph
when `edge_aware` is set. **The combination is simply `λ=0.02` + `edge_aware=True`**
— no new mechanism. `apply_spatial_preset` is extended to also set
`edge_aware_min_weight` from a preset and to record `combined`/`min_edge_weight`.

### 4. Minimal code changes
- `spatial_presets.py`: add `combined_weak_edge_smoothing` preset (λ=0.02,
  edge_aware=True, min_edge_weight=0.05); add `combined`/`min_edge_weight` to
  `SpatialPresetInfo` + metadata; honour `min_edge_weight` in `apply_spatial_preset`.
- `cli.py`: add the preset to `--spatial-preset` choices.
- New benchmark `benchmarks/diagnostics/combined_spatial_smoothing_benchmark.py`.
- Tests appended to `tests/test_in_solver_edge_aware.py`.
- **No** change to defaults, solver, soft gating, unresolved mass, or bulk.

### 5. Benchmark plan
breast + lung synthetic spatial, 5 seeds; 9 configs: default / no / weak /
edge_aware / combined(λ0.02,mew0.05) + variants (λ0.02 mew0.01/0.03; λ0.03 mew0.03;
λ0.05 mew0.03). Promotion gates vs default, weak, and edge_aware. CARD/RCTD/c2l
reused qualitatively (no new slow reruns).

### 6. Risks
- At λ=0.02 the CAR penalty is already weak, so edge weights have a small
  *absolute* effect; the combination may behave close to weak_smoothing with a
  modest boundary gain (a compromise rather than a strict improvement on both
  targets). To be decided empirically.
- Combination may inherit weak_smoothing's tissue-dependent rare-niche dip.

**Feasibility: no solver rewrite** (combination is a config of existing levers).
Proceeding.

---

## 1. Motivation
Two experimental levers each improved one axis but not the other: `weak_smoothing`
(λ=0.02) cut global oversmoothing but its boundary/rare-niche behaviour was
tissue-dependent; in-solver `edge_aware_smoothing` (λ=0.05, edge-weighted graph)
improved boundary F1 on both tissues but cut oversmoothing only ~6–9% (failed the
≥20% gate). This task asks, empirically: does running **both at once** (low λ *and*
edge weighting) get weak's oversmoothing reduction *and* edge_aware's boundary
preservation, without new machinery.

## 2. Why they were expected to be complementary
λ scales the global CAR penalty `λ·Θᵀ L Θ`; edge weighting reshapes `L` so that
expression-dissimilar (likely cross-boundary) neighbours are coupled less. In
principle λ controls *how much* smoothing and the edge graph controls *where*. The
risk flagged in Phase 0 §6 was that at λ=0.02 the penalty is already so weak that
reshaping `L` has little *absolute* effect — i.e. the combination would collapse
onto `weak_smoothing`. The benchmark tests which it is.

## 3. Implementation
Pure configuration of existing levers: `combined_weak_edge_smoothing` =
`lambda_spatial=0.02`, `edge_aware=True`, `min_edge_weight=0.05`. No new solver, no
change to defaults, soft gating, unresolved mass, or bulk. Preset + CLI choice +
metadata fields (`combined_preset`, `min_edge_weight`) added; `apply_spatial_preset`
honours `min_edge_weight`.

## 4. Datasets and scenarios
Synthetic gold-truth spatial scenarios on **breast** (`real_breast_cancer` ref) and
**HLCA/lung** (`hlca_subset.h5ad`, 49 ref types, rare = NK cells), 20×20 hex grid,
40 cells/spot, **5 seeds each** (0–4). Each scenario contains a sharp boundary,
a gradient, a rare niche, mixed spots, and collinear fine subtypes.

## 5. Configurations tested (9)
default (λ0.1) · no_smoothing (λ0.0) · weak_smoothing (λ0.02) · edge_aware_smoothing
(λ0.05, edge_aware) · **combined_l0.02_mew0.05** (named preset) · variants
combined_l0.02_mew0.01 · combined_l0.02_mew0.03 · combined_l0.03_mew0.03 ·
combined_l0.05_mew0.03. (edge_scale/composition-distance variants not run — not
exposed as stable levers; not invented for this task.)

## 6. Metrics
Full accuracy / spatial-structure / composition / runtime / edge-diagnostic panel
from `benchmarks/shared/spatial_metrics.py` + `metrics.py`, written to the 8 TSVs
under `benchmarks/outputs/combined_spatial_smoothing/`. Values below are 5-seed
means per dataset.

## 7. Results vs default (named preset combined_l0.02_mew0.05)
| dataset | metric | default | combined | Δ |
|---|---|---|---|---|
| breast | oversmoothing | 1.800 | **1.293** | −28.2% |
| breast | boundary F1 | 0.260 | **0.359** | +0.099 |
| breast | fine Pearson | 0.859 | 0.868 | +0.009 |
| breast | rare-niche sens. | 0.711 | 0.667 | **−0.044** |
| breast | conditional RMSE | 0.133 | 0.150 | **+12.8%** |
| lung | oversmoothing | 1.852 | 1.515 | −18.2% |
| lung | boundary F1 | 0.270 | 0.286 | +0.016 |
| lung | fine Pearson | 0.571 | 0.596 | +0.025 |
| lung | rare-niche sens. | 0.756 | 0.778 | +0.022 |
| lung | conditional RMSE | 0.276 | 0.284 | +3.0% |

Combined improves oversmoothing, boundary F1 and fine accuracy vs default on both
tissues, but (a) on breast it dips rare-niche sensitivity and worsens conditional
within-family RMSE >5%; (b) on lung the oversmoothing reduction is only −18%,
below the 20% gate.

## 8. Results vs weak_smoothing
| dataset | metric | weak | combined | verdict |
|---|---|---|---|---|
| breast | oversmoothing | 1.307 | 1.293 | tie (−1.1%, negligible) |
| breast | boundary F1 | 0.353 | 0.359 | tie (+0.006) |
| breast | rare-niche sens. | 0.667 | 0.667 | identical |
| lung | oversmoothing | 1.524 | 1.515 | tie (−0.6%, negligible) |
| lung | boundary F1 | 0.283 | 0.286 | tie (+0.003) |
| lung | rare-niche sens. | 0.822 | 0.778 | **worse (−0.044)** |
**Combined ≈ weak_smoothing on every axis except lung rare-niche, where it is
worse.** The edge-aware component contributes essentially nothing at λ=0.02 — the
Phase 0 §6 risk is confirmed. Combined does **not** beat weak on oversmoothing.

## 9. Results vs edge_aware_smoothing
| dataset | metric | edge_aware | combined | verdict |
|---|---|---|---|---|
| breast | oversmoothing | 1.635 | 1.293 | combined better (lower λ) |
| breast | boundary F1 | 0.338 | 0.359 | combined slightly higher |
| lung | oversmoothing | 1.740 | 1.515 | combined better |
| lung | boundary F1 | 0.299 | 0.286 | **edge_aware better (−0.013)** |
Combined's lower oversmoothing comes entirely from low λ, not edge weighting. On
lung, combined does **not** match edge_aware's boundary F1 — it fails
`boundary_f1_close_to_edge`. Combined does **not** reproduce edge_aware's distinct
boundary behaviour.

## 10. Results vs no_smoothing
no_smoothing collapses oversmoothing (breast 0.361 / lung 0.515 — far *below* truth,
i.e. under-smoothed/noisy) at a large accuracy cost (breast fine 0.773, lung 0.486;
local RMSE up to 0.038/0.051) and degrades rare-niche sensitivity (0.356 / 0.622).
Combined keeps default-level fine accuracy and local RMSE while reducing
oversmoothing — clearly preferable to no_smoothing, but no_smoothing confirms the
oversmoothing axis is dominated by λ, not the edge graph.

## 11. External comparison (CARD/RCTD/cell2location)
Reused qualitatively from prior gold-truth reports (no new reruns — runtime). CARD
still reaches substantially lower oversmoothing than any TissueResolve config here
(combined ~1.3 breast / ~1.5 lung remains well above CARD). TissueResolve (incl.
combined) remains competitive or stronger on fine accuracy / local RMSE. Combined
does **not** close the oversmoothing gap to CARD.

## 12. Edge-weight diagnostics
Edge weights are finite and bounded in [min_edge_weight, 1.0] in every run
(breast/lung mean ≈0.50–0.53, median ≈0.51–0.53; max 1.0 at zero expression
distance). Varying `min_edge_weight` 0.01/0.03/0.05 shifts the mean by ≤0.02 and
changes oversmoothing/boundary by <0.01 — the floor is nearly inert at λ=0.02. No
graph fallbacks on real scenarios.

## 13. Rare-niche behaviour
Breast: default 0.711 → combined 0.667 (**dip**, inherited from weak). Lung:
default 0.756 → weak 0.822 → combined 0.778 (combined below weak). Combined neither
matches weak's lung gain nor avoids weak's breast dip. Rare-niche precision and
rare false-positive rate are not worsened. **Net: rare-niche behaviour is
tissue-dependent and not improved.**

## 14. Effective-N and entropy
Truth effective-N: breast 16.96, lung 12.18. Default over-mixes (17.81 / 19.34);
combined is slightly closer to truth (16.57 / 16.36) — comparable to weak
(16.74 / 16.72). Entropy tracks the same ordering. Lung effective-N stays well
above truth for all smoothed configs; only no_smoothing under-shoots (10.3). Modest,
weak-equivalent improvement.

## 15. Conditional within-family accuracy
Combined worsens conditional within-family RMSE vs default on both tissues
(breast 0.133→0.150, **+12.8%**, fails the ≤5% gate; lung 0.276→0.284, +3.0%, within
gate). Expected: less smoothing leaves noisier within-family (collinear-subtype)
splits. Combined does **not** address collinearity.

## 16. Runtime
Combined 15–19 s/fit, comparable to default (11–16 s) and weak (16–20 s); edge graph
construction adds no meaningful cost. **0 failures across 90 fits.**

## 17. Promotion gates (named preset combined_l0.02_mew0.05)
**Neither dataset passes all gates → not promotable.**
- **breast 9/11** — fails: rare-niche sensitivity preserved (0.667 vs 0.711);
  conditional within-family RMSE ≤5% (+12.8%).
- **lung 9/11** — fails: oversmoothing −20% vs default (only −18.2%); boundary F1
  close to edge_aware (0.286 vs 0.299, gap 0.013 > 0.01).
- Cross-dataset gate "beats weak on oversmoothing AND edge_aware on boundary F1" is
  **not** met: combined ties weak on oversmoothing and loses to edge_aware on lung
  boundary F1. Variants (λ 0.03/0.05) only slide along the weak↔edge_aware line; none
  passes all gates either (best is combined_l0.05_mew0.03 on lung at 10/11, but it is
  edge_aware with negligible difference).

## 18. Supported claims
- Combined reduces oversmoothing and improves boundary F1 **vs default** on both
  tissues, at default-level fine accuracy and local RMSE, with no runtime cost or
  failures.
- Combined is **statistically indistinguishable from `weak_smoothing`**: the
  oversmoothing reduction and boundary gain come from λ=0.02; the edge-aware
  component is nearly inert at this λ.

## 19. Unsupported claims
- ✗ "Combined beats weak_smoothing on oversmoothing." (It ties.)
- ✗ "Combined matches edge_aware on boundary preservation." (Loses on lung.)
- ✗ "Combined improves rare-niche behaviour." (Tissue-dependent; worse than weak on
  lung, dips vs default on breast.)
- ✗ "Combined solves collinearity / improves within-family accuracy." (Worsens it.)
- ✗ Redeconve-like / state-similarity — not implemented, not claimed.
- ✗ Any change to default behaviour. (Unchanged: λ=0.1, edge_aware=False.)

## 20. Recommendation
**Reject for promotion. Retain as experimental, opt-in, NOT recommended.** It is a
**compromise dominated by the λ term** and is effectively redundant with
`weak_smoothing` — the edge-aware lever adds no measurable benefit at λ=0.02, and it
fails promotion gates on **both** datasets. Users wanting less oversmoothing should
prefer `weak_smoothing` (simpler, same result); users wanting boundary preservation
should prefer `edge_aware_smoothing` (λ=0.05). The combined preset is kept only
because it is wired, tested, and harmless; it carries no recommendation. Default is
unchanged.

## 21. Remaining limitations
- Synthetic gold-truth scenarios only; no real spatial ground truth.
- The edge-aware lever was tested only at low λ here; whether edge weighting helps
  at intermediate λ (0.05–0.1) is the `edge_aware_smoothing` question, already
  answered (modest oversmoothing reduction, failed ≥20% gate).
- Oversmoothing gap to CARD remains open and is not addressed by smoothing-strength
  tuning alone.
- Collinear/within-family resolution is not improved by any smoothing config; that
  requires a different mechanism (out of scope; no Redeconve-like regularization in
  this task).
