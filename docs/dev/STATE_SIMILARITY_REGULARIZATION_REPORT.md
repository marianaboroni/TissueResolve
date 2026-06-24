# State-Similarity / Redeconve-inspired Regularization Report

Experimental, opt-in, non-default, benchmark-gated. Tests whether regularizing
**similar cell states jointly** (a state-similarity graph + sparsity shrinkage,
applied as a constrained within-family refinement of spatial estimates) improves
the limitations that pure spatial smoothing could not fix: diffuse mass across
collinear fine states, inflated effective-N, and weak conditional within-family
recovery.

> Inspired by the *concept* of state-aware regularization (as discussed for
> Redeconve and related methods). This is **not** a copy of Redeconve's code,
> API, objective, or claims, and no equivalence to Redeconve is asserted.

## Phase 0 — Audit and design (written before any code change)

### 1. Motivation
Prior smoothing work (`weak_smoothing`, post-fit edge-aware, in-solver edge-aware,
combined) improved or traded oversmoothing/boundary but **never** improved
conditional within-family RMSE — combined smoothing actually worsened it
(breast +12.8%). The documented, safe conclusion from that work:

> Combined weak and edge-aware smoothing reduced oversmoothing relative to the
> default but did not improve over weak smoothing alone. This indicates that, at
> low spatial regularization strength, edge reweighting contributes little
> additional benefit. The combined preset is retained as experimental but is not
> recommended over weak_smoothing.

### 2. Why simple smoothing reached its limit
Spatial smoothing acts on the **spot graph** (couples neighbouring spots). It does
nothing about the *within-spot* ambiguity between biologically similar fine states
(e.g. collinear T-cell subtypes): the solver can spread a family's mass across
many near-collinear states at a single spot with little likelihood penalty. That
manifests as inflated effective-N, diffuse false-positive subtype mass, and poor
conditional (within-family) recovery — none of which a spot-graph penalty targets.
The missing lever is a **state graph** that couples similar *states*, not similar
*spots*.

### 3. Current solver architecture (insertion-point audit)
- `spatial/model.py::SpatCARModel` — alternating block coordinate descent:
  (a) multiplicative NB update per spot (guaranteed non-negative, simplex-renormalised),
  (b) optional spatial proximal mixing `π ← (1−α)π_NB + α·Aπ`, α=min(0.5,λ),
  (c) periodic mismatch update. The spatial penalty is a *proximal blend*, not an
  explicit additive Laplacian term in a global objective.
- `spatial/graph.py` — `SpatialGraph(A, L=I−A, …)` (hex or edge-weighted);
  `build_edge_aware_spatial_graph`. These are **spot** graphs.
- `spatial/pipeline.py::SpatialPipeline.run` — orchestrates: marker selection →
  graph → `SpatCARModel.fit` → `Pi` (spot×state) → QC/Moran's/result. **`Pi` is
  available as a plain matrix after step 6 and before result assembly (step 10).**
- `config.py::SpatialSolverConfig` — `lambda_spatial=0.1`, `edge_aware` flags.
  `HierarchicalConfig` carries soft gating + `allow_unresolved` (unresolved mass).
- `reference/resolution.py`, `reference/separability.py` — Bhattacharyya-based
  separability + merge recommender (reused for state grouping). Families come from
  a fine→broad `mapping` (hierarchy TSV / obs broad column).

### 4. Feasible insertion points
- **Option A (in-solver additive state penalty):** would require adding a
  `λ_state·Lstate` term to the multiplicative NB update derivation — i.e. rewriting
  the update step. The current update is a closed-form multiplicative rule, not a
  gradient step, so a state-Laplacian does not slot in cleanly. **This is a major
  solver rewrite → ruled out by the stop condition.**
- **Option B (post-fit constrained within-family refinement) — CHOSEN.** Operate
  on the fitted `Pi` after `SpatCARModel.fit`: within each broad family, couple
  similar states toward their similarity-weighted mean (shrinks diffuse spread) and
  apply sparsity shrinkage, **conserving each spot's broad-family mass exactly** and
  redistributing only among fine states within the same family. Keeps the solver,
  defaults, soft gating, and unresolved mass untouched. Clearly documented as a
  refinement step, not a re-derived objective.
- **Option C (diagnostic only):** state graph + grouping recommendations without
  changing estimates. Implemented as the `adaptive_resolution_experimental` mode.

**Decision: implement B (default mechanism) + C (diagnostic mode). A is not
attempted (would be a major solver rewrite).**

### 5. Proposed experimental design
- `experimental/state_similarity_regularization.py`:
  - `compute_state_similarity_graph(reference_profiles, cell_types, family_map,
    method, k_states, min_similarity, within_family_only)` → `StateSimilarityGraph`
    (sparse edge list over states; within-family-only by default).
  - `apply_state_regularized_refinement(proportions, state_graph, family_map,
    lambda_state, lambda_sparse, preserve_broad_mass)` — within-family graph
    shrinkage + sparsity; conserves broad mass; preserves unresolved_* columns.
  - `apply_sparsity_aware_refinement(proportions, family_map, state_graph,
    lambda_sparse, rare_protection, preserve_broad_mass)`.
  - `recommend_state_groups(state_similarity_graph, separability_metrics,
    family_map, threshold)`.
- `config.py`: new `StateRegularizationConfig` (all default-safe; `enabled=False`).
- `spatial_presets.py` + CLI: presets `state_regularized_experimental`,
  `sparsity_state_regularized_experimental`, `adaptive_resolution_experimental`.
- Pipeline hook: opt-in post-fit step in `SpatialPipeline.run`, gated by
  `cfg.state_regularization.enabled` (default off), using an optional `family_map`
  passed through `deconv_spatial(**kwargs)`. With no family_map, within-family-only
  yields no edges → exact no-op (documented limitation).

### 6. Risks
- The refinement is post-fit; like post-fit edge-aware smoothing it may be too late
  to undo solver decisions (mitigated: it changes the *within-family allocation*,
  which the spot graph never touched — a different axis).
- Aggressive shrinkage could zero true rare states → explicit rare-state protection
  and conservative default λ.
- Collinear states are collinear *because* markers don't separate them; pulling
  them toward a shared mean may reduce false spread but cannot create information —
  may improve effective-N/spillover while leaving fine Pearson ~flat. If conditional
  RMSE does not improve, this is a negative result and will be reported as such.
- Tissue-dependence (as seen for weak_smoothing) is likely.

### 7. Benchmark plan
breast + HLCA/lung synthetic gold-truth scenarios, ≥5 seeds, same structured
scenario. Compare default / no_smoothing / weak / edge_aware / combined /
state_regularized / sparsity_state_regularized / adaptive_resolution /
state+weak. CARD/RCTD/cell2location reused qualitatively (no new slow reruns).
Outputs under `benchmarks/outputs/state_similarity_regularization/`. Promotion
gates (Phase 7) centred on conditional within-family RMSE, effective-N, spillover,
false-positives, rare-niche preservation, broad-mass conservation.

### 8. Stopping criteria
- Stop (and report) if the only faithful implementation requires rewriting the NB
  multiplicative update or replacing the solver (Option A). **Confirmed: Option B
  needs no solver change, so we proceed.**
- Keep diagnostic-only / negative result if gates fail; never recommend without
  paired, replicated benchmark support; never change the default.

---

## Results (Phases 5–9)

**Headline: a negative result for the primary goal.** Coupling similar *states*
does **not** improve conditional within-family recovery. The diagnostic state graph
and grouping are useful; the refinement is retained experimental but **not
recommended**, and the diagnostic `adaptive_resolution` mode is the only part worth
surfacing. Default behaviour is unchanged.

### 1–6. Motivation / design / module
See Phase 0 above. Implemented Option B (post-fit, within-family, mass-conserving
refinement) + Option C (diagnostic grouping). Option A (in-solver state penalty)
was **not** attempted — it would require rewriting the NB multiplicative update
(stop condition). Module enters the pipeline at **step 6b** of
`SpatialPipeline.run` (after `SpatCARModel.fit`, before QC), gated on
`cfg.state_regularization.enabled` (default False), using a `family_map` passed via
`deconv_spatial(**kwargs)`.

### 7. Datasets and scenarios
breast + HLCA/lung synthetic gold-truth, 5 seeds each (0–4), 20×20 hex, 40
cells/spot; sharp boundary + gradient + rare niche + mixed + collinear fine
subtypes. 110 fits, **0 failures**.

### 8. Configurations tested (11)
default · no_smoothing · weak_smoothing · edge_aware_smoothing · combined_weak_edge ·
**state_regularized** (preset: λ_state 0.01, λ_sparse 0.001) · state_regularized_strong
(0.1/0.05, benchmark variant) · **sparsity_state** (preset: λ_sparse 0.01) ·
sparsity_state_strong (0.1, variant) · **adaptive_resolution** (diagnostic) ·
state_regularized_weak (λ=0.02 + 0.1/0.05). Strong variants characterise the
mechanism's ceiling; they do **not** change the shipped preset defaults.

### 9. Metrics
Full accuracy / spatial / composition / runtime / state-graph panel via the reused
scoring harness, written to the 9 TSVs under
`benchmarks/outputs/state_similarity_regularization/`. 5-seed means below.

### 10. Results vs default — conditional within-family RMSE (the target)
| dataset | config | cond-RMSE Δ vs default |
|---|---|---|
| breast | state_regularized | **+1.7%** (worse) |
| breast | sparsity_state | **+2.4%** (worse) |
| breast | *_strong | +29% / +39% (much worse) |
| lung | state_regularized | **+1.6%** (worse) |
| lung | sparsity_state | **+2.1%** (worse) |
| lung | *_strong | +17.5% / +21% (much worse) |

**No configuration improves conditional within-family RMSE on either dataset.** Mild
settings worsen it slightly; strong settings worsen it a lot. The central
hypothesis is refuted.

### 11–13. Results vs weak / edge_aware / combined
State-reg at λ=0.1 leaves oversmoothing ≈ default (1.77–1.80 vs default 1.80 breast;
1.78–1.80 vs 1.85 lung) — far above weak (1.31/1.52). It does not compete on
oversmoothing/boundary (those are λ-driven, and state-reg keeps the default λ). The
only config with low oversmoothing (`state_regularized_weak`, λ=0.02) gets that from
λ, while its fine Pearson and conditional RMSE are the worst of the mild group.

### 14. Results vs CARD/RCTD/cell2location
Reused qualitatively (no reruns). CARD still has substantially lower oversmoothing;
state-reg does not change TissueResolve's standing (it does not touch the
spot-graph axis where the gap lives). No regression introduced.

### 15. Conditional within-family RMSE
Never improved (section 10). This is the decisive metric for the module's premise.
Confirms the Phase 0 §6 risk: collinear states are collinear *because* markers do
not separate them, so post-hoc redistribution of their mass cannot recover the true
within-family split — it can only move mass around, which here adds error.

### 16. Effective-N and entropy
The one axis where the module does something useful: mild configs move effective-N
**modestly toward truth in both datasets** (breast err 0.85→0.29 with sparsity;
lung 7.16→5.84) while preserving fine Pearson within 1%. Strong configs drive lung
effective-N almost exactly to truth (7.16→0.23) but **overshoot breast** (0.85→3.55)
and damage fine accuracy — so the "right" strength is tissue-dependent and cannot be
fixed without per-tissue tuning. Critically, the effective-N gain does **not**
translate into better conditional RMSE.

### 17. Spillover and false positives
Essentially unchanged for mild configs (lung spillover 0.220→0.218; FP rate flat).
Strong configs reduce spillover marginally (0.220→0.211) but that is confounded with
the accuracy loss. **No meaningful, cost-free spillover/FP improvement.**

### 18. Rare-niche behaviour
Mode- and tissue-dependent. **Sparsity** preserves/improves rare-niche sensitivity
(breast 0.711→0.733; sparsity_strong 0.978; lung 0.756 flat → strong 0.822 — the
rare-protection works). **State-graph concentration** *hurts* it (breast 0.711→0.689;
lung 0.756→0.622; state_strong lung 0.356). So redundancy concentration can suppress
a true rare state when it is similar to a more-supported sibling — a real risk.

### 19. Oversmoothing and boundary
Unchanged at default λ (state-reg acts within-spot, within-family — not on the spot
graph). boundary_f1 identical to default for pure state-reg configs.

### 20. Runtime
State-reg adds no cost (≈9 s/fit). 0 failures / 110 fits. Broad-family mass
conserved exactly (deviation = 0.0) in every config; unresolved mass preserved.

### 21. Promotion gates
**No configuration passes all gates on either dataset.** The primary gate
(conditional within-family RMSE improved ≥5%) **fails for every config on both
datasets**. Best gate counts: `adaptive_resolution` 10/12 (breast), 11/12 (lung) —
but only because it leaves estimates unchanged (it cannot *improve* cond-RMSE).
`state_regularized`/`sparsity_state` reach 9–10/12 but fail the primary gate and
(for state_regularized) rare-niche preservation. Strong variants fail many gates
(fine accuracy, effN-overshoot, rare niche).

### 22. Supported claims
- The module is correct and safe: exact broad-family mass conservation, valid
  compositions, preserved unresolved mass, 0 failures, no runtime cost, default
  untouched.
- Mild sparsity moves **effective-N modestly toward truth in both datasets** while
  preserving fine Pearson within 1%.
- The **state-similarity graph + `adaptive_resolution` grouping** correctly identify
  the collinear within-family clusters (breast: 5 T/NK states, 5 myeloid states, …;
  lung: fibroblast/epithelial subsets) and are useful *interpretation* outputs.

### 23. Unsupported claims
- ✗ "State regularization improves conditional within-family recovery." (Refuted —
  never improves; usually worsens.)
- ✗ "It reduces spillover / false positives." (No meaningful, cost-free gain.)
- ✗ "It reduces oversmoothing or closes the CARD gap." (It does not act on the spot
  graph.)
- ✗ "It improves rare niches." (Sparsity may; concentration hurts — net
  mode/tissue-dependent.)
- ✗ "It solves collinearity." (It cannot create missing marker information.)
- ✗ Any equivalence to Redeconve; any change to default behaviour.

### 24. Recommendation
**Reject for promotion and reject as a recommended preset.** The estimation
refinement (`state_regularized`, `sparsity_state`) is **retained as experimental,
opt-in, NOT recommended** — it fails the primary gate on both datasets and its
only real effect (effective-N toward truth) is tissue-dependent and bought with
conditional-RMSE / fine-accuracy / rare-niche costs. **Retain the diagnostic
`adaptive_resolution` mode + state graph + group recommendations** as safe
interpretability tooling (estimates unchanged). For genuine within-family
ambiguity, prefer the existing soft-gating / `unresolved_<family>` machinery plus
this diagnostic grouping over forced mass concentration. Default unchanged.

### 25. Remaining limitations
- Synthetic gold-truth only; no real spatial within-family ground truth.
- Post-fit refinement cannot recover information the markers do not contain;
  improving conditional recovery likely needs better within-family discriminative
  features or an in-solver state penalty (a solver rewrite, deferred).
- Effective-N calibration is tissue-dependent (no single λ_sparse fits both).
- The diagnostic grouping is reporting-only; it does not (and by rule must not)
  silently collapse fine states in default mode.

