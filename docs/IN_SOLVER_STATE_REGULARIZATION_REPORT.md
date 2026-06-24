# In-Solver State-Similarity Regularization Report (Option A)

Experimental, opt-in, non-default, benchmark-gated. Tests whether putting the
state-similarity + sparsity penalties **inside a solver objective** (rather than the
post-fit refinement of the previous task) improves conditional within-family
recovery and composition calibration.

> Inspired by the *concept* of state-aware regularization. This is NOT Redeconve,
> not Redeconve-equivalent, not the default solver, and not state-of-the-art. No
> external code/API/claims were copied. Module name:
> `state_regularized_solver_experimental`.

## Phase 0 — Solver feasibility audit (written before any code)

### 1. Current solver architecture
`spatial/model.py::SpatCARModel` — alternating **block coordinate descent**:
(a) a closed-form **multiplicative NB update** per spot (guaranteed non-negative,
simplex-renormalised), (b) optional **spatial proximal mixing**
`π ← (1−α)π_NB + α·Aπ`, α = min(0.5, λ_spatial), (c) periodic mismatch update.
`spatial/graph.py::SpatialGraph` carries row-normalised adjacency `A` (csr),
Laplacian `L`, `degree`, `batches`.

### 2. Current objective
There is **no single explicit global objective** that the production solver
minimises by gradient. The NB update is a multiplicative majorise-minimise step on
the NB negative log-likelihood; the spatial step is a *proximal blend*, not a
gradient of an additive `λ·Θᵀ L Θ` term. So the production path is not written as
`min_Θ recon + λ_spatial·spatial + …`.

### 3. How `lambda_spatial` enters
Only through α = min(0.5, λ) in the proximal mixing step (b). It is **not** a
coefficient on an additive penalty in a joint objective.

### 4. Can a state penalty be inserted safely into the production solver?
**No — not without a rewrite.** Adding `λ_state·(state penalty)` to step (a) would
require re-deriving the multiplicative update (the penalty breaks the clean
NB multiplicative form) or bolting a gradient step onto a non-gradient solver.
Either destabilises a validated solver. This is exactly the "major rewrite" the
stop condition forbids.

### 5. Risks of modifying the production solver
Loss of the NB multiplicative non-negativity/monotonicity guarantees; convergence
regressions on the existing breast/lung validation; entanglement of the spatial
proximal step with a new additive term. Unacceptable for a validated default.

### 6. Proposed experimental route (CHOSEN)
Per rule 14 + the decision rule: **a separate experimental solver module** that
consumes the same inputs and returns TissueResolve-compatible `Θ`. It optimises an
explicit joint objective by **projected gradient descent** on a Gaussian
(least-squares) approximation of the reconstruction term — the standard,
well-behaved surrogate for simplex-constrained deconvolution:

```
min_Θ  ||Y_norm − Θ R||²_F
       + λ_spatial · tr(Θᵀ L_sym Θ)                     (spot smoothing)
       + λ_state   · P_state(Θ)                           (state coupling)
       + λ_sparse  · Σ_ik Θ_ik log Θ_ik                   (negative-entropy sparsity)
   s.t. Θ ≥ 0, rows on the simplex (or per-family mass if preserve_broad_mass)
```

with two state-penalty variants:
- **competition** (default): `P = Σ_spots Σ_{a≠b} s_ab Θ_ia Θ_ib = tr(Θ A_state Θᵀ)`
  — discourages co-assigning mass to *several similar states in the same spot*
  (targets diffuse effective-N / false-positive subtype mass);
- **laplacian**: `P = Σ_{a,b} s_ab ||Θ_:a − Θ_:b||² = tr(Θ L_state Θᵀ)` — smooths the
  abundance of similar states across spots (may worsen rare niches).

Gradients (all consistent with the loss, so backtracking guarantees monotone
descent): `∇recon = 2(ΘR−Yn)Rᵀ`, `∇spatial = 2λ_spatial L_sym Θ`,
`∇state = 2λ_state Θ M` (M = A_state or L_state), `∇sparse = λ_sparse(log Θ + 1)`.
Negative-entropy sparsity (not hard thresholding) avoids destroying rare states;
optional `rare_protection` floors protected states after each projection. The state
graph is **reused** from `state_similarity_regularization.py` (no duplication).

Integration is **Option 2**: the production `SpatCARModel` runs first to produce a
warm-start `Θ` (production solver untouched); the experimental optimiser then
minimises the joint objective from that init. The final estimates come from the
joint optimiser — this is genuine in-objective optimisation, not post-fit mass
redistribution. Gated on `cfg.state_regularized_solver.enabled` (default False).

### 7. Fallback plan
If the joint optimiser proves unstable, fall back to recording diagnostics only and
keep the production estimates (the pipeline wraps the call in try/except and never
breaks a run). The post-fit module from the previous task remains available.

### 8. Benchmark plan
breast + lung synthetic gold-truth, ≥5 seeds; compare spatial baselines + post-fit
state refinement + the new solver modes (competition low/med, laplacian low/med,
competition+sparsity, competition+weak-spatial) on a small λ grid. Promotion gates
centred on conditional within-family RMSE, effective-N, spillover, FP rate, rare
niche; CARD/RCTD/c2l reused qualitatively. Outputs under
`benchmarks/outputs/in_solver_state_regularization/` (gitignored).

**Feasibility decision: do NOT modify the production solver; implement a separate
experimental projected-gradient solver. Proceeding.**

---

## Results (Phases 4–7)

**Headline: a decisive negative result.** Putting the state penalty *inside* the
objective does not fix what post-fit redistribution could not. The separate
projected-gradient solver is numerically sound (monotone, converges, 0 failures)
but converges to **worse** deconvolution solutions than the production NB-CAR
solver: the Gaussian least-squares reconstruction surrogate is a weaker fit than
the NB likelihood, and the state penalties trade one pathology for another —
**competition over-concentrates** (effective-N collapses) and **laplacian
over-smooths** (effective-N inflates). No configuration passes the promotion
gates. Default unchanged; module retained experimental, **not recommended**.

### 1–3. Motivation / why post-fit failed / why Option A
Post-fit state redistribution (previous task) could not recover fine-state mixtures
once the solver had produced a diffuse solution. Option A tests the stronger
hypothesis: regularize *inside* the objective. Implemented as a separate solver
(production NB-CAR untouched), warm-started from the production fit.

### 4–7. Implementation / state graph / penalties / sparsity
Joint objective `||Yn − ΘR||² + λ_spatial·tr(ΘᵀL_symΘ) + λ_state·P_state + λ_sparse·H`,
minimised by backtracking **projected gradient** (simplex per row, or per-family
mass if `preserve_broad_mass`). `P_state` = **competition** `tr(Θ A_state Θᵀ)` or
**laplacian** `tr(Θ L_state Θᵀ)`. Sparsity = entropy `H` (smooth; no hard
thresholding). State graph **reused** from the post-fit module. Optional
`rare_protection` floors protected states.

### 8. Datasets / scenarios
breast + HLCA/lung synthetic gold-truth, 5 seeds, structured scenario (boundary +
gradient + rare niche + mixed + collinear). 130 fits, **0 failures**.

### 9. Configurations (13)
5 spatial baselines (default/no/weak/edge_aware/combined) + post-fit
state_regularized + adaptive_resolution + 6 in-solver modes (competition low/medium,
laplacian low/medium, competition+sparse, competition+weak-spatial).

### 10. Optimization diagnostics
Loss **monotone in every run**; converges in ~60–150 iters (competition+sparse the
exception: over-sparsifies and stops early, converged in 40% of lung seeds). 0
failures, runtime 6–9 s/fit. The optimiser is correct — the problem is the
*objective*, not the optimisation.

### 11–16. Results vs baselines / post-fit / external
The in-solver solver **regresses fine Pearson on both datasets vs default**
(breast 0.859 → 0.60–0.71; lung 0.571 → 0.24–0.61) — the LS surrogate underfits
relative to the NB likelihood. It does not beat weak_smoothing, edge_aware, or
combined on their strengths, and does not improve over the post-fit module's
(already negative) conditional-RMSE result except in one tissue-specific case
(below). CARD/RCTD/c2l reused qualitatively — no change to standing; no regression
introduced to production.

### 17. Conditional within-family RMSE (target metric)
| dataset | default | competition (low/med) | laplacian (low/med) |
|---|---|---|---|
| breast | 0.133 | 0.275 / 0.310 (worse) | 0.159 / 0.166 (**+19% / +25% worse**) |
| lung | 0.276 | 0.320 / 0.389 (worse) | **0.251 / 0.229 (−9% / −17% better)** |

Laplacian improves lung conditional RMSE but **regresses breast by ~19–25%** — so
the cross-dataset primary gate (improve ≥5% in one AND not worse >2% in the other)
**fails**. Competition worsens conditional RMSE everywhere.

### 18. Effective-N and entropy
**Competition over-concentrates**: effective-N collapses far below truth (breast
17→1.2–8.3; lung 12→1.4–16). **Laplacian over-smooths**: effective-N inflates above
truth (breast 17→20–22; lung 12→25). Neither moves effective-N closer to truth in
*both* datasets. The mechanism cannot be calibrated tissue-robustly.

### 19. Spillover and false positives
**Worsened**, not improved: laplacian raises lung spillover (0.220→0.230–0.238) and
FP subtype rate (0.116→0.151–0.166); competition raises lung spillover to
0.31–0.35. Only the degenerate competition+sparse (effective-N≈1) lowers spillover,
by collapsing onto one state.

### 20. Rare-niche behaviour
Rare-niche *sensitivity* reads 1.00 for all solver configs — but this is an
**artifact** of over-mixing (laplacian) / rare_protection floors plus rising false
positives (section 19), not a genuine gain. Rare *precision* drops.

### 21. Fine/broad/local accuracy
Fine Pearson regresses on both datasets (section 11); local RMSE worsens
(breast 0.024→0.035–0.145; lung 0.039→0.044–0.125); broad Pearson regresses for
competition (breast 0.935→0.57–0.79). Laplacian preserves broad better but at the
cost of inflated effective-N.

### 22. Oversmoothing and boundary
Competition lowers oversmoothing only by degenerate concentration; laplacian
*raises* oversmoothing on breast (1.80→2.1–2.3). Boundary F1 ≤ weak in all cases.

### 23. Runtime
6–9 s/fit (comparable to default); 0 failures / 130 fits.

### 24. Promotion gates
**No configuration passes all gates on either dataset.** Best: lung
`solver_laplacian_low` 9/13 (passes cond-RMSE + effN-closer on lung) — but it fails
on breast (cond-RMSE +19%, fine accuracy, effN inflation) and raises lung
spillover, so the cross-dataset requirement fails. Competition configs 5–7/13
(fail cond-RMSE, fine accuracy, effN, spillover). The cross-dataset / recommend /
promote criteria are **not** met by any config.

### 25. Supported claims
- The separate experimental solver is numerically correct: monotone loss,
  convergence, valid compositions, exact per-family mass conservation when enabled,
  0 failures, no runtime blow-up — and the production NB-CAR solver is untouched.
- On **lung only**, the laplacian penalty reduces conditional within-family RMSE
  (−9% to −17%), at the cost of inflated effective-N and higher spillover/FP.

### 26. Unsupported claims
- ✗ "In-solver state regularization improves conditional within-family recovery."
  (Tissue-specific at best; regresses breast.)
- ✗ "It improves effective-N." (Competition collapses it; laplacian inflates it.)
- ✗ "It reduces spillover / false positives." (Worsened.)
- ✗ "It preserves fine accuracy." (Regressed on both datasets.)
- ✗ "Competition is the better penalty." (It is uniformly worse than laplacian here.)
- ✗ "It solves collinearity." (It does not — conditional RMSE/spillover do not both
  improve.)
- ✗ Redeconve / Redeconve-equivalent / state-of-the-art / default solver. None claimed.
- ✗ Any change to default behaviour.

### 27. Recommendation
**Reject.** Retain the in-solver solver as **experimental, opt-in, NOT
recommended**. It does not beat weak_smoothing (or the production default) on the
target metrics, regresses fine accuracy, and its only improvement (lung conditional
RMSE via laplacian) is tissue-specific and bought with inflated effective-N and
increased spillover/false positives. Competition over-sparsifies; laplacian
over-smooths. The production NB-CAR solver remains the default; for within-family
ambiguity prefer soft gating / `unresolved_<family>` + the diagnostic state grouping
from the previous task. Default unchanged.

### 28. Remaining limitations
- The Gaussian least-squares reconstruction surrogate underfits relative to the NB
  likelihood; an NB-faithful in-solver state penalty would require the production
  solver rewrite that the stop condition forbids.
- Synthetic gold-truth only; no real within-family ground truth.
- Combined with the post-fit negative result, the evidence indicates within-family
  collinearity is an **identifiability / marker-information limit**, not a
  solver-regularization limit — no penalty on these inputs both improves conditional
  RMSE and avoids effective-N/spillover degradation across tissues.
