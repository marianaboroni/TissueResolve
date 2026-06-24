# Distribution-Aware Reference Modeling — Design (P2a/b)

**Design document only — no implementation in this task.** Motivated by the
accumulated evidence that collinear fine states are an *identifiability /
information* limit: with a **mean-only** reference (`μ_k`), pairs with Bhattacharyya
coefficient `BC > 0.97` are not separable, and neither post-fit nor in-solver
state-similarity regularization recovered them (both prior reports are negative).
The only way to raise the ceiling is to give the model **more information** than the
per-type mean. This document scopes four routes, their risks, and a benchmark plan.

## 1. Covariance-aware signatures
**Idea.** Model each type not only by its mean `μ_k` but by its gene–gene
second-order structure `Σ_k` (or a low-rank approximation). Two states with nearly
equal marginal means but different covariance (co-expression) become distinguishable
under a likelihood that scores the *joint* pattern, e.g. a per-type Gaussian/NB
factor model `x | k ~ (μ_k, U_k U_kᵀ + D_k)` with low-rank `U_k` (r ≈ 5–20).
**Use in deconvolution.** Bulk `b` is a *sum* of cell contributions, so the spot/sample
covariance is `Σ_k θ_k Σ_k` plus mixing terms; a moment-matching or composite
likelihood on (mean, covariance) lets `θ` exploit covariance differences.
**Why it can beat the ceiling.** It adds an orthogonal information channel (2nd
moment) that mean-matching discards — the only one of the four that directly attacks
`BC(μ)`-collinearity.
**Risks.** Covariance estimation is high-variance for rare types (needs P2c
reliability gating); bulk loses some covariance information by summation
(identifiability gain is partial, must be measured, not assumed); heavier compute.

## 2. Metacells by donor × state
**Idea.** Replace single mean profiles with a small set of **metacell** profiles per
(donor × state) — robust pseudo-cells that preserve within-type heterogeneity and
donor structure. Deconvolve against the metacell dictionary, then aggregate metacell
weights back to states.
**Why it helps.** Captures within-type variation and donor effects the single mean
erases; pairs separable on *some* metacells become partially identifiable; integrates
naturally with P2c (donor-aware reliability) and the NB-GLM (P1).
**Risks.** Larger dictionary → more collinear columns (κ ↑) → needs the abstention /
adaptive-resolution machinery (P4) and possibly grouped reporting; metacell
construction adds a reference-build dependency.

## 3. Within-family NMF / topic programs
**Idea.** Within each broad family, factor the cell-level expression into a few
non-negative **programs** (NMF / topic model). Deconvolve programs (well-conditioned,
fewer than the collinear subtypes), then map program loadings → states where the
mapping is identifiable; where it is not, report the program/group (ties into P4).
**Why it helps.** Programs are chosen for conditioning, sidestepping the collinear
subtype basis; this is the BayesPrism-style state↔type idea built from existing
non-negative factorization, no Gibbs sampler, no DL.
**Risks.** Program→state mapping can itself be ambiguous (must be reported, not
forced); over-/under-factorization (choosing #programs) needs validation; programs
are less directly interpretable than named subtypes.

## 4. Bayesian reference uncertainty
**Idea.** Carry the **posterior** of each `μ_k` (and dispersion) — e.g. a
Normal/Gamma posterior whose width scales with cells/donors (P2c already estimates
the inputs) — and propagate it into deconvolution (hierarchical model or
uncertainty-weighted likelihood) so unreliable rare profiles are automatically
down-weighted and their estimates widen rather than mislead.
**Why it helps.** Turns P2c's reliability into calibrated *estimate* uncertainty;
improves rare-type calibration and honest abstention rather than raw accuracy.
**Risks.** Full Bayesian inference is heavier; a cheap empirical-Bayes /
reliability-weighted GLM is the pragmatic first step (composes with P1).

## 5. How these address collinear states
- §1 and §3 add information the mean lacks (covariance / program structure) → can, in
  principle, raise the `BC(μ)` ceiling — to be proven by **conditional within-family
  RMSE**, not assumed.
- §2 and §4 mostly improve *calibration, robustness, and honest resolution* (rare
  types, donor effects, uncertainty) rather than forcing a split — they raise
  *effective* precision and pair with P4 (adaptive resolution) and P2c.

## 6. Expected risks (cross-cutting)
- Higher-moment / metacell / program models add collinear columns → must lean on
  abstention (soft gating, `unresolved_<family>`) and P4 grouping, never silent merges.
- Reference reliability (P2c) becomes a hard dependency for the rare-type cases.
- Compute and reference-build complexity rise; keep all routes opt-in and gated.
- No deep learning; no Redeconve code/claims; no equivalence claims.

## 7. Benchmark plan
Same gold-truth pseudobulk (bulk) + synthetic spatial gold-truth, breast + lung, ≥5
donor-held-out seeds, with a **collinear-pair scenario** as the primary stressor.
Primary metric: **conditional within-family RMSE** (the metric all prior negatives
turned on), plus effective-N, spillover, rare precision/recall, and uncertainty
calibration (ECE/Brier from `soft_hierarchy/confidence.py`). Promotion gates mirror
the NB-GLM gates; a route is retained only if it improves conditional within-family
recovery without degrading rare niches / spillover, replicated across seeds and both
tissues. Start with **§1 (covariance, low-rank)** as the highest-information bet and
**§4 (reliability-weighted GLM)** as the cheapest composable win on top of P1.
