# TissueResolve v2 — Technical Design (experimental, diagnosis-driven)

Status: **DESIGN ONLY.** Nothing here is implemented or wired into defaults.
Prioritisation follows the Phase-1 diagnosis (`HIERARCHICAL_FAILURE_DIAGNOSIS.md`)
and rules 1–17 (notably 3: no default change before all gates pass; 5: feature
flags; 7: no strong L1 by default; 16: prefer the simpler method; 17: keep
unsupported methods experimental).

## Target claim (not yet earned)
"TissueResolve separates shared lineage programs from subtype-specific contrasts,
jointly estimates broad and fine RNA-derived composition, and reports only the
fraction of subtype signal supported by stable, donor-consistent evidence."
Forbidden until the promotion gates in `TISSUERESOLVE_V2_VALIDATION_PLAN.md` pass.

## Architecture
All v2 code is isolated and **not imported by default**:
```
src/tissueresolve/experimental/
  soft_hierarchy/        # Challenge 1 (PRIORITY)
  rna_content/           # Challenge 3 (independent)
  distribution_alignment/# Challenge 4 (independent)
  nbcar_vi/              # Challenge 2 (GATED — only if cheap spatial fix fails)
```
Each exposes one minimal entry point (`fit_soft_hierarchy(...)`,
`estimate_rna_content_factors(...)`, `fit_distribution_alignment(...)`,
`fit_nbcar_vi(...)`), typed dataclasses (Config/Fitted/Diagnostics/Prediction),
and records: algo version, feature status, config, seed, selected genes, device,
convergence, runtime, memory, warnings, estimate type.

---

## Challenge 1 — Soft-hierarchical model (PRIORITY; diagnosis-supported)

The diagnosis shows the failure is **gating (primary)** + **broad-mass coupling
(secondary)**, NOT gene selection. So v2 Challenge 1 is implemented in this order,
each behind a flag, **ablated individually before combination** (rule 15, Phase-2
rule "do not merge before measuring").

### 1.5 Partial confidence-weighted unresolved mass — **FIRST** (biggest lever)
Replace the all-or-nothing family gate (which parked 70–86% of mass as
unresolved, collapsing eff-N 10→2.5) with calibrated partial resolution:
- per-subtype calibrated confidence `c_{k,s} ∈ [0,1]`;
- `θ_resolved = θ_k · c_{k,s}`, `u_f = π_f − Σ_k θ_resolved`;
- confidence integrates bootstrap stability, gene-masking stability, within-family
  separability, pairwise marker support, query overlap, panel-perturbation
  sensitivity, optimisation uncertainty;
- **calibrated on held-out synthetic mixtures** (never the final test set).
Minimum baseline to beat: ungated combine (fine Pearson 0.609). Target: ≥ flat
(0.852) while retaining honest abstention only where truly non-separable.

### 1.2 Soft joint broad–fine reconciliation — SECOND (secondary lever)
Replace frozen sequential (broad → freeze → fine) with soft consistency
`π_{f,s} ≈ Σ_{k∈f} θ_{k,s}`. Two formulations, **benchmarked against each other
and against the simpler one wins** (no assumption that probabilistic > constrained):
- **A (constrained, default first):** `min_{π,θ≥0} L_recon + λ_h‖π − Aθ‖² + λ_r R(θ)`,
  simplex constraints; `R` is a *mild* ridge/entropy regulariser — **NOT strong L1**
  (rule 7; predictions are already too sparse).
- **B (probabilistic):** LogisticNormal hierarchy `π ~ LN(μ_π,Σ_π)`,
  `q_f ~ LN(μ_f,Σ_f)`, `θ_k=π_{f(k)}q_k`. Only pursued if A underperforms.
Identifiability: subtype contrasts centred within family,
`Σ_{k∈f} ω_k D_{gk}=0` (ω = subtype prevalence weights) so `B` (shared lineage)
and `D` (contrast) are separable; document the chosen centering.

### 1.1 Shared lineage + contrast signatures — supporting representation
`S_{gk}=B_{g,f(k)}+D_{gk}`. Enables the soft hierarchy and confidence (contrast
magnitude ↔ separability), but per the diagnosis is **not itself the bottleneck**.

### 1.3 Contrastive within-family gene selection — **DEPRIORITISED**
Diagnosis E2: an oracle within-family DE panel gave *identical* accuracy to the
default. Kept experimental; implemented only if 1.5+1.2 plateau below flat and an
ablation shows panel-limited residual error. No strong-L1 selection.

### 1.4 Evidence modules — folded into 1.5's confidence (multi-gene, cross-donor,
mask-stable support) rather than a separate component.

## Challenge 2 — Spatial (CHEAP FIX FIRST; VI gated)
Diagnosis E4: default λ=0.1 over-smooths; λ≈0.02 already beats it on accuracy +
calibration. Therefore:
- **2.a (first, cheap):** expose + re-tune λ_spatial; add **level-specific λ**
  (`λ_broad > λ_fine`) via the existing NB-CAR — likely a config/loop change, not
  a new model. Must beat the tuned-global-λ baseline on held-out validation.
- **2.b edge-preserving adaptive weights** (expression/H&E kernels): only if
  level-specific λ leaves boundary blurring.
- **2.c NB-CAR VI prototype** (`experimental/nbcar_vi/`): **GATED** — built only
  if 2.a/2.b fail to fix over-smoothing without losing accuracy. Backend: audit
  shows the stack is **NumPy/scipy/sklearn (no PyTorch/JAX in deps)**; a VI
  prototype would add a heavy dependency, so it must clear a strict
  cost/benefit gate (materially better fine accuracy + boundary at acceptable
  runtime) before adoption. Do not claim MPS/GPU unless measured.

## Challenge 3 — RNA-content → cell fraction (independent track)
Identifiability framing: `θ̂` (RNA proportion) estimated; `P_c` (cell fraction)
identifiable only with RNA-content factors `s_c`. Estimate `s_c` via (1) external
measurements, (2) reference proxies with protocol correction, (3) empirical-Bayes
with donor/platform random effects + intervals, (4) sensitivity mode. Output BOTH
`rna_proportion` and `cell_fraction_estimate` (distinct fields, never overwrite);
if `s_c` unreliable, emit RNA proportions only and state cell fractions are not
identifiable (rules 8). Behind a flag; off by default.

## Challenge 4 — Distribution alignment (independent track)
Conservative ladder (rule: simplest first): protocol-aware gene weighting →
shared low-rank factorisation → query-adaptive reweighting → regularised OT only
if needed. Nested CV (train donors select genes, val tune thresholds, test
evaluate); held-out-donor overcorrection checks; never use test-donor info
(rule 11).

## Cross-cutting requirements (rule 15, every component)
isolated implementation · unit tests · synthetic recovery tests · runtime+memory
profiling · ablation · docs · feature-status label.
