# Conditional family estimator — design (experimental)

New hypothesis (after the distribution-aware negatives): collinear fine states may not
be recoverable through **mean-reference geometry**, but a **supervised, family-specific
conditional estimator**, trained and validated on **donor-held-out** pseudobulks, may
recover them *where a donor-robust learnable signal exists*. Fine states are reported
**only when learnability is demonstrated**; otherwise TissueResolve returns
grouped / `unresolved_<family>` mass.

This is fundamentally different from the prior failed levers (NNLS mean-matching,
discriminative reweighting, metacells, covariance): instead of inverting a mean
reference, we **learn** the map `bulk features → conditional within-family
proportions` from simulated mixtures, and we **measure** out-of-donor learnability,
abstaining when it is absent.

## Why this could work where the others failed

- Discriminative reweighting failed because it amplified **donor-variable** genes. A
  supervised model trained *across* train donors and validated on *held-out* donors
  can instead learn donor-robust combinations (and is penalised at validation if it
  relies on donor-specific signal).
- Cell-level classifiers showed some collinear T/NK states are separable → the signal
  exists; the question is whether it **survives bulk mixing and donor shift**, which is
  exactly what donor-held-out pseudobulk validation measures.

## Method

1. **Family mass** comes from existing TissueResolve output (Poisson GLM family-level
   estimate) — unchanged. The estimator only models the **conditional** split within a
   family, scaled to that family mass.
2. **Training data = simulated pseudobulks** from **train donors only**: draw
   conditional proportions `p` over the family's states (Dirichlet at several
   concentrations to span balanced / imbalanced / rare / absent), sample cells
   accordingly, sum counts → bulk; target = `p`.
3. **Features** (selected inside the train fold only): within-family marker panel
   expression; per-state marker **module scores**; **residual** after subtracting the
   pan-family mean signature; combined. (CLR/log-ratio optional.)
4. **Models** (interpretable, no deep learning): `ridge` (default), `elastic_net`
   (multitask), `pairwise_ridge` (one-vs-one contrasts). Random forest diagnostic-only.
5. **Validation = donor-held-out** pseudobulks (val/test donors never seen in
   training). Compare conditional RMSE vs the mean-NNLS baseline.
6. **Reject option / learnability classes:**
   - `learnable` — conditional RMSE improves ≥ τ (default 10%) over the Poisson-GLM /
     mean direct split, reproduced across ≥4/5 seeds, rare FPR not increased.
   - `partially_learnable` — small but positive improvement; report as grouped.
   - `not_learnable` — no improvement → `unresolved_<family>`.
   - `unstable` — high cross-donor/seed variance → abstain.
   - `diagnostic_only` — poor calibration.

Gates are **configurable**, not hardcoded.

## Outputs

Benchmark: `benchmarks/outputs/family_learnability/` (gitignored): summary,
per-family / per-seed / per-scenario metrics, model status, calibration, verdict.json.
Estimator: a fitted object carrying model + selected features + donor-held-out
validation metrics + decision; `predict_conditional_family` returns proportions summing
to `family_mass`, a confidence score, and the decision (with reject).

## Non-negotiables

Opt-in, experimental; **defaults, Poisson GLM solver, unresolved mass, soft gating,
spatial code all unchanged**; not integrated into the default pipeline in this task;
no test-seed tuning; outputs gitignored; nothing committed.
