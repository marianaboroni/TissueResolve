# Hierarchical gating: audit + learnability-gate (negative result)

**Status: experimental audit; NO change to the default gate.** Etapa 5 of the
Rectangle-driven improvement work. Evidence is donor-held-out pseudobulk (breast + lung),
Poisson solver; single-seed/few-seed probes — directionally clear, not a multi-seed proof.

## 1. Audit — is the hierarchical gate over-abstaining?

The default soft gate decides per-family resolvability from **mean-profile within-family
separability** (`evaluate_within_family_resolvability`: Bhattacharyya, #discriminating
genes, sibling Pearson). Coverage on breast/lung is ~0 (it abstains on every multi-member
family), which looked like severe over-abstention.

Cross-tabbing the gate's resolve/abstain decision against the **achievable** per-family
conditional RMSE (flat Poisson on held-out query donors, 3 seeds):

| tissue | multi-member families | correct abstention | **false abstention** |
|---|---|---|---|
| breast | 6 | 5 (achievable cond-RMSE 0.20–0.30) | **1 — T/NK (0.131)** |
| lung | 9 | 8 (0.15–0.40) | **1 — Myeloid (0.104)** |

**Finding: the gate is *not* broadly over-conservative.** 13/15 families genuinely cannot
be resolved to fine subtypes from these references — abstaining is correct; resolving them
would be false resolution that inflates spillover/FPR. Only ~1 family per tissue (the large
collinear-but-learnable T/NK / Myeloid) is a false abstention. A **blanket gate loosening
would be harmful** (it would convert the 13 correct abstentions into false resolutions).

## 2. Learnability-gate hypothesis (new) — and its falsification

Hypothesis (distinct from the prior failed `conditional_family_estimator` /
`discriminative_within_family`, which built *new estimators*): recalibrate the gate's
**confidence** using an empirical **reference-donor-held-out** recovery of the *standard*
Poisson within-family split — resolve a family only if its subtypes are recovered
out-of-donor on reference-simulated mixtures.

Test: for each breast family, donor-held-out CV (train reference donors → held-out reference
donors) conditional RMSE of the Poisson split, compared to the query-truth achievability
above.

| family | ref-CV cond-RMSE (signal) | query-truth achievable | agree |
|---|---|---|---|
| Endothelial | 0.090 (ranks #1) | 0.200 (poor) | ✗ over-optimistic |
| Myeloid | 0.122 | 0.208 (poor) | ✗ over-optimistic |
| T/NK | 0.132 | 0.131 (good) | ✓ |
| Mural | 0.136 | 0.238 (poor) | ✗ |

**Falsified.** The reference-CV signal does not predict query achievability — it ranks
Endothelial/Myeloid above T/NK, but those are correct abstentions on the query. Gating on
`ref_cv < 0.15` would resolve 4 families, 3 of them false resolutions. This reproduces the
project's central fragility: **clean same-atlas donor-held-out learnability is over-optimistic
and does not transfer** (protocol/mixture shift), the same reason the conditional family
estimator failed its real-bulk gate.

Confound noted: the ref-CV used family-restricted deconvolution (easier than the full-reference
query deconvolution), which inflates optimism uniformly; but the **ranking disagreement**
(Endothelial/Myeloid above T/NK) is the core failure and is not explained by that confound.

## 3. Recommendation

- **Do not modify the default hierarchical gate.** Its abstention is mostly correct; the
  audit shows only ~1 recoverable family per tissue, and no reference-only signal tested here
  can identify it without causing false resolutions.
- The audit metrics (coverage / false-abstention / false-resolution vs achievable
  conditional RMSE) are a useful **diagnostic** and could be surfaced as such, but they
  require query truth so they are a benchmark tool, not a predict-time gate.
- Raising the collinear-family ceiling (T/NK, Myeloid) remains the open hard problem; the
  honest levers are *more information in the reference* (covariance/programs — see
  `docs/FUTURE_WORK_DISTRIBUTION_AWARE_REFERENCE.md`), not gate recalibration.
