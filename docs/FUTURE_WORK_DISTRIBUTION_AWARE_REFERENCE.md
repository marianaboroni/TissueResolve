# Future work — distribution-aware reference modeling (P2a/b)

**This is future research and is not part of the alpha release.** No implementation is
included. It extends the design in `docs/dev/DISTRIBUTION_AWARE_REFERENCE_DESIGN.md` into a
prioritized roadmap, motivated by the alpha's central negative result: with a
**mean-profile** reference, collinear fine states (Bhattacharyya > 0.97) are not
identifiable, and no regularization recovers them. Raising that ceiling needs *more
information* in the reference, not more regularization.

## Directions

1. **Donor × state metacells.** Replace single mean profiles with a few robust
   metacell profiles per (donor, state) that preserve within-type heterogeneity and
   donor structure. Deconvolve against the metacell dictionary, then aggregate to
   states where the mapping is identifiable.
2. **Robust centroid estimation.** Outlier-resistant / shrinkage centroids per state,
   weighted by the reference-reliability score already computed
   (`experimental/reference_uncertainty.py`), so noisy rare profiles are down-weighted.
3. **Covariance-aware signatures.** Model per-type gene–gene covariance (low-rank
   `Σ_k`), not just the mean `μ_k`. States with equal means but different co-expression
   become separable under a likelihood that scores the joint pattern — the one route
   that directly attacks `BC(μ)` collinearity.
4. **Within-family NMF / topic programs.** Factor each family's cells into a few
   non-negative programs (better conditioned than the collinear subtype basis),
   deconvolve programs, and map to states where identifiable (else report the program/
   group). Built from existing NMF; no Gibbs sampler, no deep learning.
5. **Bayesian reference uncertainty.** Carry a posterior over `μ_k` (width ∝ cells/
   donors) and propagate it so unreliable rare estimates widen rather than mislead;
   a cheap empirical-Bayes / reliability-weighted GLM is the pragmatic first step,
   composing with the Poisson GLM.
6. **Improved dispersion estimation.** Per-cell-type / shrunk donor-level dispersion.
   (Alpha evidence: donor-level dispersion did **not** make NB beat Poisson; revisit
   only if a regime with new evidence appears.)

## How these address collinear fine states

§3 and §4 add information the mean discards (covariance / program structure) and could,
in principle, raise the identifiability ceiling — to be **proven** by conditional
within-family RMSE, not assumed. §1, §2, §5 mainly improve calibration, robustness, and
honest resolution (rare types, donor effects, uncertainty), composing with the existing
adaptive-resolution reporting rather than forcing splits.

## Risks

Higher-moment / metacell / program models add collinear columns → must lean on
abstention (soft gating, `unresolved_<family>`) and adaptive grouping, never silent
merges. Reference reliability becomes a hard dependency for rare types. Compute and
reference-build complexity rise → keep all routes opt-in and benchmark-gated. No deep
learning; not Redeconve; no equivalence claims.

## Benchmark plan

Same gold-truth pseudobulk + synthetic spatial, breast + lung, ≥5 donor-disjoint seeds,
with a **collinear-pair scenario** as the primary stressor. Primary metric: conditional
within-family RMSE (the metric all prior negatives turned on), plus effective-N,
spillover, rare precision/recall, and uncertainty calibration (ECE/Brier). A route is
retained only if it improves conditional within-family recovery without degrading rare
niches/spillover, replicated across seeds and both tissues. Start with **§3 (covariance,
low-rank)** and **§5 (reliability-weighted GLM)** as the highest-value, composable bets.
