# Stage 1.6 — adaptive-signature audit (separability & resolvability)

Answers to the §3 questions, verified against the code.

1. **Separability metrics that exist:** `reference/separability.py::compute_separability` —
   pairwise **Bhattacharyya coefficient** (+ Jeffreys divergence, Pearson r) on L1-normalised
   CPM **mean profiles**; `SeparabilityReport` with problematic pairs (BC>threshold).
2. **Duplicate resolvability:** YES — `reference/resolution.py` (`classify_resolvability`,
   scalar-score → class) and `hierarchy.evaluate_within_family_resolvability` (per-family,
   calls `compute_separability` + discriminating-gene + spillover thresholds).
3. **Which is wired to the pipeline:** the **hierarchical** path uses
   `evaluate_within_family_resolvability` (via `assemble_hierarchical_estimates`). The bulk
   **flat/default** path (GeneSelector + wNNLS) does not use either. `classify_resolvability`
   is used by the reference `resolution.py` analysis layer (reporting), not the solver.
4. **Separability computed from:** **mean profiles (reference)** — NOT donor pseudobulks, NOT
   individual cells, NOT the query.
5. **Global or pairwise:** **pairwise** (all cell-type pairs); a per-type summary is derived.
6. **Donor-aware:** **NO** (mean profiles only; ignores cross-donor variability).
7. **Profile distance vs mixture performance:** **profile distance only** (BC on means); it does
   NOT measure real donor-held-out deconvolution/mixture performance.
8. **Closest confounders:** highest-BC pair per type (from the pairwise matrix).
9. **Spillover:** `benchmark/spillover.py` — simulate mixtures, estimate a spillover matrix
   (mixture-based, not from mean profiles).
10. **Varies with #genes:** YES — BC is computed on whatever gene set is passed, so separability
    is panel-dependent (a lever for adaptive allocation, but also a confound if uncontrolled).
11. **Data leakage:** `compute_separability` itself uses only the reference; leakage risk is in
    how it's used (must be train-donor-only when guiding budget).
12. **Metric that can guide budget without the test set:** donor-held-out (train-donor)
    **incremental deconvolution gain** per type + pairwise BC to the closest confounder +
    donor stability — none of which need the final query truth.
13. **Current budget application:** `GeneConfig.n_genes` (fixed 500) for GeneSelector;
    `fine_global` = `top_n_per_type` union or a stratified round-robin cap.
14. **Order-favouring truncation:** the flat `genes[:N]` bug was fixed in 1.5 (stratified
    round-robin); GeneSelector caps by composite score (not positional).
15. **Duplicate genes across populations:** deduped by `dict.fromkeys`; a gene shared by two
    types counts once — adaptive allocation must redistribute the freed budget.
16. **Coverage of eligible types:** stratified round-robin guarantees ≥ min per eligible type.
17. **Non-resolvable types:** hierarchical gating routes them to `unresolved_<family>`; the
    reference resolution layer classifies them; the flat path does not abstain.
18. **No-donor references:** 1.5 fallback chain (donor→batch→pooled→inadequate) with provenance +
    capped status.
19. **Query gene availability limiting the decision:** NOT currently used in selection (genes are
    intersected at solve time, but selection does not plan around query retention).
20. **Reuse vs consolidate:** REUSE `compute_separability`, `donor_aware_de`/`select_donor_aware_genes`,
    the donor-held-out CV harness from 1.5, and the Poisson solver. CONSOLIDATE the two
    resolvability layers eventually (do NOT add a third separability implementation — Stage 1.6
    adds an *allocator* on top of the existing separability + donor-held-out gain, not a new
    distance metric).

## Consequence for Stage 1.6
Existing separability is **mean-profile, pairwise, not donor-/deconvolution-aware** — insufficient
as the sole budget driver. The allocator will combine (a) pairwise BC to the closest confounder
(reuse) with (b) a **donor-held-out per-type incremental deconvolution-gain curve** (the
leakage-safe, deconvolution-aware signal §12 asks for), and stop at diminishing returns. No third
separability metric is created; no default changed.
