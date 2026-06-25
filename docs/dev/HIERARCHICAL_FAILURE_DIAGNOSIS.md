# Hierarchical-Failure Diagnosis (Phase 1)

**Mandate.** Diagnose *why* the hierarchical path underperforms the flat path on
bulk and spatial ground truth, **before** changing any algorithm. Diagnosis only
— no algorithm was added or modified. Baseline frozen in
`benchmarks/baselines/current_baseline/`.

**Method.** Four oracle experiments (`benchmarks/diagnostics/hierarchical_oracle.py`,
deterministic) on held-out-donor synthetic data with known broad+fine mRNA truth
(bulk: 3 seeds × {imbalanced, similar_subtypes}; spatial: 2 seeds × λ-sweep).
Raw table: `benchmarks/outputs/hierarchical_oracle_experiments.tsv`. Reference
built from reference donors only; oracle gene panels derived from the reference
(not test donors) — rule 11 respected.

---

## Bulk — mean fine Pearson over seeds × scenarios

| Variant | fine Pearson | within-family conditional r | pred eff-N (true 10.2) | complexity err↓ |
|---|---|---|---|---|
| flat_nnls (baseline) | **0.852** | 0.641 | 9.96 | 0.68 |
| hier_oracle_broad (true family mass + est. conditional, ungated) | 0.707 | **0.755** | 9.73 | 0.75 |
| hier_pregating (all-subtype combine, no gating) | 0.609 | 0.571 | 9.97 | 1.03 |
| hier_no_gating (`allow_unresolved=False`) | 0.609 | 0.571 | 9.97 | 1.03 |
| hier_oracle_panel_nogating (within-family DE panel, ungated) | 0.609 | 0.571 | 9.97 | 1.03 |
| **hier_standard (current gated default)** | **0.128** | 0.317 | **2.55** | **7.62** |

`hier_standard` parks **70–86% of mass as unresolved** (5 families gated every
run), collapsing effective-N from the true ~10.2 to **2.5**.

## Bulk — dominant failure modes (ranked by evidence)

1. **Confidence gating / unresolved-mass allocation is the PRIMARY failure
   (Experiment 3).** Removing gating lifts fine Pearson from **0.128 → 0.609**
   and restores effective-N from 2.5 → ~10.0 (complexity error 7.62 → 1.03). The
   gate marks all 5 multi-member families unresolved and dumps 70–86% of mass
   into `unresolved_*`, which is both the accuracy collapse and the
   "over-concentration / 4–5 populations" symptom. **This single component
   accounts for the majority of the hierarchical deficit.**
2. **Broad-stage mass error is the SECONDARY failure (Experiment 1).** With gating
   removed, injecting the *true* family mass lifts fine Pearson **0.609 → 0.707**
   and within-family conditional r **0.571 → 0.755**. So broad→fine mass error
   propagates, but it is a smaller effect than gating.
3. **Within-family gene selection is NOT a bottleneck (Experiment 2).** An oracle
   within-family DE panel (top discriminative genes from the reference) gave
   **identical** accuracy to the default panel (0.609 = 0.609). The fine
   signatures / panel are *not* the limiting factor for this failure.
4. **Residual factorization gap.** Even the best hierarchical variant
   (oracle-broad, ungated: 0.707) still trails flat (0.852). The sequential
   broad→conditional factorization is inherently less accurate than the flat
   joint solve, independent of gating and broad error.

## Spatial — λ (smoothing) sweep, mean over 2 seeds

| λ_spatial | fine Pearson | oversmoothing (1=ideal) | local RMSE↓ | domain ARI |
|---|---|---|---|---|
| 0.0 (none) | 0.689 | 0.32 (under) | 0.042 | 0.219 |
| **0.02 (weak)** | **0.793** | **1.39** | **0.029** | **0.250** |
| 0.1 (current default) | 0.775 | 1.78 | 0.030 | 0.241 |
| 0.5 (strong) | 0.655 | 1.49 | 0.014 ARI | 0.035 |

## Spatial — dominant failure mode (Experiment 4)

5. **The default λ=0.1 over-smooths; a weaker λ≈0.02 dominates it** — higher fine
   Pearson (0.793 vs 0.775), less over-smoothing (1.39 vs 1.78), equal/better
   local RMSE, and better domain recovery. λ=0 under-smooths (fine 0.689), so
   *some* smoothing helps, but the current default is simply **too strong**. The
   spatial over-smoothing is a **tunable-λ problem**, not evidence that a new
   model is required.

---

## What the diagnosis means for the proposed v2 changes

The evidence **reprioritises** the committee's proposals and, per rules 16–17
(prefer the simpler method; keep unsupported methods experimental):

| Proposed change | Justified by diagnosis? | Priority |
|---|---|---|
| **Partial confidence-weighted unresolved mass** (Ch. 1.5) — replace all-or-nothing gating | **YES — strongest lever** (0.128→0.609) | **Phase 2, first** |
| **Soft joint broad–fine reconciliation** (Ch. 1.2) — reduce broad-error propagation | **YES — secondary lever** (0.609→0.707 with oracle broad) | **Phase 2, second** |
| **Weaker / level-specific spatial λ** (subset of Ch. 2.2) | **YES — and simple** (λ=0.02 already beats default) | **Phase 3, first (cheap)** |
| Contrastive within-family gene selection (Ch. 1.3) | **NO** — oracle DE panel gave identical accuracy | keep experimental / deprioritise |
| Full NB-CAR VI rewrite (Ch. 2.3–2.6) | **NOT YET** — weak-λ already fixes over-smoothing | gate behind: does it beat tuned-λ baseline? |
| RNA-content correction (Ch. 3) | orthogonal (identifiability, not the hierarchical failure) | independent track |
| Distribution alignment (Ch. 4) | orthogonal (cross-modality) | independent track |

**Headline.** The hierarchical deficit is **mostly a post-processing (gating)
artifact**, secondarily a **broad-mass-coupling** issue — both addressable with
the *cheap* changes (recalibrated partial confidence + soft joint reconciliation).
The *expensive* proposals (contrastive gene modules, NB-CAR VI) are **not
supported as priorities** by this diagnosis and must clear a "beats the cheap
fix" gate before any investment.

## Limitations of this diagnosis

- One tissue, 3 bulk seeds × 2 scenarios, 2 spatial seeds — directional, not
  CI-grade (Phase 2 adds the full 25-replicate statistics).
- "Oracle" gene panel = reference within-family DE (a strong but imperfect proxy
  for truly discriminative genes); a different oracle could change Experiment 2.
- Spatial λ sweep used the global (non-level-specific) λ; level-specific λ
  (Ch. 2.2) is untested and is the recommended Phase-3 experiment.
- Gating is evaluated by toggling `allow_unresolved`; the *calibrated partial*
  alternative is not yet implemented (that is the Phase-2 prototype).

---

*Generated by `benchmarks/diagnostics/hierarchical_oracle.py`. Source:
`benchmarks/outputs/hierarchical_oracle_experiments.tsv`. No core algorithm
modified; nothing committed. Phase 1 is complete — do not proceed to Phase 2
without committee sign-off.*
