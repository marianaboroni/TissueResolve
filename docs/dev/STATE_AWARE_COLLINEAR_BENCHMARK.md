# State-aware 3-level vs 2-level hierarchical — collinear-type recovery

Experimental test (negative result). Question: does the opt-in **state-aware**
3-level deconvolution (`broad → cell type → state`) improve prediction of
**collinear fine types** vs the standard 2-level hierarchical path?

Setup (so collinear subtypes become "states" inside a well-separated "cell type"):
`state = fine cell type`, `cell type = broad family` (T/NK, Myeloid, Endothelial…),
`broad = super-compartment` (Immune / Stromal / Epithelial). Breast pseudobulk,
donor-disjoint, 5 seeds × 3 scenarios (similar_subtypes / imbalanced / rare).
Script: `benchmarks/dev/state_aware_collinear_benchmark.py` (outputs gitignored).

## Result — no improvement on collinear types

**Conditional within-family RMSE for the collinear T/NK family (11 subtypes) is
byte-identical between the two methods in every seed × scenario (max |Δ| = 0.0000).**

| metric (5-seed mean) | 2-level hierarchical | state-aware 3-level |
|---|---|---|
| conditional within-family RMSE (all) | 0.275 | 0.262 |
| **cond-RMSE [T/NK] (collinear)** | **identical (Δ=0.000)** | **identical (Δ=0.000)** |
| cond-RMSE [Myeloid] (collinear) | identical | identical |
| fine Pearson | 0.245 | 0.352 |
| rare precision | 0.015 | 0.001 (worse) |
| pairwise spillover | 0.0056 | 0.0013 |
| unresolved mass | ~0.80 | ~0.97 (abstains more) |

## Interpretation

- **The collinear split is unchanged by design.** State-aware Stage 2
  (`cell type → state`) for a family uses the *same* within-family resolution
  machinery (`assemble_hierarchical_estimates` + soft gating + within-family panels)
  as the 2-level family→fine step. So for T/NK, Myeloid, etc. the conditional
  `P(state | family)` estimate is the *identical* computation — hence Δ = 0.0.
- **What state-aware changes is coarser-level allocation and abstention**, not the
  collinear resolution: it adds a compartment level and pushes *more* mass to
  `unresolved_*` (0.80 → 0.97). The small overall conditional-RMSE gain (0.275 →
  0.262) and the higher fine Pearson come from a **non-collinear** family
  (Epithelial) and better broad-level allocation — **not** from resolving collinear
  states. Rare precision actually drops (heavier abstention on the collinear Treg).

## Conclusion

**State-aware deconvolution does not improve prediction of collinear fine types.**
This is the same identifiability ceiling seen across the project: with a mean-profile
reference, collinear states (BC > 0.97) are not identifiable, and adding hierarchy
levels (like adding regularization) does not create the missing information — it only
shifts where mass is abstained. State-aware remains a useful *structuring/abstention*
tool (correct unresolved mass at the right level), **not** a fix for collinearity.
Recommended path for collinear states is unchanged: interpret them at the family/
group level (adaptive resolution) or pursue distribution-aware reference modeling
(`docs/FUTURE_WORK_DISTRIBUTION_AWARE_REFERENCE.md`), not more hierarchy/regularization.
Default behaviour unchanged; state-aware stays experimental/opt-in.
