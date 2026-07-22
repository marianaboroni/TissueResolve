# Conditional family estimator — results (experimental)

**First positive result on collinear fine-state recovery.** A supervised, family-specific
conditional estimator, trained and validated on **donor-held-out** pseudobulks,
improves conditional within-family RMSE over mean-reference geometry on **6 families
across both datasets**, with reproducibility across seeds and **without increasing rare
FPR** — while the reject option correctly abstains where improvement would cost rare
precision/calibration.

Benchmark: `benchmarks/dev/family_learnability_benchmark.py` (breast + HLCA/lung, 5
seeds, models ridge / elastic_net / pairwise_ridge, feature_mode=combined). Baseline =
mean-reference conditional split (uniform NNLS on the within-family panel). Outputs
gitignored under `benchmarks/outputs/family_learnability/`.

## Best model per family (5-seed means; Δ = RMSE improvement vs mean-NNLS)

| dataset | family | states | best model | baseline→model RMSE | Δ% | seeds↑ | rare FPR Δ | decision |
|---|---|---|---|---|---|---|---|---|
| lung | Alveolar epithelium | 3 | elastic_net | 0.235→0.089 | **+57.8** | 5/5 | −0.07 | **learnable** |
| lung | Airway epithelium | 7 | elastic_net | 0.196→0.084 | **+56.5** | 5/5 | −0.01 | **learnable** |
| lung | Blood vessels | 5 | elastic_net | 0.210→0.094 | +55.8 | 5/5 | +0.22 | diagnostic_only |
| lung | Myeloid | 12 | elastic_net | 0.116→0.075 | **+33.6** | 5/5 | −0.03 | **learnable** |
| breast | Myeloid | 6 | elastic_net | 0.162→0.109 | **+32.4** | 5/5 | −0.16 | **learnable** |
| breast | Epithelial | 3 | ridge | 0.259→0.138 | +32.2 | 4/5 | +0.07 | diagnostic_only |
| lung | Fibroblast lineage | 6 | elastic_net | 0.171→0.114 | **+31.7** | 4/5 | +0.04 | **learnable** |
| lung | Lymphoid | 6 | elastic_net | 0.161→0.111 | +31.5 | 5/5 | +0.23 | diagnostic_only |
| breast | Endothelial | 5 | elastic_net | 0.184→0.135 | +25.4 | 5/5 | +0.16 | diagnostic_only |
| lung | Smooth muscle | 3 | elastic_net | 0.319→0.244 | +17.6 | 2/5 | −0.42 | partially_learnable |
| **breast** | **T/NK** | 10 | elastic_net | 0.120→0.100 | **+14.8** | 4/5 | −0.15 | **learnable** |
| lung | Submucosal Gland | 4 | pairwise_ridge | 0.272→0.263 | +4.0 | 1/5 | +0.40 | diagnostic_only |

## Findings

- **6/12 families are `learnable`** (≥10% RMSE gain, reproduced in ≥4/5 seeds, rare FPR
  not increased): lung Alveolar/Airway epithelium, lung & breast Myeloid, lung
  Fibroblast, **breast T/NK** (the original collinear target, +14.8%, FPR −0.15).
- **The reject option works as intended:** families where RMSE improves but **rare FPR
  rises** (lung Blood vessels +0.22, lung Lymphoid +0.23, breast Endothelial +0.16) are
  correctly demoted to `diagnostic_only` — not promoted. The estimator does **not**
  force resolution where it would cost rare precision.
- **`elastic_net` is the best model** in 9/12 families (sparse, donor-robust). Random
  forest was not used. 0 failures.
- **Why this works where mean-geometry failed:** training across donors and validating
  on held-out donors forces the model onto **donor-robust** combinations — precisely
  the failure mode (donor-variable discriminating genes) that sank discriminative
  reweighting. The supervised map also absorbs the simulation's mixing/dilution.

## Honest scope / limitations

- Baseline is the **mean-reference NNLS conditional split** on the within-family panel
  (the geometry that failed), not the full-pipeline Poisson GLM family split; the
  improvement is over mean-reference geometry, holding features fixed.
- Train and validation pseudobulks are **simulated** from real single cells
  (donor-disjoint). This controls donor variation but **not** real-bulk confounders
  (protocol, true mRNA content, technical batch). Donor-held-out learnability is
  necessary, not sufficient, for real bulk — real-data validation is the next gate.
- Calibration error is low (~0.04–0.13) but CI coverage was not yet implemented.
- Not all families are learnable, and some only "diagnostic" — the feature is
  **family-specific**, not universal.

## Real-bulk-realism gate (cross-platform + technical confounders) — does NOT pass

The clean same-platform learnability above was re-tested under real-bulk confounders
(`benchmarks/dev/family_realbulk_validation.py`, breast, 5 seeds). Donors are
platform-disjoint (40/41), so cross-platform = donor-held-out + true protocol shift.

| gate | family | best Δ RMSE | seeds↑ | rare FPR Δ | decision |
|---|---|---|---|---|---|
| A cross-platform (3′v3→3′v2) | T/NK | +49% | 5/5 | **+0.21** | diagnostic_only |
| A cross-platform | Endothelial | +44% | 5/5 | **+0.28** | diagnostic_only |
| A cross-platform | Myeloid | +33% | 5/5 | +0.16 | diagnostic_only |
| A cross-platform | Epithelial | +3.5% | 2/5 | −0.17 | partially_learnable |
| B confounders (lib/noise/dropout) | Myeloid | +16% | 5/5 | −0.33 | **learnable** |
| B confounders | T/NK | +15% | 5/5 | +0.06 | diagnostic_only |
| B confounders | Endothelial / Epithelial | <8% | ≤2/5 | + | diagnostic_only |

**Finding: the learnability does NOT robustly survive real-bulk realism.**
- **Under a protocol shift (cross-platform), NO family remains `learnable`:** RMSE still
  drops on paper, but the model **inflates the rare-state FPR** (+0.16 to +0.28) and
  T/NK Pearson collapses (0.29–0.43 vs ~0.62 clean) → the reject option correctly
  demotes all to `diagnostic_only`.
- **Under technical confounders, only breast Myeloid survives** as `learnable` (1
  family, 1 dataset); T/NK falls to `diagnostic_only`.
- The reject option worked exactly as designed — it caught the FPR inflation and
  abstained rather than report fine states.
- TCGA-TNBC real bulk: conditional accuracy **deferred** (no fine-state ground truth;
  isolating a family signal needs full deconvolution — not fabricated).

## Verdict: `needs_more_data` (fails the real-bulk gate; do NOT integrate)

Revised from the clean-pseudobulk `promising_experimental` **down** after the
real-bulk-realism gate:

- On **clean, same-platform** donor-held-out pseudobulk: genuinely promising (6
  learnable families, incl. collinear T/NK).
- On **cross-platform / confounded** bulk (the realistic gate): the gain does **not**
  generalize — it survives only for **breast Myeloid under technical noise**, and a
  protocol shift breaks every family (rare-FPR inflation). So the signal is real but
  **protocol/donor-fragile** — exactly the failure mode that sank discriminative
  reweighting, reappearing across platforms.

It improves conditional within-family RMSE on **≥2 families across ≥2 datasets**
(actually 6, incl. the collinear T/NK target) with acceptable FPR and calibration, and
abstains correctly elsewhere. Per the task rules this is **promising_experimental**.

- **Remain experimental: yes.** Opt-in; not wired into the default pipeline; defaults,
  Poisson GLM solver, unresolved mass, soft gating, spatial code all unchanged.
- **Integrate into pipeline now: no.** Needs real-bulk validation, CI coverage, and a
  decision on how the per-family `learnable/diagnostic/unresolved` routing composes
  with the existing soft gating before any integration.
- **Public README: no** (future-work note only).
