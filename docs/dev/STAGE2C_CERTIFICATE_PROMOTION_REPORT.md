# Stage 2C — Certificate Promotion to (opt-in) Core Output

## OBJECTIVE
Promote the Stage 2B calibrated identifiability certificate toward a core TissueResolve output:
(1) validate a **fixed default `shift_scale`** so it calibrates without a per-run calibration split;
(2) add **`recommended_merge` guidance** for confounded (GROUP_ONLY/UNRESOLVABLE) clusters;
(3) **wire it into the bulk pipeline as an opt-in reported output that never changes estimates**;
(4) re-confirm the third-tissue status. No new solver / adaptive / cohort work; deconvolution
estimates and all existing defaults are unchanged.

## DEFAULT SHIFT_SCALE VALIDATION
`run_stage2c_default_shift.py` measures held-out predicted-interval coverage (nominal 90%) at a
*fixed* `shift_scale` (no fitting) on breast + lung, full-depth and 1000-gene panel:

| fixed shift_scale | mean held-out coverage (target 0.90) |
|---|---|
| 1.0 | 0.881 (under-covers; lung full 0.82) |
| **2.0** | **0.930** (breast 0.915/0.951, lung 0.890/0.963) |
| 3.0 | 0.950 (over-covers) |

`shift_scale = 2.0` is closest to nominal across both tissues and both regimes (all conditions
0.89–0.96, slightly conservative on panels). Set as `DEFAULT_SHIFT_SCALE = 2.0` — the certificate now
gives ~nominal donor-aware coverage **without requiring the user to hold out a calibration split**.
(Stage 2B's per-run calibration remains available and also lands at ≈2.)

## RECOMMENDED_MERGE
`calibrated_identifiability_certificate` now emits `recommended_merge` per type: for
GROUP_ONLY/UNRESOLVABLE types (individual split not recoverable) it lists the confusable-cluster
members that should be read out as an **aggregate** instead of individually; empty for RESOLVABLE
types. This is the actionable guidance a user needs when a cluster is flagged.

## PIPELINE WIRING (opt-in, estimate-preserving)
- `tissueresolve.bulk_identifiability(bulk, ref, *, shift_scale=None, min_detect_frac=0.10)` — new
  public function. Derives query-detectable genes and library size from the bulk matrix, warns if the
  reference lacks `donor_cv` (donor-level uncertainty then unavailable), and returns the calibrated
  certificate. Uses `DEFAULT_SHIFT_SCALE` when `shift_scale is None`.
- `deconv_bulk(..., identifiability=False, identifiability_shift_scale=None)` — new opt-in. When
  `True`, the certificate is attached to the result's `identifiability` attribute across every
  pipeline path (flat / auto / hierarchical / state-aware / explicit-solver) via a `_finish` wrapper;
  failures degrade to a warning and never break the estimate.
- `BulkPipelineResult.identifiability: Optional[Any] = None` — new optional field (default None →
  fully backward compatible). Exported `bulk_identifiability` from the package root.

**Estimate invariance** is guaranteed by construction (the certificate is computed from the
reference + query genes only, never touching the solve) and asserted by a test that
`deconv.proportions` is byte-identical with and without `identifiability=True`.

## THIRD TISSUE
Still **NOT EXECUTED** — CRC = CMS molecular subtypes (not cell types); no other donor-annotated
cell-typed reference offline. `run_stage2b_third_tissue.py` provides the reproducible drop-in path;
the pipeline is tissue-agnostic.

## TESTS ADDED
`tests/test_bulk_identifiability_integration.py` (5): default off (backward compatible); opt-in
attaches a certificate AND leaves `deconv.proportions` byte-identical; `recommended_merge` populated
for the confounded A;C cluster and empty for RESOLVABLE types; `DEFAULT_SHIFT_SCALE == 2.0` used;
warning when the reference has no `donor_cv`. (Plus the 8 Stage 2B calibration tests.)

## TEST RESULTS
13/13 new (5 integration + 8 calibration) pass. Full default suite: **1524 passed, 1 skipped, 0
failed** (+5 vs Stage 2B; no regressions). Production files touched: `api.py`, `bulk/pipeline.py`,
`__init__.py` — all additive/opt-in; defaults and estimates unchanged.

## DECISIONS (restated)
- **Certificate → B, now shipped as a first-class opt-in core output.** It is calibrated (Stage 2B:
  swap r≈0.99, zero false-merge/false-resolution, ~nominal coverage) and now (a) calibrates with a
  validated fixed default (no user split) and (b) is reachable from the bulk pipeline
  (`deconv_bulk(identifiability=True)` / `bulk_identifiability`). Full **default-on (A)** remains
  gated on a third tissue and on validating `DEFAULT_SHIFT_SCALE` beyond two tissues. Per standing
  instruction, no default was flipped (the output is opt-in and additive).
- **Cohort covariance → E** (unchanged; out of scope).

## LIMITATIONS
`DEFAULT_SHIFT_SCALE=2.0` validated on two tissues only; donor-level uncertainty needs a
donor-built reference (else a warning + counting-noise-only intervals); one solver family / Poisson-NB
noise; third tissue outstanding.

## NEXT RECOMMENDED STAGE
Stage 2D — flip to default-on core once a third tissue confirms `DEFAULT_SHIFT_SCALE` and the
swap-based classes; surface the certificate in `generate_report` (recoverability class + swap risk +
donor-aware interval + recommended_merge per type); optionally let `resolution_mode="auto"` consume
`recommended_merge` to aggregate GROUP_ONLY clusters (a real default change — separate gated stage).
