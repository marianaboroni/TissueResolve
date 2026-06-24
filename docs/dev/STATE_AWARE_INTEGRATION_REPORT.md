# State-aware integration report (Part 9)

Integration of the granular-signature + state-aware components behind one
explicit experimental flag, with honest reporting. No deferred features added.

## 1. Files created
- `docs/STATE_AWARE_INTEGRATION_PLAN.md`, `docs/FEATURE_STATUS.md`,
  `docs/STATE_AWARE_INTEGRATION_REPORT.md`
- `tests/test_state_aware_integration.py`
- `benchmarks/tests/test_state_aware_benchmark_method.py`
- `examples/real_breast_cancer/tests/test_state_aware_report.py`

## 2. Files modified
- `src/tissueresolve/api.py` — `deconv_bulk(..., state_aware=False, reference_adata,
  state_to_celltype, broad_col, cell_type_col, state_col)`; `_run_state_aware_bulk`
  routes to the 3-level solver with metadata + two-level fallback; gating kwargs
  now override cfg defaults (fixes a duplicate-keyword collision).
- `src/tissueresolve/cli.py` — `--state-aware` flag (default off; only effective
  with `--resolution-mode hierarchical`); records `hierarchy_mode`,
  `state_aware_enabled`, `state_aware_feature_status`, fallback reason in
  `analysis_plan.json`; `_execute_bulk` writes state-aware outputs.
- `benchmarks/bulk/methods/tissueresolve.py` + `shared/method_registry.py` —
  `TissueResolve_state_aware` registered as a distinct internal bulk method.
- `examples/real_breast_cancer/scripts/07_generate_reports.py` — graceful
  "State-aware / multi-granularity deconvolution" report subsection.

## 3. Tests added
api routing (default unchanged → `HierarchicalBulkResult`; `state_aware=True` →
`StateAwareBulkResult` + experimental metadata; no-state-labels → two-level
fallback + `fallback_reason`; with state mapping → 3-level); CLI flag + plan
metadata (incl. fallback when not hierarchical); benchmark method registered +
distinct + internal; report subsection states status + experimental + no-state
note. **28 + 8 new tests pass.**

## 4. pytest result
Affected areas green: state-aware/granular/three-level (28), state-aware report
(8), cli/api/hierarchical/config (88), reference/hierarchy area (57). Full
~35-min suite not re-run this step (recommended as the commit gate).

## 5–8. Integration status
- **Granular signatures generated in the real workflow?** Wired: the api builds
  them from a per-cell `reference_adata` when provided; for the breast-cancer
  reference (only broad + fine labels, AnnData not loaded at report time) the
  panels are **not generated**, and the report says so. So: *capable and wired,
  not exercised on the current real run.*
- **State-aware solver wired into benchmark?** **Yes** —
  `TissueResolve_state_aware` is a distinct registered bulk row.
- **State-aware solver wired into report?** **Yes** — graceful subsection.
- **Two-level fallback works?** **Yes** — verified (no state labels → broad→
  cell-type, mass preserved, `fallback_reason` recorded).

## 9–15. Real-data effect (honest)
**Not measured.** The state-aware solver was not run on the real breast-cancer
reference because (a) it is not the default and (b) that reference has **no
state labels** (only broad + fine), so there is no third level to resolve and no
state-level accuracy to measure. The two-level fallback is, by construction,
equivalent in structure to the existing broad→fine path, so **no change to
fine-level / family-level accuracy or unresolved mass is claimed**. A genuine
before/after requires a reference with real state labels (or a synthetic one) —
deferred.

## 16. Should the feature remain experimental?
**Yes.** It is behind `--state-aware` / `state_aware=True`, labelled experimental
in the report, CLI help, and `FEATURE_STATUS.md`, and is not validated across
datasets.

## 17. Remaining limitations
- No state labels in the current real reference → state level untested on real
  data; only synthetic-data validation exists.
- Granular panels require per-cell AnnData; absent at report time for the saved
  reference, so the real run uses the two-level fallback (clearly stated).
- Reference adaptation and expression reconstruction remain **not implemented**
  (and not claimed).
- Full benchmark run (`benchmarks/run_all.py`) not executed here (slow / needs
  data); the state-aware row will populate when the benchmark is run.

## 18. Safe to commit?
Yes — additive, default behavior unchanged, all affected tests green, deferred
features untouched. Recommend running the full suite as the commit gate.

## 19. Suggested commit message
```
Integrate state-aware deconvolution behind an experimental flag

Expose the 3-level state-aware solver via deconv_bulk(state_aware=True) and a
--state-aware CLI flag (default off; experimental; two-level fallback when no
state labels, with fallback_reason recorded in analysis_plan). Register
TissueResolve_state_aware as a distinct benchmark method, add a graceful
state-aware report subsection, and document feature status. Additive; default
path and all prior tests unchanged. Reference adaptation and expression
reconstruction remain deferred.
```
