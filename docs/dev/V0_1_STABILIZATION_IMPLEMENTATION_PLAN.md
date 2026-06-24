# v0.1 stabilization implementation plan (Part 1)

Plan (no code yet) for the low/medium-risk stabilization actions from
`PRODUCT_SCOPE_AUDIT.md` / `V0_1_PRUNING_PLAN.md`.

## Preflight findings
- **Core→benchmarks back-edge = exactly one import:**
  `src/tissueresolve/reference/suitability.py:140`
  `from benchmarks.shared.batch_effects import compute_celltype_batch_confounding`.
  This is the only `benchmarks` import under `src/tissueresolve/`.
- **Expression reconstruction** is mentioned only as *deferred/not implemented*
  (README "deferred" list + audit docs). No overclaim to fix — just keep it that
  way and add a guard test.
- **State-aware** is already labelled experimental in `FEATURE_STATUS.md`,
  `V0_1_SCOPE.md`, the `--state-aware` CLI help, and README Status. Add the exact
  required wording + a report-text note + tests.
- Report already split into concise `report.html` + `technical_appendix.html`
  (last step); keep it.

## Low-risk changes to apply now
1. **Docs/label alignment (Part 2/7/8):** add the exact "Experimental: …" and
   "Cell-type-specific expression reconstruction is planned/deferred…" wordings
   to README + FEATURE_STATUS + V0_1_SCOPE; verify CLI `--state-aware` help wording;
   keep external-benchmark caveat. Files: `README.md`, `docs/FEATURE_STATUS.md`,
   `docs/V0_1_SCOPE.md`, `src/tissueresolve/cli.py` (help only).
2. **Architecture docs (Part 4/5):** append "Current naming caveats and planned
   cleanup" + "Canonical report path vs deprecated modules" to
   `docs/ARCHITECTURE_DEPENDENCY_MAP.md`. Deprecation docstrings already on
   `report/{sections,templates,assets}.py`.
3. **Stabilization report (Part 11):** `docs/V0_1_STABILIZATION_REPORT.md`.

## Medium-risk change to apply now
4. **Remove core→benchmarks dependency (Part 3):**
   - Create `src/tissueresolve/diagnostics/__init__.py` + move the batch-effect
     utility to `src/tissueresolve/diagnostics/batch_effects.py` (copy verbatim —
     no algorithm change).
   - Repoint `suitability.py` import to the new core location.
   - Make `benchmarks/shared/batch_effects.py` re-export from the core location
     (back-compat) so benchmark code keeps working.
   - Files: new `src/tissueresolve/diagnostics/`, `reference/suitability.py`,
     `benchmarks/shared/batch_effects.py`.

## Tests to add
- import-guard: no file under `src/tissueresolve/` imports from top-level
  `benchmarks` (AST/text scan).
- suitability still computes the batch/confounding component (regression).
- benchmark batch-effect utilities still import (via re-export).
- README/FEATURE_STATUS: state-aware experimental; expression reconstruction not
  implemented; CLI `--state-aware` help says experimental.
- tuning still quarantined (existing test).

## Defer (not in this step)
- Move `validation/gene_masking.py` → `solver/` and dissolve empty `validation/`.
- Rename `src/tissueresolve/benchmark/` → `diagnostics/` (would now collide with
  the new `diagnostics/batch_effects.py`; do the spillover-module rename later as
  its own refactor — for now `diagnostics/` holds only batch_effects).
- Delete `report/templates.py` / consolidate the legacy per-modality report path.
- 10-section main-report renumber.
- Validate state-aware on real 3-level data.

## Rollback risk
Low overall. The batch_effects move is the only behavioral-adjacent change; it is
a verbatim copy + re-export, covered by a regression test, and reversible by
restoring the import. Doc/label changes are zero-risk.
