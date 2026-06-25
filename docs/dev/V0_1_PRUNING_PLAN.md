# TissueResolve v0.1 Pruning Plan

A conservative, planning-only set of actions to tighten v0.1 scope. **No code is
changed and no test is removed by this document.** Severity is the priority of
acting; risk is the danger of the action itself. Prefer quarantine over removal.

## Immediate cleanup (low risk, label/doc only)

| Action | Severity | Files affected | Risk | Tests needed | Expected benefit |
|---|---|---|---|---|---|
| Mark state-aware as experimental everywhere it surfaces (README, CLI help, report); state plainly that standard HTML report does not render state-aware results yet | High | `README.md`, `src/tissueresolve/cli.py` (help text), `docs/FEATURE_STATUS.md` (already done) | Low | Assert `feature_status=='experimental'` in `tests/test_state_aware_integration.py` (exists) | Removes the only material overclaim |
| Reword "expression reconstruction" → "cell-type-level deconvolution estimates" in README differentiator | High | `README.md` | Low | None | Aligns claim with `FEATURE_STATUS.md` line 21 (not implemented) |
| Document the `validation/` vs `io/validation.py` split and the `benchmark/` vs `benchmarks/` naming in a short architecture note | Medium | `docs/ARCHITECTURE_DEPENDENCY_MAP.md` (this set) | Low | None | Removes reader confusion |
| Add missing doc pages for surfaced modules lacking docs (`uncertainty/bootstrap.py`, `plotting/export.py`, `report/components.py`) | Low | `docs/` | Low | None | Closes documented-but-thin gaps |

## Safe quarantine (keep file + tests; remove from default path / claim)

| Action | Severity | Files affected | Risk | Tests needed | Expected benefit |
|---|---|---|---|---|---|
| Keep `tuning/__init__.py` quarantined (already raises `NotImplementedError`); do not resurrect tuning for v0.1 | High | `src/tissueresolve/tuning/__init__.py` | Low | `tests/test_tuning_quarantine.py` (exists) | Prevents return of fabricated-metric code |
| Quarantine the "state-aware reports work in the standard harness" claim until the report gap is fixed | High | `README.md`, `src/tissueresolve/report/html.py` (doc), `examples/.../07_generate_reports.py` | Low | Existing example test stays | Honest reporting |
| Flag `benchmarks/shared/protocol_detection.py` as unused-at-v0.1 (wire in or quarantine from default import graph) | Medium | `benchmarks/shared/protocol_detection.py` | Low | Add a small unit test if retained | Shrinks untested surface |
| Quarantine `report/templates.py` only after report path consolidation (it is untested and superseded) | Medium | `src/tissueresolve/report/templates.py` | Medium | Cover replacement path in `tests/report/` | One report path |

## Needs refactor (behavior-preserving structural work, do with regression tests)

| Action | Severity | Files affected | Risk | Tests needed | Expected benefit |
|---|---|---|---|---|---|
| Relocate `benchmarks/shared/batch_effects.py` into `src/tissueresolve/` (e.g. `reference/batch_effects.py` or a `diagnostics/`); update `reference/suitability.py` import | High | `benchmarks/shared/batch_effects.py`, `src/tissueresolve/reference/suitability.py` | Medium | Regression test: `suitability` still scores batch component; benchmark tests still import | Core stops depending on dev benchmark tree |
| Consolidate the two report builders (`sections.py`+`templates.py`+`assets.py` vs `unified.py`+`components.py`) onto one canonical path | High | `report/html.py`, `report/sections.py`, `report/unified.py`, `report/templates.py`, `report/assets.py` | Medium | Golden-file/HTML-content tests in `tests/report/` before refactor | Removes duplicate logic |
| Integrate state-aware result into the standard report harness so `report/html.py` can render `StateAwareBulkResult` (e.g. adapter exposing `.deconv`-like view) | High | `src/tissueresolve/report/html.py`, `src/tissueresolve/bulk/state_aware_hierarchical.py` | Medium | New failing-first test: standard harness renders state-aware section; then port example test | Removes DRY violation in example script |
| Add a dedicated unit test for `spatial/hierarchical.py` (currently tested only indirectly) | Medium | `tests/spatial/` (new test) | Low | The new test itself | Direct coverage for a core module |
| Add unit tests for untested core-ish modules: `benchmarks/shared/{io,report,plotting,environment,imported,prepare_external_inputs,composite_score}` | Medium | `benchmarks/tests/` (new) | Low | The new tests | Coverage for benchmark sub-product |
| Relocate `validation/gene_masking.py` next to its importers (`solver/`) or document the split; remove the empty `validation/__init__.py` shell only after relocation | Medium | `src/tissueresolve/validation/`, `src/tissueresolve/solver/auto.py`, `solver/ensemble.py` | Medium | Regression: `auto`/`ensemble` still select via gene-masking CV | Clearer layout |
| Consider splitting `reference/hierarchy.py` (1044 lines, legacy + current APIs) into legacy vs current modules | Low | `src/tissueresolve/reference/hierarchy.py` | Medium | Existing `test_hierarchical.py` must stay green | Maintainability |
| Rename in-package `src/tissueresolve/benchmark/` (holds `spillover.py`) to avoid collision with top-level `benchmarks/` | Low | `src/tissueresolve/benchmark/`, importers of `benchmark.spillover` | Medium | Import smoke test | Removes naming confusion |

## Defer to v0.2 (do not build for v0.1)

| Action | Severity | Files affected | Risk | Tests needed | Expected benefit |
|---|---|---|---|---|---|
| Validate state-aware (broad→cell-type→state) on real data, then promote from experimental | High | state-aware stack | n/a | real-data validation under explicit flag | Honest promotion path |
| Wire spatial benchmark extras (marker-recovery, synthetic-truth, dashboard) | Medium | `src/tissueresolve/spatial/benchmark.py` | n/a | benchmark tests | Completes "partial" feature |
| Hyper-parameter tuning (replacement for quarantined `tuning`) | Low | new module | n/a | full tests + no fabricated metrics | Real tuning without integrity risk |
| Reference adaptation, cell-type expression reconstruction | Low | new modules | n/a | full tests | New deferred features (currently must not be claimed) |

## Remove after confirmation (only once v0.1 freezes and importers verified)

| Action | Severity | Files affected | Risk | Tests needed | Expected benefit |
|---|---|---|---|---|---|
| Remove `validation/__init__.py` empty shell — only after `gene_masking.py` relocated and no importers remain | Low | `src/tissueresolve/validation/__init__.py` | Medium | Import scan shows zero references | Cleaner tree |
| Remove `report/templates.py` — only after report consolidation and confirming `sections.py` no longer used | Low | `src/tissueresolve/report/templates.py` | Medium | Report tests green without it | Less dead code |
| Remove `tuning/` directory — only after v0.1 ships, with quarantine test migrated | Low | `src/tissueresolve/tuning/` | Medium | Migrate `test_tuning_quarantine.py` expectations | Removes neutralised stub |

> Removal rule (from project policy): never delete tests; always confirm zero
> importers with a scan before deleting any module; prefer quarantine.

---

## Proposed simplified v0.1 package structure

Target layout for `src/tissueresolve/`. "Belongs" / "does not belong" / moves are
proposals; nothing is moved by this document.

### `reference/`
- **Belongs:** `build.py`, `gene_filters.py`, `markers.py`, `pairwise_markers.py`,
  `within_family_markers.py`, `separability.py`, `resolution.py`, `suitability.py`,
  `hierarchy.py`, `hierarchical_build.py`, `__init__.py`.
- **Does not belong (v0.1):** `three_level_hierarchy.py`, `granular_signatures.py`
  → move under an `experimental/` namespace or keep but clearly gate as
  state-aware-only.
- **Moves/merges:** receive `batch_effects.py` from `benchmarks/shared/`; consider
  splitting `hierarchy.py` legacy vs current.

### `bulk/`
- **Belongs:** `pipeline.py`, `solver.py`, `qc.py`, `hierarchical.py`, `__init__.py`.
- **Does not belong (v0.1 default):** `state_aware_hierarchical.py` → `experimental/`
  or keep but flagged experimental and excluded from default claims.
- **Moves/merges:** none.

### `spatial/`
- **Belongs:** `pipeline.py`, `model.py`, `graph.py`, `qc.py`, `neighbourhood.py`,
  `auto_params.py`, `hierarchical.py`, `__init__.py`.
- **Does not belong:** `benchmark.py` (partial, self-consistency only) → move to a
  `benchmark`/`diagnostics` scope or keep but mark partial.
- **Moves/merges:** add dedicated test for `hierarchical.py`.

### `hierarchy/` (proposed new home for cross-modality hierarchy arithmetic)
- **Belongs:** the level-agnostic hierarchy assembly currently in
  `reference/hierarchy.py` (build/validate/aggregate/assemble), shared by bulk and
  spatial.
- **Does not belong:** marker selection (stays in `reference/`).
- **Moves/merges:** optional — only if `hierarchy.py` is split; otherwise leave in
  `reference/`.

### `solver/`
- **Belongs:** `base.py`, `nnls.py`, `weighted_nnls.py`, `marker_nnls.py`,
  `ridge_nnls.py`, `ensemble.py`, `auto.py`, `__init__.py`.
- **Does not belong:** model-selection engine should be co-located here.
- **Moves/merges:** move `validation/gene_masking.py` → `solver/gene_masking.py`.

### `diagnostics/` (proposed; replaces in-package `benchmark/` name)
- **Belongs:** `spillover.py` (from `benchmark/`), `batch_effects.py` (from
  `benchmarks/shared/`), and could host separability/resolution if a cross-cutting
  diagnostics home is preferred.
- **Does not belong:** the dev benchmark harness (stays at top-level `benchmarks/`).
- **Moves/merges:** rename `src/tissueresolve/benchmark/` → `diagnostics/`.

### `plotting/`
- **Belongs (core):** `style.py`, `export.py`, `palette.py`, `captions.py`,
  `bulk_plots.py`, `spatial_plots.py`, `qc_plots.py`, `separability_plots.py`,
  `spillover_plots.py`, `resolution_plots.py`, `reference_qc_plots.py`,
  `signature_qc_plots.py`, `histology.py`, `summary_figures.py`, `__init__.py`.
- **Does not belong (core):** `benchmark_plots.py`, `benchmark_report_plots.py`,
  `bulk_benchmark_plots.py`, `spatial_benchmark_plots.py` → group under a
  `plotting/benchmark/` subpackage tied to the benchmark sub-product.
- **Moves/merges:** group the four benchmark plot modules into one subfolder.

### `report/`
- **Belongs (one canonical path):** `__init__.py`, `html.py`, `components.py`,
  `style.py`, `glossary.py`, `interpretation.py`, `methods_text.py`, `figures.py`.
- **Does not belong / consolidate:** `sections.py`, `templates.py`, `assets.py`,
  `unified.py` — pick ONE builder path; quarantine the rest after migration.
- **Moves/merges:** fold the example state-aware report logic
  (`examples/.../07_generate_reports.py`) back into `html.py`.

### `benchmark/` (top-level `benchmarks/`)
- **Belongs:** `run_all.py`, `run_real_external_benchmark.py`,
  `import_external_results.py`, `shared/*`, `bulk/*`, `spatial/*`, `synthetic/*`,
  `tests/*`.
- **Does not belong:** `batch_effects.py` (move into `src` core) — it is the one
  module the shipped library imports.
- **Moves/merges:** none beyond `batch_effects.py`.

### `io/`
- **Belongs:** `__init__.py`, `validation.py`, `autodetect.py`, `reference.py`,
  `bulk.py`, `spatial.py`.
- **Does not belong:** nothing extra. This is the correct home for input
  validation (clarify vs the `validation/` package).
- **Moves/merges:** none.

### `validation/`
- **Belongs (after refactor):** nothing — the only real module
  (`gene_masking.py`) moves to `solver/`, and input validation already lives in
  `io/validation.py`.
- **Does not belong:** the empty `__init__.py` shell (remove after relocation).
- **Moves/merges:** dissolve the package; relocate `gene_masking.py` → `solver/`.
