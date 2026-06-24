# Report Path Consolidation Plan

Status: **IMPLEMENTED** (see "Implementation outcome" below). Goal: converge
TissueResolve's report generation onto a single canonical core path
(`report/unified.py` + a new core orchestration module), and retire the
duplicated `report/html.py` renderers safely and in stages.

## Implementation outcome (executed)

- **Stage A:** added `report/orchestration.py` — the single
  `generate_report(modality, source, out)` entry + `render_sections()` unified
  renderer. Unit-tested in `tests/report/test_orchestration.py`.
- **Stage C:** `report/__init__`, `api.py`, and the CLI report paths now route
  through `orchestration`. `tissueresolve run`, `tissueresolve report`,
  `bulk/spatial report`, and `api.generate_report` all emit the unified report.
- **Stage D:** no test migration needed — orchestration reuses the existing
  section *content* builders (`sections.py` for results-dir;
  `html.{bulk,spatial}_result_sections` for in-memory) and renders them through
  the unified shell, so all title/substring assertions still hold.
- **Stage E/F:** the duplicate page renderers were removed from `html.py`
  (`_generate_*_from_result`, `_generate_report_from_dir`, `_write`, `_section`,
  `_CSS`). The in-memory section builders + their HTML helpers were then
  **relocated** to a new `report/result_sections.py` (the in-memory counterpart
  of `sections.py`), so `html.py` is now a **pure compatibility shim**: three
  deprecated public functions that delegate to `orchestration`, plus
  backward-compat re-exports of the relocated builders. `sections.py`/`assets.py`/
  `interpretation.py` are retained as the **results-dir content layer** for the
  canonical path; `templates.py` provides HTML fragments for that content layer
  and still backs the harness's combined validation report. Their docstrings
  were updated to reflect this role (no longer "deprecated, do not extend").

`html.py` now contains no report-building logic of its own — only shims and
re-exports — so it can be deleted in a future cleanup once external callers stop
importing `tissueresolve.report.html` directly (in-repo callers already use the
package entry / `orchestration`).

The staged plan as originally written follows.

Date: 2026-06-02. Prerequisite met: alpha stabilization pass complete, full
offline suite green (1070 passed, 1 skipped).

---

## 0. TL;DR

There are **three** report renderers today, not two:

| # | Renderer | Entry | Section style | Consumers |
|---|---|---|---|---|
| A | `html.py` **in-memory** (`_generate_{bulk,spatial}_report_from_result`) | `report.generate_report(modality, result)` | "Input summary", "Deconvolution estimates (mRNA proportions)", "QC metrics", "Methods" | `tissueresolve run` (bundle), `api.generate_report` |
| B | `html.py` **results-dir** (`_generate_report_from_dir` → `sections/templates/assets`) | `report.generate_report(modality, results_dir)` | "Executive summary", "Key findings", "Single-cell reference quality", "Deconvolution predictions", "Prediction QC", "Detailed outputs" | CLI `report` / `bulk report` / `spatial report` |
| C | `unified.py` (`build_unified_report` + `Section`) | `07_generate_reports.py::generate_unified_report()` | "1. Executive decision summary", "2. Reference quality", … (QC-first decision workflow) | the real-data harness only |

**Canonical target = C** (the QC-first unified report). The consolidation moves
the *general* orchestration that currently lives in the 2095-line harness script
into the core package, repoints A and B at it, migrates the dependent tests, and
only then deprecates/removes `html.py` + `sections/templates/assets`.

Out of scope: **benchmark reports** use a self-contained `benchmarks/shared/report.py`
(no dependency on the core report package) — they are *not* part of this work.

---

## 1. Report entry points (mapped)

### 1.1 CLI (`src/tissueresolve/cli.py`)
- `tissueresolve report --modality bulk|spatial --results-dir` → `report.generate_report` (cli.py:987) → **renderer B**.
- `tissueresolve bulk report` (cli.py:699) and `tissueresolve spatial report` (cli.py:867) → `report.generate_report` → **renderer B**.
- `tissueresolve run` → `_write_run_report_bundle` (cli.py:584) → `report.generate_report(modality, result)` → **renderer A**; also writes `methods.txt` (cli.py:564 `methods_text`) and `warnings.json` (cli.py `_collect_run_warnings`).

### 1.2 API (`src/tissueresolve/api.py`)
- `api.generate_report(result, ...)` (api.py:388) → `html.generate_bulk_report` / `html.generate_spatial_report` → **renderer A**.

### 1.3 Core report package (`src/tissueresolve/report/`)
- `report/__init__.py` re-exports `generate_report`, `generate_bulk_report`, `generate_spatial_report` **from `html.py`** — the public surface.
- `report/html.py` dispatches A vs B on whether `source` is a path or an in-memory result; B pulls in `sections.py` (→ `interpretation.py`, `templates.py`) and `assets.py`.
- `report/unified.py` (renderer C primitives) is imported only by the harness and by two unified-based tests.

### 1.4 Real-data harness (`examples/real_breast_cancer/scripts/07_generate_reports.py`, ~2095 lines)
- `generate_bulk_outputs` / `generate_spatial_outputs` (07:249, 07:303) → `generate_report(modality, rdir)` → **renderer B** (per-modality sub-reports).
- `generate_unified_report()` (07:1044) → builds `Section`s using `components`, `glossary`, `figures` and calls `build_unified_report` → **renderer C** (the canonical `report.html` + `technical_appendix.html`).
- Supporting orchestration (general vs harness-specific):
  - **General, portable:** `_figure_cards`, `_caption_for`/`_CAPTIONS`, `_df_collapsible`, `_collect_warnings`, `_component_status`, `_decision_statuses`, `_hierarchy_summary`, `_state_aware_subsection`, `_gene_overlap_by_modality`, figure-routing sets (`_SUMMARY_FIGS`, `_BULK_APPENDIX_FIGS`, `_SPATIAL_APPENDIX_FIGS`, `_SIGNATURE_APPENDIX_FIGS`).
  - **Harness-specific (stays in 07):** `_load_reference`, `generate_diagnostic_figures`, `generate_benchmark_figures`, `_he_overlays`, hard-coded `examples/.../outputs` paths, dataset-specific captions, `_write_output_index`.

### 1.5 Benchmarks — OUT OF SCOPE
- `benchmarks/shared/report.py` (`build_modality_report`, `build_summary_report`) emits its own HTML; `run_all.py` uses it. No core-report dependency.

---

## 2. Files that depend on `report/html.py` (must be handled before removal)

### 2.1 Source
- `src/tissueresolve/report/__init__.py` — re-exports the three functions from `html.py`.
- `src/tissueresolve/api.py:388` — imports `report.html`.
- `src/tissueresolve/cli.py` — 4 call sites via `report.generate_report` (+ `methods_text`).
- `src/tissueresolve/report/html.py` itself imports `sections`, `templates`, `assets`, `methods_text`.
- `src/tissueresolve/report/sections.py` imports `assets`, `interpretation`, `templates`, `methods_text`.

### 2.2 Tests bound to html.py headings/behaviour (renderers A & B)
- `tests/report/test_report.py` — in-memory (A): asserts "Warnings", "Methods", "mRNA proportions"/"cell fractions", "did NOT converge", "Smoothing parameters", "Convergence", "poorly-separable"/"CRITICAL", "λ_spatial".
- `tests/report/test_results_dir_report.py` — A+B: heading list ("Single-cell reference quality", "Input data summary", "Deconvolution predictions", "Prediction QC", "Warnings and limitations", "Methods"), "Analysis plan", "Recommended merge families", "Family-level estimates", figure iframes.
- `tests/report/test_report_redesign.py` — B: "Executive summary", "Key findings", "Main publication figure", "Detailed outputs", "Main results interpretation", "mRNA-derived", "Bootstrap uncertainty was not computed", "No warnings." absence.
- `tests/test_compliance.py` — A: "Warnings", "Methods", "mRNA"/"cell fractions", "did NOT converge".
- `tests/test_report_stability.py` *(added this pass)* — A+B: "Single-cell reference quality", "Deconvolution predictions", "mRNA-derived", "Output files"/"Detailed outputs", `generate_spatial_report` failed-convergence.
- `tests/test_end_to_end_workflows.py` *(added this pass)* — A: report from in-memory run result ("mRNA" present); also `generate_report(bulk, results_dir)` ("Methods", "low R").
- `tests/test_import.py` — imports `tissueresolve.report` (smoke).

### 2.3 Tests already on the unified path (renderer C) — low/no migration
- `tests/report/test_report_redesign_ordering.py`, `tests/report/test_report_components.py` — use `unified.Section`/`build_unified_report` + `components`/`glossary`/`figures` directly.
- `examples/.../tests/test_report_figures.py`, `test_benchmark_split.py`, `test_report_simplification.py`, `test_state_aware_report.py` — use `figures.FigureManifest`/`components` (unified side).

---

## 3. Canonical future architecture

```
report/
  unified.py          # KEEP — Section, build_unified_report (renders report.html + appendix)
  components.py       # KEEP — figure_card, metric_grid, collapsible_table, nav, sections
  style.py            # KEEP — REPORT_CSS
  figures.py          # KEEP — FigureManifest / FigureRecord / status
  glossary.py         # KEEP — term definitions incl. no-ground-truth caveat
  methods_text.py     # KEEP — compose_bulk_methods / compose_spatial_methods
  interpretation.py   # KEEP (or fold relevant bits into orchestration) — data-driven prose
  warnings.py         # NEW  — single warning collector (see §3.2)
  orchestration.py    # NEW  — build_report_sections() + generate_report() (renderer C from any source)
  __init__.py         # re-export generate_report from orchestration (was html)
  html.py             # DEPRECATE → thin compatibility wrapper → eventually remove
  sections.py         # DEPRECATE → remove (folded into orchestration)
  templates.py        # DEPRECATE → remove (replaced by components/style/unified)
  assets.py           # DEPRECATE → remove (figure embedding via figures.py/components)
```

### 3.1 `report/orchestration.py` (the heart of the migration)
- `build_report_sections(modality, source, *, run_metadata=None, warnings=None, figures=None) -> list[Section]` — accepts a **results-dir path OR an in-memory pipeline result** and returns the ordered QC-first sections (+ a parallel appendix section list). Ports the *general* logic from `07.generate_unified_report` and its portable helpers (§1.4).
- `generate_report(modality, source, out=None) -> Path` — builds sections and calls `build_unified_report` to write `report.html` (+ `technical_appendix.html`). This becomes the single public entry, replacing both html.py renderers.
- A small `_load_from_results_dir()` and `_load_from_result()` adapter normalise the two source shapes into one section-builder input — eliminating the A/B split.

### 3.2 `report/warnings.py`
Unify the three warning collectors that exist today:
- `html.py::_collect_bulk_warnings` / `_collect_spatial_warnings`
- `cli.py::_collect_run_warnings` (added this pass)
- `07.py::_collect_warnings`
into one `collect_warnings(modality, *, result=None, results_dir=None) -> list[Warning]` so `warnings.json` and the report agree by construction.

---

## 4. Staged migration plan (each stage keeps the suite green)

### Stage A — Add core orchestration (no behaviour change)
- Create `report/orchestration.py` and `report/warnings.py` by **porting** the portable parts of `07.generate_unified_report` + helpers. Do **not** touch existing entry points yet.
- `report/__init__.py` keeps exporting the html-based `generate_report` (unchanged).
- New unit tests for `orchestration.build_report_sections` + `generate_report` against tiny synthetic results dirs and in-memory results (assert QC-first ordering, populated predictions, warnings surfaced, appendix written).
- **Exit criterion:** new module fully tested; nothing else changed; suite green.

### Stage B — Make the example script call core orchestration
- Repoint `07.generate_unified_report()` to call `report.orchestration.generate_report(...)`, leaving its harness-specific figure generation in place (it passes figures/manifest in).
- Keep `07.generate_bulk_outputs/generate_spatial_outputs` (renderer B sub-reports) untouched for now.
- **Exit criterion:** `examples/.../tests` green; the harness still produces `report.html` + `technical_appendix.html` byte-comparable in structure (heading set unchanged).

### Stage C — Make CLI/API use core orchestration
- Point `report.generate_report` (in `__init__.py`) at `orchestration.generate_report`.
- Update `api.generate_report` and `cli._write_run_report_bundle` to the unified entry (in-memory source).
- Update the CLI `report`/`bulk report`/`spatial report` commands (results-dir source).
- **Exit criterion:** CLI/API produce the unified QC-first report from both source shapes; `tissueresolve run` bundle still emits `report.html`/`methods.txt`/`warnings.json`.

### Stage D — Migrate tests
- Update the renderer-A/B tests in §2.2 to assert the **unified** section titles/strings (e.g. "Reference quality", "Bulk deconvolution") instead of html.py headings. Preserve the *intent* of every assertion (QC-first ordering, estimate-type disclaimer, warnings surfaced, failed-convergence shown, no empty figure cards, source-data links).
- Update the two tests added this pass (`test_report_stability.py`, `test_end_to_end_workflows.py`) accordingly.
- **Exit criterion:** full suite green with no `html.py` heading assertions remaining.

### Stage E — Deprecate `sections.py`/`templates.py`/`assets.py`
- Remove their last importers (only `html.py` after Stage C/D). Delete the three modules (or reduce to no-ops) once nothing imports them.
- **Exit criterion:** `grep` shows zero importers; suite green.

### Stage F — Remove `html.py` (or make it a compatibility wrapper)
- Option F1 (preferred): replace `html.py` with a thin shim whose `generate_report`/`generate_bulk_report`/`generate_spatial_report` delegate to `orchestration.generate_report`, emitting a `DeprecationWarning`. Keeps any external importers working.
- Option F2: delete `html.py` entirely once no in-repo importer remains and the public re-exports come from `orchestration`.
- **Exit criterion:** one canonical renderer; `html.py` is either a shim or gone; suite green; `docs/reporting.md` updated to drop the "two paths" language.

---

## 5. Migration risk

| Area | Risk | Mitigation |
|---|---|---|
| Porting ~500 lines of `generate_unified_report` + helpers into core | **High** — easy to subtly change section content/ordering | Stage A in isolation with dedicated ordering/content tests before any rewiring |
| ~7 test files assert html.py headings | **High** — heading strings change under unified | Stage D migrates intent, not literal strings; do it only after C |
| `run` bundle + `api.generate_report` now emit unified instead of A | Medium — output shape changes for users | Documented as the consolidation goal; covered by e2e tests |
| Harness-specific logic entangled with portable logic in 07 | Medium | Explicit portable/harness split in §1.4; only portable parts move |
| Hidden importers of `sections/templates/assets` | Low | `grep` gate at Stage E (only `html.py` imports them today) |
| Benchmark reports | None | Separate module, explicitly out of scope |

Overall: **high-risk refactor**, justified only because it removes a genuine
duplication. The staging keeps the suite green at every step and never deletes a
module while it still has importers.

---

## 6. Safest first implementation step

**Stage A, sub-step 1:** create `report/orchestration.py` with
`build_report_sections(modality, source, …)` that handles the **results-dir
source only** (the simpler of the two shapes), porting the portable section
builders from `07.generate_unified_report`, plus `report/warnings.py`. Add
unit tests asserting QC-first ordering, populated predictions, surfaced
warnings, and a written appendix — **without touching `__init__.py`, `cli.py`,
`api.py`, `html.py`, or `07`.** This introduces the canonical path with zero
blast radius and becomes the foundation everything else repoints to.

---

## 7. Explicitly NOT in this plan
- No code changes, deletions, feature additions, or commits (this is a plan).
- No change to benchmark report generation.
- No new report figures or sections beyond what renderer C already produces.
</content>
