# TissueResolve v0.1 Product-Scope Audit

This document classifies every implemented feature against the **core product
claim** and recommends a defensible v0.1 scope. It is a planning document only:
no code is changed, no file is deleted, and no test is removed.

## The core claim (what v0.1 must defend)

From `README.md` and `docs/output_interpretation.md`, the honest core claim is:

> TissueResolve estimates **RNA-derived composition** from **bulk RNA-seq** and
> **10x Visium** using a **shared single-cell reference**, and **reports which
> level of cellular resolution can be trusted** (hierarchy, separability,
> spillover, reference suitability, unresolved mass), with **publication-quality
> reports whose every figure has saved source data**.

A feature is **CORE** only if it satisfies ALL of:

1. supports the core claim,
2. is integrated into a default or documented public pathway (`api.py` / `cli.py`),
3. is tested in the default offline suite,
4. is surfaced in a report or saved output,
5. is documented, and
6. does not overclaim (estimate-type honesty preserved).

Anything that supports an *extended* claim but is gated, unvalidated, or not
surfaced in the standard report is **EXPERIMENTAL**. Anything that is a stub,
fabricated, unused, duplicated, undocumented-but-claimed, or untested is a
**QUARANTINE / REMOVE** candidate.

---

## (a) v0.1 CORE set

These features meet all six criteria. Paths are relative to repo root.

### Reference layer
- `src/tissueresolve/reference/build.py` — `ReferenceBuilder`, the single
  reference entry point. Tested, documented (`docs/tutorial.md`), used by
  `api.py`/`cli.py`. Implements "shared single-cell reference".
- `src/tissueresolve/reference/markers.py`, `gene_filters.py`,
  `pairwise_markers.py`, `within_family_markers.py` — protocol-aware composite
  marker selection. Tested; reached through the reference scope and hierarchy.
- `src/tissueresolve/reference/separability.py` — Bhattacharyya separability +
  `SeparabilityWarning` + `merge_nonseparable_types`. Surfaced in report
  (`plotting/separability_plots.py`), drives gating. Directly supports the
  "which resolution can be trusted" claim.
- `src/tissueresolve/reference/resolution.py` — resolvability classes and family
  inference; consumes separability metrics. Surfaced in report.
- `src/tissueresolve/reference/suitability.py` — reference suitability score
  (PASS/CAUTION/WARNING/FAIL). Surfaced in report; supports the trust claim.
- `src/tissueresolve/reference/hierarchy.py` — broad→fine arithmetic and
  unresolved-mass assembly. Default pathway. Core to hierarchy + unresolved-mass
  claims.
- `src/tissueresolve/reference/hierarchical_build.py` — hierarchical reference
  persistence. Supports the hierarchical claim; used in example stage 9.

### Bulk
- `src/tissueresolve/bulk/pipeline.py`, `qc.py`, `solver.py`,
  `bulk/__init__.py` — default bulk path (wNNLS, QC, optional explicit mRNA
  correction). Tested, surfaced, documented. `estimate_type='mRNA_proportion'`
  enforced (no overclaim).
- `src/tissueresolve/bulk/hierarchical.py` — two-level bulk hierarchical
  (default when broad/fine labels exist). `FEATURE_STATUS.md` line 13 = stable.

### Spatial
- `src/tissueresolve/spatial/pipeline.py`, `model.py`, `graph.py`, `qc.py`,
  `neighbourhood.py`, `auto_params.py`, `spatial/__init__.py` — default NB-CAR
  spatial path with recorded smoothing parameter (`auto_params`), spot QC and
  Moran's I. `estimate_type='spot_rna_composition'` (no overclaim).
- `src/tissueresolve/spatial/hierarchical.py` — two-level spatial hierarchical.
  Stable; **lacks a dedicated unit test file** (tested indirectly) → refactor item.

### Solver
- `src/tissueresolve/solver/*` (base, nnls, weighted_nnls, marker_nnls,
  ridge_nnls, ensemble, auto, `__init__`) — solver backbones and the default
  `auto` (gene-masking CV). Tested. `auto` is stable (`FEATURE_STATUS.md` line 12).
- `src/tissueresolve/validation/gene_masking.py` — the data-driven model-selection
  engine behind `auto`/`ensemble`. Tested. (Mislocated; see refactor notes.)

### Diagnostics shared with core
- `src/tissueresolve/benchmark/spillover.py` — spillover matrix and per-type
  risk; surfaced in report. Supports spillover claim.
- `src/tissueresolve/protocol/detect.py`, `mismatch.py` — protocol-aware checks.
- `src/tissueresolve/uncertainty/bootstrap.py` — bootstrap CIs (no doc page → refactor).

### I/O, config, results, plotting, report (core)
- `io/*` (facade, validation, autodetect, reference, bulk, spatial),
  `config.py`, `results.py`, `presets.py`, top-level `api.py`, `cli.py`,
  `__init__.py`.
- Plotting core: `style.py`, `export.py`, `palette.py`, `captions.py`,
  `bulk_plots.py`, `spatial_plots.py`, `separability_plots.py`,
  `spillover_plots.py`, `resolution_plots.py`, `qc_plots.py`,
  `reference_qc_plots.py`, `signature_qc_plots.py`, `histology.py`,
  `summary_figures.py`. All enforce "every figure has saved data".
- Report core: `report/__init__.py`, `report/components.py`, `report/style.py`,
  `report/glossary.py`, `report/interpretation.py`, `report/methods_text.py`.

### Benchmark harness (core for the optional benchmark claim)
- `benchmarks/run_all.py`, the `benchmarks/shared/*` utilities, internal
  baselines, registry, and offline tests. These support the explicitly-optional
  benchmark claim and are not on the default deconvolution path, but they are
  tested and honest, so they remain core to the benchmark sub-product.

---

## (b) EXPERIMENTAL set (behind flags, labelled, NOT part of the claim)

| Module | Why experimental | Gate / label |
|---|---|---|
| `src/tissueresolve/bulk/state_aware_hierarchical.py` | 3-level broad→cell-type→state; not validated on real data; report path not wired (see audit table) | `--state-aware`; `feature_status='experimental'`; `FEATURE_STATUS.md` line 15 |
| `src/tissueresolve/reference/three_level_hierarchy.py` | only used on state-aware path | imported by `api.py` only |
| `src/tissueresolve/reference/granular_signatures.py` | only used on state-aware path | imported by `api.py` only; `FEATURE_STATUS.md` line 14 |
| `benchmarks/synthetic/state_hierarchy.py`, `run_state_aware_synthetic.py` | offline validation of the experimental state feature | not in `run_all.py` |
| `benchmarks/tests/test_state_aware_*` , `tests/.../test_state_aware_*` | tests of the experimental feature | keep, never weaken |
| `src/tissueresolve/spatial/benchmark.py` | self-consistency only; marker-recovery / synthetic-truth / dashboard not wired | `FEATURE_STATUS.md` line 18 = partial |
| `src/tissueresolve/report/unified.py`, `report/figures.py` | report-redesign work; not the single canonical report path yet | keep; consolidate in v0.2 |
| `src/tissueresolve/plotting/benchmark_plots.py` | explicitly "synthetic benchmarks only, not real-data validation" | benchmark scope |

**Honesty requirement:** the state-aware path must keep
`feature_status='experimental'` in metadata and must NOT be presented as a
validated default in README/CLI help/report.

---

## (c) QUARANTINE / REMOVE candidates

Conservative stance: prefer **quarantine** over removal; never delete tests.

| Module | Issue | Recommendation |
|---|---|---|
| `src/tissueresolve/tuning/__init__.py` | Stub that previously emitted **fabricated metrics** (scientific-integrity violation); now raises `NotImplementedError`. Already quarantined, tested (`tests/test_tuning_quarantine.py`). | **Keep quarantined.** Do not resurrect tuning for v0.1; point users to `solver='auto'` and `spatial.auto_params`. Remove only after v0.1 ships and tests confirm no importers. |
| `src/tissueresolve/validation/__init__.py` | Empty (60 bytes) stub; misleading because real validation lives in `io/validation.py` and the only real module here is `gene_masking.py`. | **Quarantine the directory shape**: either document the split or relocate `gene_masking.py` (see refactor). Do not delete `gene_masking.py`. |
| `src/tissueresolve/report/templates.py` | Bare HTML f-string fragments, **untested**, imported only by `sections.py`; likely superseded by `components.py` + `unified.py`. | **Quarantine after the report path is consolidated.** Tie its fate to `sections.py`. |
| `src/tissueresolve/report/sections.py` + `report/assets.py` | `sections.py` overlaps `unified.py`; `assets.py` is untested. Two parallel report-building paths. | **Refactor/consolidate** to one report path; quarantine the loser. Keep until consolidation lands. |
| `src/tissueresolve/bulk/state_aware_hierarchical.py` (report side) | Deconvolution works but the **standard report harness cannot consume `StateAwareBulkResult`** (`report/html.py` expects `result.deconv`); report is custom-coded in `examples/.../07_generate_reports.py`. | Keep the module **experimental**; quarantine the *claim* that state-aware reports work in the standard harness until the gap is fixed. |
| `examples/real_breast_cancer/scripts/07_generate_reports.py` | Custom report code duplicating `tissueresolve.report.html` (DRY violation; a workaround for the gap above). | **Simplify** after state-aware report is integrated into `src`. Do not delete (it is validation tooling). |
| `benchmarks/shared/protocol_detection.py` | Not imported by current runners; untested; "available for future integration". | **Quarantine candidate** if still unused at v0.1 freeze. Wire it in or move out of the default benchmark import graph. |

**No fabricated-metric code remains live** other than the already-neutralised
`tuning` stub. No module is a dead unused orphan except the empty
`validation/__init__.py` shell. There is one true duplicate-logic pair worth
acting on: the report builders (`sections.py`/`templates.py`/`assets.py` vs
`unified.py`/`components.py`) and the example report workaround.

---

## Documentation-claim audit

| Claim | Source | Classification | Honest rewrite / action |
|---|---|---|---|
| RNA-derived composition from bulk + Visium via shared sc reference | README 3-22 | **supported** | Keep. `estimate_type` enforced in `results.py`. |
| Hierarchical broad→fine is default when labels exist | README 12-39; `cli.py` | **supported** | Keep. `_select_resolution_mode` falls back to flat with warning. |
| Reports which resolution can be trusted (hierarchy/separability/spillover/suitability/unresolved mass) | abstract; README 28-72 | **supported** | Keep. All five integrated into report. |
| State-aware broad→cell-type→state | README 6; `api.py`; `cli.py --state-aware` | **partially supported** | Rewrite: "Experimental: state-aware deconvolution runs behind `--state-aware`; **standard HTML report does not yet render state-aware results** (custom example script only). Not validated on real data." Do not present as a finished default. |
| Separability (Bhattacharyya) diagnostics | README 69 | **supported** | Keep. |
| Spillover matrix + per-type risk | README 70 | **supported** | Keep. |
| Publication HTML reports with source data per figure | README 29,73-74 | **supported** | Keep. `plotting/export.py` writes `.data.tsv`. |
| Real-data breast-cancer validation harness | README 44-45 | **supported** | Keep. Tooling only; does not weaken core. |
| Optional benchmark vs internal baselines + external tools | README 43; `docs/benchmarking.md` | **supported** | Keep. Caveat: only executed tools ranked; exported-only not counted. |
| Cell-type-specific **expression reconstruction** | README differentiator (inferred) | **partially supported / overclaim risk** | Rewrite to "cell-type-level deconvolution **estimates**" (proportions/composition), NOT expression reconstruction. `FEATURE_STATUS.md` line 21 marks reconstruction **not implemented**; must not be claimed. |
| Superiority over other tools | README 55-76 differentiator table | **unsupported (correctly caveated)** | Keep current framing: "high-level orientation, not a claim that TissueResolve outperforms." |
| Spatial accuracy vs other spatial tools | README differentiator (implied) | **unsupported** | Keep caveat: "No absolute accuracy claimed without ground truth" (`docs/benchmarking.md`). |
| Graph-aware NB-CAR spatial model | README 26 | **supported** | Keep. Smoothing param recorded. |
| Normalization/protocol/library/batch-aware diagnostics | README 40 | **supported** | Keep. |
| Family-aware handling of non-separable types | README 37 | **supported** | Keep. |
| Reference adaptation (residual-driven gene-weight update) | not in README; `FEATURE_STATUS.md` line 20 | **not implemented** | Must not be claimed anywhere. Keep absent. |

---

## Is the package too broad for v0.1?

**Yes, mildly.** The deconvolution core (reference + bulk + spatial + solver +
diagnostics + report + io + config/results) is coherent, tested, and matches the
claim. The breadth risk comes from: (1) the experimental **state-aware** stack
that is not report-integrated, (2) **two parallel report-builder code paths**,
(3) a large multi-tool **external benchmark** surface with many untested wrapper
modules, and (4) structural confusion (`validation/` shell, `tuning/` stub,
`benchmark/` vs `benchmarks/` name collision, a core module importing from
`benchmarks/shared/batch_effects.py`).

The correct v0.1 move is **scope clarity, not deletion**: clearly label
experimental/partial, consolidate the report path, and fix the
core→benchmarks dependency and naming — before adding any new feature.
