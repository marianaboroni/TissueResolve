# TissueResolve Architecture & Dependency Map

Derived from the recorded import edges and module inventory. Indented trees show
"imports / depends on". This is a read-only map; no code is changed.

## (a) Core workflow dependency chain

### Entry points
```
__init__.py            (re-exports the 5 public functions)
  └─ api.py
cli.py
  └─ api.py
api.py / cli.py
  ├─ config.py
  ├─ results.py
  └─ io/__init__.py
       ├─ io/validation.py
       ├─ io/autodetect.py   (cli)
       ├─ io/reference.py
       ├─ io/bulk.py
       └─ io/spatial.py
```

### Reference construction
```
api.build_reference / cli
  └─ reference/build.py
       └─ io/reference.py, io/validation.py
reference/__init__.py  (aggregator)
  └─ build, gene_filters, markers, pairwise_markers,
     separability, resolution, hierarchy
markers.py
  └─ gene_filters.py
pairwise_markers.py
  └─ (used by) within_family_markers.py
hierarchy.py
  ├─ separability.py
  ├─ within_family_markers.py
  └─ resolution.py
suitability.py
  ├─ separability.py
  └─ benchmarks/shared/batch_effects.py   (CROSS-SCOPE, see section c)
```

### Bulk deconvolution (default path)
```
api.deconv_bulk
  ├─ bulk/pipeline.py
  │    ├─ bulk/solver.py
  │    └─ bulk/qc.py
  ├─ bulk/hierarchical.py        (when broad/fine labels exist; default)
  │    ├─ bulk/pipeline.py
  │    ├─ reference/hierarchy.py
  │    ├─ reference/separability.py
  │    └─ reference/within_family_markers.py
  └─ (optional) solver/auto.py | solver/ensemble.py
       └─ validation/gene_masking.py
uncertainty/bootstrap.py         (CI when enabled)
```

### Spatial deconvolution (default path)
```
api.deconv_spatial
  ├─ spatial/pipeline.py
  │    ├─ spatial/model.py
  │    │    └─ protocol/mismatch.py
  │    ├─ spatial/graph.py
  │    ├─ spatial/neighbourhood.py
  │    ├─ spatial/qc.py
  │    └─ spatial/auto_params.py     (when lambda_spatial not given; default)
  └─ spatial/hierarchical.py         (when mapping available)
       ├─ spatial/pipeline.py
       └─ reference/hierarchy.py
```

### Report (default path)
```
api.generate_report / cli
  └─ report/__init__.py
       └─ report/html.py
            ├─ report/methods_text.py
            └─ report/sections.py
                 ├─ report/assets.py
                 ├─ report/interpretation.py
                 ├─ report/templates.py
                 └─ report/methods_text.py
report/components.py
  └─ report/style.py
```

### Plotting (default path, via api.plot_results)
```
plotting/bulk_plots.py, spatial_plots.py
  ├─ plotting/captions.py
  ├─ plotting/style.py
  └─ plotting/export.py
plotting/{reference_qc,signature_qc,resolution,separability,spillover,
          summary_figures,histology}_plots.py
  ├─ plotting/export.py
  ├─ plotting/palette.py   (most)
  └─ plotting/style.py
```

**Core foundation modules (most-depended-on):** `plotting/export.py`,
`plotting/style.py`, `plotting/palette.py` (each imported by 13–15 plotting
modules); `results.py`, `config.py`, `reference/separability.py`,
`reference/hierarchy.py`, `benchmarks/shared/base.py`.

## (b) Optional / experimental dependencies

### State-aware (experimental, gated by `--state-aware`)
```
api.deconv_bulk (state_aware=True)
  └─ bulk/state_aware_hierarchical.py
       ├─ reference/hierarchy.py
       └─ reference/three_level_hierarchy.py
api (state-aware reference panels)
  └─ reference/granular_signatures.py
       └─ reference/within_family_markers.py
(report) examples/real_breast_cancer/scripts/07_generate_reports.py
       └─ custom state-aware subsection  (NOT report/html.py)
```

### Benchmark harness (optional sub-product)
```
benchmarks/run_all.py
  ├─ benchmarks/bulk/run_bulk_benchmark.py
  │    └─ shared/{metrics,normalization,method_registry,method_selection,
  │              io,report,plotting,environment,analysis,synthetic,imported}
  ├─ benchmarks/spatial/run_spatial_benchmark.py
  │    └─ shared/{...,spatial_multimetric_ranking}
  └─ shared/{report,io}
benchmarks/run_real_external_benchmark.py
  └─ shared/{io,environment,composite_score,prepare_external_inputs}
method_registry.py
  └─ shared/base.py  (+ all method wrappers -> shared/base, environment, _export)
```

### Experimental benchmark (state-aware synthetic)
```
benchmarks/synthetic/run_state_aware_synthetic.py
  └─ benchmarks/synthetic/state_hierarchy.py
```

## (c) Unclear or circular dependencies

1. **Core → benchmarks back-edge (worst architectural smell).**
   `src/tissueresolve/reference/suitability.py` imports
   `benchmarks/shared/batch_effects.py`. A shipped library module depends on the
   (notionally optional, dev-only) benchmark tree. Not circular, but it inverts
   the intended layering and means the benchmark code is effectively a runtime
   dependency of the core. → relocate `batch_effects.py` into
   `src/tissueresolve/` (e.g. `reference/` or a `diagnostics/` module) in v0.2.

2. **Two parallel report builders.** `report/html.py` → `sections.py` →
   {`assets`, `interpretation`, `templates`, `methods_text`} coexists with
   `unified.py` → `components.py` → `style.py`. Both build HTML; the boundary is
   unclear and `templates.py`/`assets.py` are untested. Not circular but
   ambiguous ownership. → consolidate to one path.

3. **Name collision `benchmark/` vs `benchmarks/`.** `src/tissueresolve/benchmark/`
   (in-package, holds `spillover.py`) versus top-level `benchmarks/` (dev harness).
   Confusing for imports and discovery. → rename in-package to `diagnostics/` or
   merge `spillover.py` near separability.

4. **Self-edges in recorded graph** (`normalization→normalization`,
   `synthetic→synthetic`, `spatial_multimetric_ranking→self`,
   `composite_score→benchmarks/composite_score.py`) are parse artifacts of the
   edge extractor, not real cycles. The `composite_score` edge points at a
   non-`shared` path that does not exist in the inventory → verify the import
   target is `benchmarks/shared/composite_score.py`.

5. **`validation/` shape vs reality.** `validation/__init__.py` is empty and the
   directory name implies it owns validation, but input validation is in
   `io/validation.py` and only `gene_masking.py` lives under `validation/`. Layout
   misleads readers. → document or relocate.

No true import cycles were found in the core deconvolution graph.

## (d) Modules removable without breaking core

"Core" = the default bulk/spatial/reference/report/plot path reachable from
`api.py`/`cli.py`. The following can be removed/quarantined without breaking it
(tests for them stay):

- `src/tissueresolve/tuning/__init__.py` — quarantined stub; nothing on the
  default path imports it.
- The entire `benchmarks/` tree — optional dev harness. **Caveat:** core
  currently imports `benchmarks/shared/batch_effects.py` via
  `reference/suitability.py`; that one file must be relocated first, after which
  the rest of `benchmarks/` is fully detachable from the shipped library.
- State-aware stack (`bulk/state_aware_hierarchical.py`,
  `reference/three_level_hierarchy.py`, `reference/granular_signatures.py`) —
  only reached when `state_aware=True`. Disabling the flag removes them from the
  default path.
- `report/templates.py` and `report/assets.py` — only reached via `sections.py`;
  if the report path is consolidated onto `unified.py`/`components.py`, these
  drop out of core (quarantine, do not delete yet).
- Benchmark-only plotting (`plotting/benchmark_plots.py`,
  `benchmark_report_plots.py`, `bulk_benchmark_plots.py`,
  `spatial_benchmark_plots.py`) — not imported by `api.plot_results`; tied to the
  benchmark sub-product.
- `benchmarks/shared/protocol_detection.py` — not imported by any current runner.

**Not removable** (load-bearing for core): `results.py`, `config.py`, all of
`io/*` used by pipelines, the reference marker/separability/resolution/hierarchy
chain, `bulk/*` (pipeline/solver/qc/hierarchical), `spatial/*`
(pipeline/model/graph/qc/neighbourhood/auto_params/hierarchical), `solver/*`,
`validation/gene_masking.py`, `uncertainty/bootstrap.py`,
`benchmark/spillover.py`, `protocol/*`, and the report/plotting foundation
modules (`export`, `style`, `palette`, `captions`, `components`, `methods_text`,
`interpretation`, `glossary`, `html`).
