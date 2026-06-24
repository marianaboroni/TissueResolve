# Project structure audit

Classification of every top-level area, a per-script review of the breast-cancer
harness, and a proposed simplified workflow. Goal: reduce folder/script
confusion without deleting anything users rely on.

## A. Core package — keep

`src/tissueresolve/` — the importable library and CLI.

- `api.py`, `cli.py`, `config.py`, `presets.py`, `results.py`, `utils.py`
- `io/` (bulk, spatial, reference, validation, autodetect, normalization-aware checks)
- `reference/` (build, markers, gene_filters, separability, resolution, hierarchy,
  hierarchical_build, pairwise_markers)
- `protocol/` (metadata, risk, mismatch, detect)
- `bulk/` (solver, pipeline, qc, hierarchical)
- `spatial/` (model, graph, pipeline, qc, neighbourhood, benchmark, hierarchical)
- `uncertainty/`, `plotting/`, `report/` (incl. new `unified.py`), `tuning/`, `benchmark/`

All required for the API/CLI. No removals.

## B. Tests — keep

`tests/` (unit + integration, offline), `examples/real_breast_cancer/tests/`
(harness), `benchmarks/tests/` (benchmark framework). All collected via
`pyproject.toml` `testpaths`.

## C. Examples — keep (harness scripts reviewed below)

`examples/real_breast_cancer/` — real-data validation harness (validation
tooling only; never alters core algorithms).

## D. Benchmarks — keep

`benchmarks/` — optional, offline-first benchmark harness (internal baselines +
optional external tools + import path). See `docs/benchmarking.md`.

## E. Documentation — keep

`README.md`, `docs/` (tutorial, quickstart, input_formats, advanced_parameters,
reporting, output_interpretation, qc_thresholds, benchmarking,
normalization_and_protocols, batch_effects, library_type_references,
development_workflow, this audit).

## F. Generated / local — must NOT be committed (git-ignored)

- `examples/real_breast_cancer/data/` (downloaded `.h5ad`, derived TSVs)
- `examples/real_breast_cancer/outputs/`
- `benchmarks/outputs/`
- `data/`, `.venv/`, `demo_demo/`, `demo_wiz/` (local scratch — should stay
  local; verify they are git-ignored, see note below)

> Note: `demo_demo/` and `demo_wiz/` look like local wizard/demo scratch
> directories. They are not referenced by the package or tests and should be
> git-ignored or removed locally (not done automatically here — no deletion
> without explicit instruction).

## Per-script review — `examples/real_breast_cancer/scripts/`

| Script | Purpose | Still needed? | Can merge? | Dev-only? | Notes |
|---|---|---|---|---|---|
| `00_download_data.py` | Download reference + Visium (opt-in) | Yes | No | No | Network only with `--run-real-data` |
| `01_prepare_reference.py` | Build TissueResolve reference | Yes | No | No | |
| `02_make_pseudobulk.py` | Make pseudobulk mixtures + truth | Yes | Could merge into bulk validation | No | Produces ground truth |
| `03_run_bulk_validation.py` | Flat bulk deconv vs truth | Yes | Could pair with 02 | No | |
| `04_run_spatial_validation.py` | Spatial deconv + QC | Yes | No | No | |
| `05_summarize_results.py` | Text summary of runs | Optional | Folds into 07 | Partly | Superseded by unified report |
| `06_resolution_spillover_analysis.py` | Separability/spillover diagnostics | Yes | Overlaps 08 | No | Consider merging with 08 |
| `07_generate_reports.py` | HTML reports + **unified `report.html`** | Yes | — | No | Now the main report entry point |
| `08_resolution_analysis.py` | Resolution recommender + family aggregation | Yes | Overlaps 06 | No | Merge candidate with 06 |
| `09_run_hierarchical_deconvolution.py` | Broad→fine bulk+spatial + outputs | Yes | No | No | |
| `_download_utils.py` | Download helpers | Yes | No | Yes (support) | Imported by 00 |
| `_harness.py` | Shared paths/IO/metrics | Yes | No | Yes (support) | Imported by all |

### Overlaps identified

- **06 vs 08**: both produce resolution/separability/spillover outputs. They can
  be merged into a single `resolution_analysis` step.
- **05** (text summary) is largely superseded by the unified `report.html` from
  07; keep as a quick CLI summary or fold into 07.

## Proposed simplified user-facing workflow

The current numbered scripts are kept (changing numbers risks breaking docs and
muscle memory), but the **recommended path** is now:

| Step | Script | Does |
|---|---|---|
| download | `00_download_data.py` | fetch reference + Visium (opt-in) |
| reference | `01_prepare_reference.py` | build the reference |
| bulk | `02_make_pseudobulk.py` → `03_run_bulk_validation.py` | pseudobulk + bulk validation |
| spatial | `04_run_spatial_validation.py` | spatial deconvolution + QC |
| hierarchical | `09_run_hierarchical_deconvolution.py` | broad→fine bulk + spatial |
| resolution | `08_resolution_analysis.py` (merge 06 in future) | separability/spillover/families |
| report | `07_generate_reports.py` | **builds the unified `outputs/report.html`** |

If a flatter scheme is later desired, the target names are:
`00_download_data`, `01_prepare_reference`, `02_run_bulk_validation`,
`03_run_spatial_validation`, `04_run_hierarchical_analysis`, `05_generate_report`
(02 would absorb pseudobulk generation; 04 the hierarchical run; 05 the unified
report). This is documented here rather than applied, to avoid breaking the
existing harness and its tests in this pass.

## Where are my results?

- **Main:** `examples/real_breast_cancer/outputs/report.html` (unified).
- Benchmarks: `benchmarks/outputs/benchmark_summary_report.html`.
- Detailed sub-reports under `outputs/{bulk,spatial,validation_summary}/` and
  `outputs/hierarchical/`. See each directory's `index.html` / `README.md`.
