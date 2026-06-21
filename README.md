# TissueResolve

**TissueResolve:** resolution-aware deconvolution for bulk RNA-seq and 10x Visium
spatial transcriptomics.

TissueResolve estimates RNA-derived cell-type and cell-state composition using a
shared single-cell (or single-nucleus) reference layer for both bulk and spatial
workflows. It combines protocol-aware gene selection, robust QC, separability
and spillover diagnostics, normalization/library/batch awareness, and
publication-ready HTML reports.

When the reference carries **broad** and **fine** cell-type annotations,
TissueResolve uses **hierarchical broad→fine** deconvolution by default
(`--resolution-mode auto`): it estimates broad cell-type families first, then
fine subpopulations within each family, and reports `unresolved_<family>` mass
where subtypes are not separable. Flat (fine-only) deconvolution remains
available but must be requested explicitly (`--resolution-mode flat`).

Within-family resolution uses **partial confidence-weighted soft gating** by
default (`--hierarchical-gating soft`): each subtype's estimated mass is scaled
by a calibrated confidence in [0,1] and the residual family mass goes to
`unresolved_<family>`, conserving total mass. Soft gating was validated on two
tissues (breast and lung benchmarks). The previous binary threshold gate is
still available as `--hierarchical-gating hard` (**legacy**; it over-abstains in
collinear families), and `--hierarchical-gating ungated` is a diagnostic-only
mode (no abstention). These options affect hierarchical mode only — flat and
auto solver behaviour are unchanged.

TissueResolve follows an **evidence-based decision order**: it evaluates reference
quality and query/reference compatibility, then a **Resolution Decision Layer**
(`src/tissueresolve/resolution.py`) classifies each broad family as `broad_only`,
`selected_fine`, or `full_fine` **before** fine predictions are interpreted —
driven by full-panel deconvolution reliability and within-family separability, not
by cell-level classification AUROC. Fine deconvolution is trusted only where
supported; `broad_only` families show broad mass as trusted and any fine split as
diagnostic only. Soft gating remains the final uncertainty layer. The QC-first
report shows reference quality, query compatibility, and the trusted-resolution
table before any fine predictions. TissueResolve is unified but **modality-aware**:
bulk and spatial share the reference, hierarchy, resolution decision, and soft
gating, but differ in prediction unit (sample vs spot), noise model, and
gene-weighting (bulk protocol-aware weighted NNLS; spatial marker + NB-CAR).

A compact fine-level refinement layer (`FineGranularityRefiner`) was evaluated to
push high-granularity resolution in collinear families and **failed promotion on
both breast and lung** — it lowered conditional RMSE in some settings only by
increasing spillover and false-positive subtype detection. It is **experimental, a
documented negative result, and not used by default**. Soft gating remains the
final hierarchical layer. See
[docs/FINE_GRANULARITY_REFINER_REPORT.md](docs/FINE_GRANULARITY_REFINER_REPORT.md).

## What TissueResolve does

- Estimates RNA-derived composition from bulk RNA-seq and 10x Visium data.
- Builds a reference from single-cell `.h5ad` files or saved reference
  directories.
- Supports bulk deconvolution with protocol-aware weighting, bootstrap
  uncertainty, and compatibility checks.
- Supports spatial deconvolution on Visium data using a graph-aware NB-CAR
  model.
- Reports gene overlap, separability, spillover, warnings, and methods text.
- Generates publication-style reports with figures and source data.

## Main differentiators

- Unified bulk + spatial workflow with the same reference abstraction.
- Protocol-aware reference and query compatibility checks.
- Explicit bulk mRNA-proportion output warnings.
- Spatial graph-aware Visium modeling and neighborhood statistics.
- Family-aware handling of non-separable cell types.
- Hierarchical broad→fine deconvolution (default when broad/fine labels exist)
  with explicit unresolved-mass reporting.
- Normalization-, protocol-, library-type-, and batch-aware diagnostics.
- Reproducible, family-aware colour map shared across all figures.
- Plotly HTML reports with saved figure source data.
- Optional benchmark harness comparing against external bulk and spatial tools.
- Real-data breast cancer validation harness kept separate from the default
  offline test suite.

## How TissueResolve differs from existing tools

TissueResolve is **not intended to replace every specialized method**. Instead,
it provides a unified, report-oriented, resolution-aware framework that makes
protocol compatibility, normalization, batch effects, separability, spillover,
and uncertainty **explicit** — rather than always forcing fine subtype
estimates. It is not claimed to be universally more accurate.

| Feature | TissueResolve | Bulk tools (MuSiC / Bisque / DWLS / CIBERSORTx) | Spatial tools (RCTD / cell2location / stereoscope / SPOTlight) | BayesPrism-like broad→fine |
|---|---|---|---|---|
| Bulk RNA-seq support | ✅ | ✅ | — | partial |
| Spatial transcriptomics support | ✅ | — | ✅ | partial |
| Same reference for bulk + spatial | ✅ | — | — | — |
| Protocol/library-aware decisions | ✅ | partial | partial | — |
| Normalization-aware workflow | ✅ | varies | varies | varies |
| scRNA/snRNA/mixed reference diagnostics | ✅ | — | — | — |
| Batch-effect diagnostics | ✅ | — | partial | — |
| Donor/batch-aware marker stability | ✅ | partial | — | — |
| Protocol-aware gene filtering/weighting | ✅ | partial | — | — |
| Explicit mRNA-proportion warning | ✅ | rarely | n/a | — |
| Spatial graph-aware modeling | ✅ | — | ✅ | — |
| H&E overlay reporting | ✅ | — | partial | — |
| Separability diagnostics | ✅ | — | partial | — |
| Spillover matrix/report | ✅ | — | partial | — |
| Hierarchical broad→fine mode | ✅ | — | — | ✅ |
| Unresolved family mass / abstention | ✅ | — | — | partial |
| Publication HTML report | ✅ | partial | partial | — |
| Source data for every figure | ✅ | — | — | — |
| Real-data validation harness | ✅ | — | — | — |
| Optional external benchmarking | ✅ | — | — | — |

(“partial” = available in some tools/configurations; “varies” = depends on the
specific tool. This table is a high-level orientation, not a claim that
TissueResolve outperforms these tools on any given dataset.)

## Installation

```bash
git clone https://github.com/marianaboroni/TissueResolve.git
cd TissueResolve
python3 -m venv .venv
source .venv/bin/activate
python -m pip install -U pip
python -m pip install -e ".[all]"
python -m pytest -q
```

For finer dependency control:

```bash
python -m pip install -e ".[spatial,report,realdata,benchmark]"
```

## Quickstart

**One `tissueresolve run` processes one modality** (bulk *or* spatial). To
analyse both with the same reference, run them into **separate output
directories** and then merge the summaries with `tissueresolve combine-report`.
Writing a second, different-modality run into the same `--out` is refused unless
you pass `--force` (it would overwrite the first run).

### Bulk

```bash
tissueresolve run --reference reference.h5ad --query bulk_counts.tsv \
  --out results/bulk --mode bulk --preset standard
```

### Spatial

```bash
tissueresolve run --reference reference.h5ad --query visium.h5ad \
  --out results/spatial --mode spatial --preset publication
```

Or use the spatial subcommand:

```bash
tissueresolve spatial run --visium visium.h5ad \
  --reference reference.h5ad --output results/spatial
```

### Combined bulk + spatial report

After a separate bulk run and spatial run (sharing one reference), merge them
into one report:

```bash
tissueresolve combine-report \
  --bulk-dir results/bulk --spatial-dir results/spatial \
  --out results/combined
```

This reads the two existing run directories (it does **not** re-run
deconvolution) and writes `results/combined/report.html` (+ `methods.txt`,
`warnings.json`, `run_metadata.json`). Bulk and spatial sections — and any
benchmark summaries — are kept **separate**, and the report states clearly that
it summarises two separate runs sharing a reference, **not** a single joint
bulk+spatial model.

### Hierarchical (broad → fine) deconvolution

Hierarchical mode first estimates **broad cell-type families**, then estimates
**fine subpopulations within each family**, and only trusts a subtype split
when the subtypes are demonstrably separable within their family.  Families
that are not separable are reported at the broad level as
`unresolved_<family>` rather than split into unsupported subtypes.

It works for both bulk and spatial and is enabled with
`--resolution-mode hierarchical`.  Provide explicit broad/fine annotation
columns in the reference, or a fine→broad mapping file:

```bash
# bulk, explicit broad/fine columns in the reference
tissueresolve run --reference reference.h5ad --query bulk_counts.tsv \
  --out results/bulk_hier --mode bulk \
  --resolution-mode hierarchical \
  --broad-cell-type-col broad_cell_type --fine-cell-type-col sub_cell_type \
  --preset publication

# spatial, with a mapping file (reference has only fine labels)
tissueresolve run --reference reference.h5ad --query visium.h5ad \
  --out results/spatial_hier --mode spatial \
  --resolution-mode hierarchical \
  --cell-type-hierarchy mapping.tsv --preset publication
```

`mapping.tsv` is a two-column table (`fine_cell_type`, `broad_cell_type`).
See [docs/tutorial.md](docs/tutorial.md) for the recommended reference
annotation structure and how to interpret unresolved family mass.

### Generate reports

```bash
tissueresolve report --modality bulk --results-dir results/bulk

tissueresolve report --modality spatial --results-dir results/spatial
```

## Supported CLI commands

- `tissueresolve run`
- `tissueresolve spatial run`
- `tissueresolve report`
- `tissueresolve combine-report` — merge an existing bulk run + spatial run into one report
- `tissueresolve bulk report`
- `tissueresolve spatial report`

> Note: `tissueresolve bulk run` is present in the CLI tree but not yet
> implemented; use `tissueresolve run --mode bulk` instead.
>
> One `run` processes one modality. Use **separate output directories** for bulk
> and spatial; do not write both to the same `--out` (a cross-modality write is
> refused unless you pass `--force`). Use `combine-report` for a unified summary.

## Outputs

A `tissueresolve run` writes an analysis bundle:

- `analysis_plan.json` — detected inputs, resolved mode/preset, resolution mode
- `run_metadata.json` — run provenance and resolved parameters
- `deconv/` — prediction tables (`proportions.tsv`), selected genes,
  reconstruction QC, and (hierarchical mode) family/conditional/unresolved tables
- `qc/` — per-sample/per-spot QC, Moran's I (spatial), and QC recommendations
- `methods.txt` — auto-generated methods text for the run
- `warnings.json` — surfaced warnings (estimate type, QC, non-convergence,
  protocol risk)
- `figures/` — interpretive figures (PNG + PDF/SVG + `.data.tsv` source data):
  bulk composition / proportion heatmap / reconstruction QC; spatial abundance
  maps / dominant-type map / Moran's I
- `report.html` — publication-style, **figure-driven** report generated from the
  run result (predictions, QC, methods, warnings, and the figures embedded)

`run` renders `report.html` from the in-memory result with the figures embedded.
You can also (re)generate a report from a results directory with `tissueresolve
report --modality bulk|spatial --results-dir <dir>` (see below). The richer
publication report with the full diagnostic figure set + technical appendix is
produced by the real-data validation harness
(`examples/real_breast_cancer/scripts/07_generate_reports.py`).

## Where are my results?

Open **one** file:

- **`outputs/report.html`** — the unified report with section navigation
  (executive summary, reference quality, input, bulk, spatial, hierarchical,
  resolution/spillover, benchmark, warnings, figures, methods, output files).

Supporting locations (you usually don't need to open these directly):

- `outputs/tables/`, `outputs/figures/`, `outputs/methods.txt`,
  `outputs/run_metadata.json`
- detailed sub-reports: `outputs/bulk/report.html`,
  `outputs/spatial/report.html`, `outputs/hierarchical/`
- benchmarks: `benchmarks/outputs/benchmark_summary_report.html`

Each output directory also has an `index.html` / `README.md` pointing to the
main report. You should not have to hunt through folders.

## Reports and figures

The report is split into two files (QC-first, concise):

- **`report.html`** — the concise, user-facing report: executive summary →
  reference/signature QC → input compatibility → resolution/uncertainty → final
  bulk predictions → final spatial predictions → bulk benchmark → spatial
  benchmark → warnings → methods/source-data links. One main figure per idea;
  tables collapsed; warnings summarized.
- **`technical_appendix.html`** — full separability/spillover heatmaps, spillover
  networks, spot pies, the large signature heatmap, and the complete source-data
  listing. The main report links here and back.

Figure outputs include Plotly HTML and, when `kaleido` is installed,
PDF/SVG/PNG static exports. Every figure writes a `.data.tsv` source data file,
and `figures/figure_manifest.tsv` maps each figure to its data. See
[`docs/reporting.md`](docs/reporting.md) and
[`docs/output_interpretation.md`](docs/output_interpretation.md).

## Interpretation

- Bulk outputs are RNA-derived **mRNA proportions**, not absolute cell counts.
- Spatial outputs are spot-level RNA-derived composition estimates, not direct
  single-cell counts.
- Non-separable cell-type pairs may be unreliable; family-level interpretation
  is provided when appropriate.
- In hierarchical mode, `unresolved_<family>` columns mean the family was
  detected but its subtypes could not be reliably separated — interpret those
  at the family level only, never as confident subtype fractions.
- Bootstrap uncertainty is available in publication and diagnostic modes.
- Figure colours are deterministic and consistent across all panels of a run.
  In hierarchical mode each broad family gets a distinct base colour and its
  fine subtypes use related shades; the mapping is saved to
  `cell_type_color_map.tsv` (and `color_map.json`) so reports are reproducible.

## Real-data breast cancer validation

A dedicated validation harness lives in
`examples/real_breast_cancer/` and is designed to run only when explicitly
requested.

Scripts include data download, reference building, pseudobulk generation,
bulk validation, spatial validation, and report generation.

## Benchmarking

An optional, offline-first benchmark harness lives in `benchmarks/`. It compares
TissueResolve (flat + hierarchical) against **five executable internal
baselines** (NNLS, weighted NNLS, marker-only NNLS, ridge NNLS, correlation
matcher — these are internal baselines, *not* published external tools).

**External tools require installation or imported results.** External methods
(bulk: MuSiC, BisqueRNA, DWLS, SCDC, CIBERSORTx, BayesPrism; spatial: RCTD,
cell2location, stereoscope, SPOTlight, Tangram, CARD, DestVI) are included in
the comparison **only** when they are installed locally or when you provide
their predictions via `benchmarks/import_external_results.py`. Tools that are
merely *exported* (inputs written for manual/web execution) are clearly labelled
and are **not** counted as executed benchmarks. Missing external tools are
skipped gracefully with an install hint; one missing tool never fails the run.

```bash
# run an external tool yourself, then import its predictions for a fair comparison
python benchmarks/import_external_results.py --method RCTD --modality spatial \
    --predictions rctd_predictions.tsv
python benchmarks/run_all.py --use-existing-real-data --include-imported
```

```bash
# plan only (which methods run / are skipped)
python benchmarks/bulk/run_bulk_benchmark.py --dry-run
python benchmarks/spatial/run_spatial_benchmark.py --dry-run

# toy synthetic data with ground truth (fast, fully offline)
python benchmarks/run_all.py --toy

# existing real breast-cancer harness outputs (no download)
python benchmarks/bulk/run_bulk_benchmark.py --use-existing-real-data
python benchmarks/spatial/run_spatial_benchmark.py --use-existing-real-data
```

The harness reports normalization decisions, protocol/library detection,
scRNA/snRNA/mixed reference comparisons, batch diagnostics, and flat-vs-
hierarchical comparisons, and writes a single linked
`benchmark_summary_report.html`. Bulk pseudobulk and synthetic spatial have
ground truth (accuracy is reported); real Visium has none (concordance,
spatial structure, stability, runtime, and failure modes are reported instead).
See [docs/benchmarking.md](docs/benchmarking.md). Benchmark outputs are
git-ignored and never committed.

## Documentation

See the documentation pages in `docs/`:

- `docs/tutorial.md`
- `docs/quickstart.md`
- `docs/input_formats.md`
- `docs/advanced_parameters.md`
- `docs/reporting.md`
- `docs/output_interpretation.md`
- `docs/benchmarking.md`
- `docs/normalization_and_protocols.md`
- `docs/batch_effects.md`
- `docs/library_type_references.md`

## Status

**TissueResolve is currently alpha / early-access research software. It is not
yet a fully publication-ready method until external benchmarks, multi-dataset
validation, and final API stabilization are complete.** The priority is to
stabilize what exists, not to expand the tool.

Research software / pre-release (v0.1). Scope is classified in
[`docs/V0_1_SCOPE.md`](docs/V0_1_SCOPE.md) and
[`docs/FEATURE_STATUS.md`](docs/FEATURE_STATUS.md):

- **Core (default-safe):** reference-based bulk & spatial deconvolution,
  broad/fine hierarchy with unresolved mass, solver `auto`, reference
  suitability, separability/spillover diagnostics, the QC-first report, and an
  honest basic benchmark (bulk and spatial kept separate).
- **Experimental (behind explicit flags; not default; labelled experimental):**
  state-aware three-level hierarchy (`--state-aware`), granular signatures,
  spatial multi-metric ranking, external-tool benchmark runners, the composite
  scorecard, and the synthetic state-aware benchmark.
- **Deferred / not implemented (do not assume available):** reference
  adaptation, cell-type-specific expression reconstruction, hyperparameter
  tuning, and a full BayesPrism-like Bayesian model.

> **Experimental: state-aware deconvolution runs behind `--state-aware`. It is
> not part of the default v0.1 workflow and has not been validated across real
> datasets.** State-aware outputs are experimental; full standard-report
> integration is limited.
>
> **Cell-type-specific expression reconstruction is planned/deferred and not
> implemented in v0.1.** TissueResolve produces RNA-derived composition estimates
> (cell-type-level deconvolution), not reconstructed per-cell-type expression.

## License

TissueResolve is released under the MIT License. See the top-level
[`LICENSE`](LICENSE) file for the full text; the same license is declared in
`pyproject.toml`.
