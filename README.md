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
- `tissueresolve bulk report`
- `tissueresolve spatial report`

> Note: `tissueresolve bulk run` is present in the CLI tree but not yet
> implemented; use `tissueresolve run --mode bulk` instead.

## Outputs

Runs produce:

- `tables/` — result tables and QC outputs
- `figures/` — Plotly figures and static exports
- `report.html` — publication-style report
- `run_metadata.json` — run provenance
- `analysis_plan.json` — selected mode and preset
- `warnings.json` — warnings and issues
- `methods.txt` — methods text for reports

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

The reporting layer embeds:

- executive summary cards
- key findings and interpretation
- a main publication figure
- detailed collapsible outputs
- warnings and QC notes

Figure outputs include Plotly HTML and, when `kaleido` is installed,
PDF/SVG/PNG static exports. Every figure writes a `.data.tsv` source data file.

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
TissueResolve (flat + hierarchical) against internal baselines and, when
installed, external bulk tools (MuSiC, Bisque, DWLS, CIBERSORTx-export) and
spatial tools (RCTD, cell2location, stereoscope, SPOTlight, Tangram). Missing
external tools are skipped gracefully with an install hint; one missing tool
never fails the run.

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

Research software / pre-release.

## License

License metadata is declared as BSD-3-Clause in `pyproject.toml`.
A top-level `LICENSE` file is not currently present in this repository.
