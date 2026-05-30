# CLAUDE.md

## Project

You are helping develop **TissueResolve**, a scientific Python package for robust cell-type and cell-state deconvolution in:

1. bulk RNA-seq
2. 10x Visium spatial transcriptomics

TissueResolve is being built from two legacy codebases:

1. **CHIMERA**: protocol-aware bulk RNA-seq deconvolution.
2. **SpatCAR**: spatial deconvolution for 10x Visium.

The goal is not to mechanically merge two packages. The goal is to build a stable, reproducible, scientifically robust tool with:

* a unified reference layer
* explicit protocol-awareness
* clear separation between bulk and spatial assumptions
* robust QC
* uncertainty estimation
* separability diagnostics
* real-data validation
* publication-quality outputs
* complete tests

---

## Core principles

Always prioritize:

1. Scientific validity
2. Reproducibility
3. Clear software architecture
4. Explicit assumptions
5. Tests
6. Robust error handling
7. Transparent uncertainty and QC
8. Publication-quality outputs
9. Real-data validation without weakening the core algorithms

---

## Legacy scientific ideas to preserve

Preserve the best scientific ideas from both packages.

From **CHIMERA**:

* protocol-aware gene selection
* weighted NNLS
* bootstrap confidence intervals
* compatibility QC
* explicit RNA-proportion warning
* mRNA-content correction as an explicit, optional step
* marker recall and spillover diagnostics

From **SpatCAR**:

* Visium input handling
* spatial graph construction
* negative-binomial CAR model
* mismatch correction
* spatial QC
* separability diagnostics
* neighbourhood statistics
* spatial benchmark framework

---

## Non-negotiable scientific rules

Never violate these rules:

1. Never silently remove genes.
2. Never hide warnings.
3. Never report bulk estimates as absolute cell fractions unless explicit mRNA-content correction was requested and valid input was provided.
4. Never silently convert RNA proportions into cell fractions.
5. Never smooth spatial estimates without reporting the smoothing parameter.
6. Never report confident estimates for non-separable cell types without warnings.
7. Never mix bulk and spatial assumptions in a way that makes the model unclear.
8. Never generate plots without saving the underlying data.
9. Never remove, weaken, or skip tests to make the package pass.
10. Never change core algorithms during validation unless a real bug is demonstrated by a failing regression test.

---

## Development workflow

Work in stages. At each checkpoint:

1. Read the relevant project files first.
2. Inspect the current source code and tests.
3. Make only the changes requested for the current stage.
4. Run the relevant tests.
5. Report:

   * files created
   * files modified
   * tests added
   * pytest result
   * limitations or assumptions
   * whether it is safe to proceed
6. Stop. Do not automatically proceed to the next stage.

Before editing code in any new session, read:

1. `DESIGN_SPEC.md`
2. `AUDIT_AND_MIGRATION_PLAN.md`
3. this `CLAUDE.md`
4. relevant source files
5. relevant tests

Do not assume that previous context is sufficient.

---

## Codebase discovery rule

Always discover the real public API from the current codebase.

Never assume:

* function names
* class names
* argument names
* input formats
* matrix orientation
* gene ordering
* whether inputs are raw counts, CPM, log-normalized values, or AnnData objects
* whether outputs are DataFrames, dataclasses, AnnData, TSV files, or directories

Before calling or modifying any function:

1. Read the source.
2. Read the tests.
3. Check existing examples if available.
4. Prefer using existing public APIs rather than reaching into internals.
5. If the API is unclear, stop and report the ambiguity.

---

## Testing rules

The default test suite must be fast, deterministic, and offline.

Rules:

1. Default tests must not access the internet.
2. Default tests must not download public datasets.
3. Real-data tests must run only behind an explicit flag, for example:

   * `--run-real-data`
   * `TISSUERESOLVE_RUN_REAL_DATA=1`
4. Unit tests must use toy synthetic data.
5. Integration tests may use tiny local fixtures.
6. Real-data validation scripts may download data, but only when explicitly run by the user.
7. Do not mark failing tests as skipped unless there is a documented reason.
8. If a bug is found, first add a failing regression test, then fix the bug.

---

## Environment rules

Use the project virtual environment.

Before running tests or scripts, activate:

```bash
source .venv/bin/activate
```

Prefer:

```bash
python -m pytest -v
```

over:

```bash
pytest -v
```

When adding dependencies:

1. Prefer minimal dependencies.
2. Add them to `pyproject.toml`.
3. Pin or constrain versions when needed for reproducibility.
4. State the Python version targeted by the workflow.
5. Do not add large dependencies for convenience if the task can be done with existing packages.

Current target Python version: Python 3.9+ unless the project configuration specifies otherwise.

---

## Real-data validation harness

We are building a real-data validation harness under:

```text
examples/real_breast_cancer/
```

This is validation tooling only.

The goal is to use one real single-cell breast cancer reference for both:

1. bulk validation through pseudobulk mixtures with known ground truth
2. spatial validation using a real 10x Visium breast cancer dataset

The validation harness must not weaken or alter the core algorithms.

---

## Rules for real-data validation

For work under `examples/real_breast_cancer/`:

1. Do not modify core algorithms unless validation exposes a real bug.
2. If validation exposes a bug:

   * stop
   * explain the bug
   * add a failing regression test
   * fix the bug
   * rerun tests
   * then continue
3. Do not make the default test suite depend on internet access.
4. Do not hard-code fragile temporary URLs.
5. Prefer stable dataset identifiers, manifests, or official download mechanisms.
6. Save a `download_manifest.json` with:

   * dataset name
   * source
   * URL or dataset ID
   * access date
   * package versions
   * filters applied
   * downsampling performed
7. If download fails, fail with a clear, actionable error.
8. Never fail silently.
9. If a manual download is required, print exact instructions and expected file paths.
10. Keep real-data examples small enough for a first validation run when possible.

---

## Real-data validation expected structure

Use this structure:

```text
examples/real_breast_cancer/
├── README.md
├── data/
│   ├── reference/
│   ├── spatial/
│   └── derived/
├── scripts/
│   ├── 00_download_data.py
│   ├── 01_prepare_reference.py
│   ├── 02_make_pseudobulk.py
│   ├── 03_run_bulk_validation.py
│   ├── 04_run_spatial_validation.py
│   └── 05_summarize_results.py
└── outputs/
    ├── reference/
    ├── bulk/
    ├── spatial/
    └── validation_summary/
```

---

## Real-data validation outputs

The validation workflow should produce:

### Reference outputs

* `reference_summary.tsv`
* `cell_type_counts.tsv`
* `selected_annotation_column.txt`
* saved TissueResolve reference object
* warnings and metadata

### Bulk pseudobulk outputs

* `pseudobulk_counts.tsv`
* `pseudobulk_true_proportions.tsv`
* `pseudobulk_metadata.tsv`

### Bulk validation outputs

* `bulk_estimated_proportions.tsv`
* `bulk_validation_metrics.tsv`
* `bulk_per_celltype_metrics.tsv`
* `bulk_qc.tsv`
* `bulk_warnings.json`

### Spatial validation outputs

* `spatial_spot_proportions.tsv`
* `spatial_qc.tsv`
* `morans_i.tsv`
* `spatial_warnings.json`
* `spatial_run_metadata.json`

### Summary outputs

* `validation_summary.md`
* `validation_summary.json`

---

## Real-data validation metrics

For pseudobulk validation with known ground truth, compute:

* Pearson correlation
* Spearman correlation
* RMSE
* MAE
* bias
* per-cell-type RMSE

For spatial validation without ground truth, report:

* gene overlap
* number of spots
* number of genes
* number of cell types
* spot-level QC summary
* Moran’s I
* proportion entropy
* low-quality spot count
* non-separable cell-type warnings
* whether proportions are degenerate or biologically plausible

Do not claim spatial accuracy without ground truth.

---

## Bulk-specific rules

For bulk RNA-seq:

1. Outputs are RNA-derived/mRNA proportions by default.
2. Do not call them absolute cell fractions.
3. mRNA-content correction must be explicit.
4. If mRNA correction is requested but valid mRNA content is unavailable, raise a clear error or emit a prominent warning.
5. Always report:

   * gene overlap
   * selected genes
   * excluded genes
   * protocol risk
   * reconstruction metrics
   * uncertainty if bootstrap is enabled
   * warnings

---

## Spatial-specific rules

For spatial transcriptomics:

1. Outputs are spot-level RNA-derived composition estimates.
2. Do not call them direct cell counts unless explicitly calibrated.
3. Spatial smoothing must be explicitly parameterized and recorded.
4. Always report:

   * gene overlap
   * spatial graph parameters
   * smoothing parameters
   * convergence status
   * spot QC
   * residual or likelihood metrics when available
   * Moran’s I when computed
   * warnings
5. If a spatial method is not yet implemented, validation scripts must exit with a clear message rather than pretending to run.

---

## Plotting rules

When plotting is implemented:

1. Save every figure as PDF, SVG, and PNG when feasible.
2. Save the underlying plotting data as TSV or CSV.
3. Use consistent color palettes.
4. Use readable labels.
5. Include caption templates.
6. Do not prioritize visual appeal over scientific validity.
7. Never generate a figure that cannot be reproduced from saved data and metadata.

---

## Reporting rules

Reports must include:

* input summary
* reference summary
* gene overlap
* protocol or mismatch warnings
* deconvolution estimates
* QC metrics
* uncertainty where available
* separability diagnostics
* output file list
* reproducibility metadata
* methods text suitable for a manuscript

Reports must not hide warnings or failed checks.

---

## Git and safety rules

1. Do not modify legacy folders unless explicitly instructed.
2. Do not delete user-created files.
3. Do not commit generated large datasets.
4. Do not commit downloaded real-data files.
5. Keep examples, scripts, manifests, and small fixtures under version control.
6. Large data should stay under ignored directories such as:

   * `examples/real_breast_cancer/data/`
   * `examples/real_breast_cancer/outputs/`
7. Before finishing a stage, report whether `.gitignore` needs updates.

---

## Final instruction

Be critical. Do not assume that a function is scientifically valid just because it runs.

If something is scientifically ambiguous, statistically unsafe, or insufficiently tested, say so clearly and stop rather than building on a fragile assumption.