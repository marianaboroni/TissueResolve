# TissueResolve

Unified cell-type and cell-state deconvolution for bulk RNA-seq and 10x Visium
spatial transcriptomics.

TissueResolve integrates two complementary algorithms:

- **Bulk deconvolution** — protocol-aware weighted NNLS with bootstrap CIs,
  derived from CHIMERA.
- **Spatial deconvolution** — negative-binomial CAR model for 10x Visium,
  derived from SpatCAR.

Both workflows share a common reference construction pipeline, gene selection
framework, separability diagnostics, QC system, and reporting layer.

## Important: what the outputs represent

**Bulk** — `BulkDeconvResult.proportions` contains **mRNA proportions**, not
cell fractions.  In tissues where cell types differ substantially in mRNA
content per cell (plasma cells, neurons, hepatocytes), these quantities differ.
Use `bulk.solver.MRNAContentCorrector` to convert if mRNA content data are
available.

**Spatial** — `SpatialDeconvResult.proportions` contains **spot-level
RNA-derived cellular composition estimates**, not direct single-cell counts
unless explicitly calibrated.

## Status

Stages 0–5 are implemented: shared reference layer, protocol layer, the full
bulk workflow (wNNLS + bootstrap + QC + CLI), the full spatial workflow
(NB-CAR model + graph + QC + neighbourhood + benchmark + CLI), and the
unification layer (plotting, HTML reports, methods text, docs, compliance
tests). Real-data validation is a **separate, offline-by-default** harness and
is not part of the default test suite.

## Installation

```bash
# Core (bulk only)
pip install tissueresolve

# With spatial dependencies (scanpy, scikit-learn)
pip install "tissueresolve[spatial]"

# Full install including plotting (matplotlib) and report generation
pip install "tissueresolve[all]"
```

## Development setup

```bash
git clone <repo>
cd TissueResolve
python -m venv .venv && source .venv/bin/activate
python -m pip install -e ".[all]"
```

## Running tests

The default suite is fast, deterministic, and **offline** — it never accesses
the network or downloads datasets.

```bash
source .venv/bin/activate
python -m pytest -v                 # full suite
python -m pytest tests/spatial -v   # one area
python -m pytest -k compliance -v   # scientific-rule compliance tests
```

## Bulk quickstart (toy data)

```python
import numpy as np, pandas as pd
import tissueresolve as tr
from tissueresolve.results import ReferenceSignature

genes = [f"GENE_{i:03d}" for i in range(60)]
cell_types = ["Tcell", "Bcell", "Myeloid"]

# Toy L1-normalised reference signature (genes x cell types)
rng = np.random.default_rng(0)
phi = rng.exponential(1.0, (len(genes), len(cell_types)))
phi /= phi.sum(axis=0, keepdims=True)
ref = ReferenceSignature(gene_names=genes, cell_types=cell_types, phi=phi)

# Toy bulk counts (genes x samples)
bulk = pd.DataFrame(
    rng.poisson(50, (len(genes), 4)),
    index=genes, columns=[f"sample_{i}" for i in range(4)],
)

result = tr.deconv_bulk(bulk, ref)          # BulkPipelineResult
print(result.deconv.proportions)            # mRNA proportions (NOT cell fractions)
print(result.deconv.ESTIMATE_TYPE)          # 'mRNA_proportion'
```

## Spatial quickstart (toy data)

```python
import numpy as np
import tissueresolve as tr
from tissueresolve.spatial.benchmark import simulate_visium

ds = simulate_visium(n_spots=80, n_types=3, n_genes=40, seed=0)
ref = ds.to_reference()
Y = ds.dense_counts()

result = tr.deconv_spatial(
    Y, ref, ds.array_row, ds.array_col, ds.lib_sizes, ds.gene_names,
    marker_genes=list(ds.gene_names),
)
print(result.deconv.proportions.head())     # spot-level RNA-derived composition
print(result.deconv.lambda_spatial)         # smoothing parameter (always recorded)

# Figures (each saves its underlying data) and an HTML report:
figs = tr.plot_results(result, "out/figs", array_row=ds.array_row, array_col=ds.array_col)
tr.generate_report(result, "out/report.html", figures=figs)
```

Or via the CLI: `tissueresolve spatial run --visium ... --reference ... --output ...`

## Reports & figures (publication layer)

TissueResolve produces an interactive HTML report plus publication figures for
both modalities. Figures are **Plotly** (interactive HTML); static PDF/SVG/PNG
are written when `kaleido` is installed, otherwise HTML + source data are still
saved and a warning is recorded. **Every figure writes its source data as a
`.data.tsv`** next to it.

```bash
pip install "tissueresolve[report]"          # plotly + kaleido + jinja2

# Generate a report from a results directory (tables/ + figures/):
tissueresolve bulk report --results-dir results/bulk --out results/bulk/report.html
tissueresolve spatial report --results-dir results/spatial
tissueresolve report --modality bulk --results-dir results/bulk
```

In Python: `from tissueresolve.report import generate_report` and
`tissueresolve.plotting.{bulk_plots,spatial_plots,separability_plots,spillover_plots}`.
Reports surface all warnings (estimate type, low confidence, non-separability,
spillover, non-convergence) and never hide failed checks. **Static export
requires kaleido**; without it you still get interactive HTML and source data.

## Resolution-aware handling of confusable cell types

Fine cell-type panels contain pairs that aren't reliably separable. Instead of
telling users to "merge manually", TissueResolve generates a **machine-readable
recommendation**: it groups confusable types into named merge families
(`recommended_merges.tsv`, `cell_type_families.tsv`) and can produce
**family-level estimates** as a safer interpretation. Fine estimates are never
overwritten and merging is always explicit and recorded.

```python
from tissueresolve.reference.resolution import (
    recommend_cell_type_merges, assign_resolution_families)
from tissueresolve.reference.hierarchy import aggregate_predictions_by_family
mapping = assign_resolution_families(separability_report, cell_types)
family_props = aggregate_predictions_by_family(fine_props, mapping)  # mass-preserving
```

`resolution_mode` (`none`/`suggest`/`auto`/`hierarchical`, default `suggest`)
controls whether merges are only recommended or applied. See
`docs/output_interpretation.md`.

## Architecture

See `DESIGN_SPEC.md` for the full architecture specification and
`AUDIT_AND_MIGRATION_PLAN.md` for the migration plan from the CHIMERA and
SpatCAR legacy packages.

## Non-negotiable output rules

- Bulk estimates are always labelled as **mRNA proportions**.
- Spatial estimates are always labelled as **spot-level composition estimates**.
- Warnings are never hidden.
- Genes are never silently removed.
- Smoothing parameters are always recorded in output metadata.
- Estimates for non-separable cell-type pairs are always flagged.

## License

BSD 3-Clause
