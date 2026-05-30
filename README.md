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

## Installation

```bash
# Core (bulk only)
pip install tissueresolve

# With spatial dependencies
pip install "tissueresolve[spatial]"

# Full install including plotting and report generation
pip install "tissueresolve[all]"
```

For development:

```bash
git clone <repo>
cd TissueResolve
pip install -e ".[all]"
```

## Quick start

```python
# Bulk
import tissueresolve as tr

# Spatial
import tissueresolve as tr
```

Full tutorials in `examples/`.

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
