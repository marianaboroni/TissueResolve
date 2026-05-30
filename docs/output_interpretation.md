# Interpreting TissueResolve outputs

The single most important rule: **know what the numbers mean before you use
them.** TissueResolve labels its estimate types explicitly and refuses to
silently convert between them.

## Bulk: mRNA proportions ≠ cell fractions

`BulkDeconvResult.proportions` contains **RNA-derived mRNA proportions**
(`ESTIMATE_TYPE = "mRNA_proportion"`). Each row sums to 1. These are the
fraction of *mRNA* attributable to each cell type, **not** the fraction of
*cells*.

Cell types differ in mRNA content per cell (e.g. plasma cells, neurons and
hepatocytes are mRNA-rich), so mRNA proportions can diverge substantially from
cell fractions. To obtain cell fractions you must apply explicit mRNA-content
correction:

```python
from tissueresolve.bulk.solver import MRNAContentCorrector
corrector = MRNAContentCorrector(mrna_content_per_type)
result = deconv_bulk(bulk, ref, mrna_corrector=corrector)
# result.cell_fractions is now populated (estimate_type: cell_fraction)
# result.proportions is unchanged (still mRNA proportions)
```

Saved files carry a `# estimate_type:` header so the distinction survives
export.

## Spatial: spot-level composition ≠ cell counts

`SpatialDeconvResult.proportions` contains **spot-level RNA-derived cellular
composition estimates** (`ESTIMATE_TYPE = "spot_rna_composition"`). Each spot's
row sums to 1. These describe the RNA-derived composition of each Visium spot,
**not** the number of cells of each type.

### Smoothing is always recorded

The CAR spatial prior smooths estimates across neighbouring spots. The strength
`lambda_spatial` (and the derived mixing weight `alpha`) are recorded in
`run_metadata` and the saved `proportions.tsv` header. Niche detection's
`n_smooth` is likewise recorded. Smoothing is **never** applied without
reporting the parameter.

### Dominant type ≠ pure spot

The dominant-type map shows each spot's argmax cell type. A spot can be
"dominated" by a type at 35% — the label does not imply purity. Always read it
alongside the abundance maps and the dominant-fraction column.

## QC and separability

- QC thresholds are **heuristic** (see `qc_thresholds.md`). Flags add warnings;
  they never delete data.
- Poorly separable cell-type pairs (high Bhattacharyya coefficient) produce
  warnings, and estimates for those types should be treated as unreliable.

## Reports and figures

- Every figure produced by `tissueresolve.plotting` saves its underlying data
  as a TSV next to the image.
- HTML reports surface all warnings (estimate type, non-convergence, protocol
  risk, separability) in a box at the top and never hide failed checks.
