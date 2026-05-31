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

## Resolution, spillover, and families (Stage 6)

Fine cell-type panels (e.g. 32 breast-cancer subtypes) contain pairs that are
not reliably distinguishable. TissueResolve makes this explicit rather than
pretending otherwise.

- **Cell-type proportion vs mRNA proportion** — bulk estimates are *mRNA*
  proportions (RNA contributed per type), not cell fractions; see above. The
  resolution layer does not change this — it annotates it.
- **Separability** — how distinct two reference profiles are (Bhattacharyya
  coefficient; `1` = identical, `0` = orthogonal).
- **Spillover** — measured cross-type leakage: deconvolving a *pure* type's
  pseudobulk and seeing mass land on a *different* type. The **spillover
  matrix** (rows sum to 1) and per-type **spillover risk** (`1 − self-retention`)
  plus the **main leaking partner** are reported as additional outputs.
- **Resolvability classes** — `resolved`, `partially_resolved`,
  `poorly_resolved`, `unresolved` (see `qc_thresholds.md`).
- **Families** — confusable types are grouped; a family is the level at which
  estimates are trustworthy when its members are not individually resolvable.
- **Unresolved mode** — when a family is not resolvable (low separability, high
  spillover, high uncertainty), its mass is reported as `unresolved_<family>`
  instead of being confidently split into subtypes. Total mass is preserved.
- **Hierarchical view** — `broad_proportions` (group level),
  `conditional_subtype_proportions` (within a group), `absolute_subtype_proportions`
  (broad × conditional), and `unresolved_family_mass` (kept at the broad level).

### Fine vs family-level estimates, and resolution modes

- **Fine estimates** are per cell type (the deconvolver's raw output). **Family
  estimates** sum confusable subtypes into a recommended family
  (`recommended_merges.tsv`), giving a more reliable read when subtypes are not
  separable. Fine estimates are **never overwritten** — family estimates are an
  additional, safer interpretation.
- **`resolution_mode`** (`ResolutionConfig`, and `--resolution-mode` in the
  resolution analysis script) controls behaviour:
  - `none` — report only (improved warning), estimates unchanged;
  - `suggest` *(default)* — compute recommended merges, estimates unchanged;
  - `auto` — aggregate to families;
  - `hierarchical` — broad families first, subtypes only when resolvable.
- **Post-hoc aggregation vs pre-deconvolution merging.** Post-hoc aggregation
  (the default) deconvolves at fine resolution then sums subtype proportions
  into families — mass-preserving and reversible. Pre-deconvolution merging
  (`merge_reference_cell_types`) aggregates the *reference* before fitting.
  Either way the merge mapping and stage are recorded; **merging is always
  explicit, never silent.**

**Why similar cell types may be reported as a family rather than split:** if the
reference cannot separate two subtypes, a confident split would be fabricated
precision. Reporting the family (or `unresolved_<family>`) is the honest result;
use a coarser reference, a better-resolved reference, or pairwise marker
augmentation if subtype resolution is required.

## QC and separability

- QC thresholds are **heuristic** (see `qc_thresholds.md`). Flags add warnings;
  they never delete data.
- Poorly separable cell-type pairs (high Bhattacharyya coefficient) produce
  warnings, and estimates for those types should be treated as unreliable.

## Reports and figures

- **Every figure** produced by `tissueresolve.plotting` saves its underlying
  data as a `.data.tsv` (clustered barplots also save `.sample_order.tsv` and
  `.cell_type_order.tsv`). A figure is never produced without its source data.
- Figures are interactive **Plotly** HTML; static **PDF/SVG/PNG require
  `kaleido`** (`pip install kaleido`). Without kaleido, HTML + source data are
  still written and a warning is recorded — nothing fails silently.
- Generate reports with `tissueresolve bulk report` / `tissueresolve spatial
  report` / `tissueresolve report --modality ...`, or
  `tissueresolve.report.generate_report(modality, results_dir, out)`.
- Every figure carries the **estimate type** in its subtitle: bulk =
  "mRNA-derived proportion, not absolute cell fraction"; spatial = "spot-level
  RNA-derived composition, not cell counts".
- HTML reports surface all warnings (estimate type, low-confidence/uncertainty,
  non-separability, spillover, non-convergence, protocol risk) in a box and
  never hide failed checks (e.g. `converged = False` is shown in red).
- **Citing/interpreting:** report bulk values as RNA-derived mRNA proportions
  and spatial values as spot-level RNA-derived composition; read poorly
  separable / high-spillover types at the family level (see above).
