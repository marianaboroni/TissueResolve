# Library-type references: scRNA, snRNA, and mixed

The single-cell reference's **library type** can systematically affect
deconvolution, and it is a common, real failure mode that TissueResolve makes
explicit.

## scRNA-seq vs snRNA-seq

- **scRNA-seq** profiles whole cells; cytoplasmic and 3'-biased transcripts are
  well represented.
- **snRNA-seq** profiles nuclei; it captures more intronic/nascent transcripts
  and depletes some cytoplasmic mRNAs, so per-gene composition differs even for
  the same cell type.

Deconvolving a query with a reference of a different library type can bias
estimates. The size of the effect is gene- and cell-type-dependent.

## Mixed scRNA/snRNA references

Mixing scRNA and snRNA cells in one reference is sometimes unavoidable, but it
can confound cell type with library type (e.g. if one cell type comes mostly
from snRNA). TissueResolve flags this rather than hiding it.

## Providing library metadata

Add a library/suspension column to `adata.obs` (any of: `library_type`,
`assay`, `technology`, `suspension_type`, `protocol`, `cell_or_nucleus`) and,
for the benchmark, point to it with `--library-type-col`. Detection is
automatic from these columns when present.

## How TissueResolve evaluates library-type effects

The benchmark (`benchmarks/shared/library_type.py`) detects the reference
library type and, where metadata are missing, builds **synthetic** scRNA /
snRNA / mixed references as controlled stress tests (clearly labelled as
synthetic, not biological validation — see `simulate_library_shift`). It then
reports how library type affects gene overlap, marker stability, separability,
spillover, family- and fine-level accuracy, unresolved mass, and runtime, and
compares library-aware vs library-naive runs.

## Interpreting library-type results

A method that is robust to library type shows small accuracy/marker-stability
changes across scRNA/snRNA/mixed references. Large changes indicate the
reference library type matters for your query and should be matched (or handled
with protocol-aware weighting). For mixed references, check the
cell-type-by-library-type confounding table before trusting subtype estimates.
