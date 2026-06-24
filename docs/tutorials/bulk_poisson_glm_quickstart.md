# Tutorial — Bulk deconvolution with the Poisson GLM solver

TissueResolve estimates **RNA-derived (mRNA) proportions** by default, not absolute
cell fractions. This tutorial shows the default `wNNLS` solver and the recommended
experimental **Poisson GLM** solver.

## 1. Default run (wNNLS)

```bash
tissueresolve run \
  --reference reference.h5ad \
  --query bulk_counts.tsv \
  --out out_wnnls/ \
  --mode bulk
```

`bulk_counts.tsv` is genes × samples (raw or normalised counts). The default solver is
weighted NNLS.

## 2. Recommended experimental: Poisson GLM

```bash
tissueresolve run \
  --reference reference.h5ad \
  --query bulk_counts.tsv \
  --out out_poisson/ \
  --mode bulk \
  --bulk-solver poisson_glm_experimental
```

Python / config:

```python
from tissueresolve.config import TissueResolveConfig
from tissueresolve.api import deconv_bulk
cfg = TissueResolveConfig()
cfg.bulk_solver.method = "poisson_glm_experimental"   # default is "wNNLS"
result = deconv_bulk(bulk_df, reference, config=cfg)
```

## 3. Output files

Under `--out`:

```
deconv/                 estimated proportions (mRNA-derived) + coverage R²
qc/                     per-sample QC, reconstruction, warnings
resolution/             adaptive_resolution.tsv + reference reliability (if reference is h5ad)
run_metadata.json       resolved parameters, solver, likelihood, convergence
report.html             human-readable summary
```

`run_metadata.json` records `solver`, `likelihood` (`poisson`/`nb`), iterations,
convergence, and `estimate_type = "mRNA_proportion"`.

## 4. Interpreting RNA-derived proportions

Values are the fraction of **mRNA** attributable to each cell type, **not** the
fraction of cells. Cell types with high per-cell mRNA content are over-represented
relative to their cell count. To convert to cell fractions you must apply explicit
mRNA-content correction (`MRNAContentCorrector`) with valid per-type content data —
otherwise keep the RNA-derived interpretation.

## 5. Uncertainty / confidence

- Enable bootstrap CIs via the preset/config (`bootstrap.n_bootstrap > 0`). Note: with
  a GLM point solver, bootstrap resampling currently uses wNNLS and is recorded as
  `bootstrap_solver="wNNLS"` in metadata.
- Read `resolution/adaptive_resolution.tsv` to see which fine states are trustworthy
  vs better interpreted as grouped/broad/unresolved (see the adaptive-resolution
  tutorial).

## 6. When to use Poisson GLM vs wNNLS

- **Use Poisson GLM** when subpopulation accuracy matters (fine/broad correlation,
  conditional within-family RMSE, rare recall) — it improved all of these over wNNLS on
  donor-disjoint breast + lung benchmarks (see `docs/PERFORMANCE_BENCHMARK_REPORT.md`).
- **Keep wNNLS** when you need the validated package default, or for reproducibility
  with prior runs. The default has not been changed pending broader validation.
- Prefer `poisson_glm_experimental` over `nb_glm_experimental`: NB did not improve over
  Poisson (the gain is the count likelihood, not overdispersion).
