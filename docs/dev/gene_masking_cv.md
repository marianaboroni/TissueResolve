# Gene-masking cross-validation

Gene-masking CV lets TissueResolve choose parameters (solver backbone, marker
count, regularization) **without external ground truth**, by testing how well a
configuration reconstructs *held-out* genes.

## Idea

1. Hide a fraction (default 15%) of the shared genes.
2. Deconvolve using only the remaining genes.
3. Reconstruct the held-out genes as `proportions × reference_signature`.
4. Score reconstruction (Pearson / RMSE) against the observed held-out genes.

A configuration that reconstructs unseen genes well *generalises*; one that only
fits the training genes does not. This avoids selecting parameters that overfit
one panel or inflate spillover.

## Used by `solver="auto"`

`solver=auto` scores each backbone (NNLS, weighted NNLS, marker NNLS, ridge
NNLS) by gene-masking CV and picks the best by a multi-objective score
(`CV reconstruction − conditioning penalty`) — **not** in-sample reconstruction
alone. The selected solver, its CV score, and the comparison table are recorded
in `selected_solver.json` / `solver_comparison.tsv` and in the report.

On the real breast-cancer pseudobulk, auto selects full-gene **NNLS**
(masked-gene Pearson ≈ 0.99) and reaches fine-level Pearson **0.768**, matching
the strongest baseline while keeping the resolution-aware layer available.

## API / CLI

```python
from tissueresolve.validation.gene_masking import compare_solvers_by_gene_masking
```

```bash
tissueresolve run ... --solver auto          # default
tissueresolve run ... --solver nnls          # force a backbone
```

Functions: `split_genes_for_masking`, `reconstruct_masked_genes`,
`score_masked_gene_reconstruction`, `run_gene_masking_cv_bulk`,
`compare_solvers_by_gene_masking`.

## Limitations

- Reconstruction quality is necessary but not sufficient for proportion
  accuracy; it is one signal among several.
- Spatial gene-masking CV and full parameter sweeps (lambda_spatial, marker
  count) are planned extensions.
