# Validation plan

## Objective

Demonstrate, on real data, that TissueResolve recovers known composition in a
controlled bulk setting and produces spatially coherent, well-QC'd estimates on
a real Visium section — without altering the core algorithms.

## Design

### One reference, two uses
A single human breast-cancer single-cell atlas (CZ CELLxGENE) is built into one
`ReferenceSignature` with NB overdispersion enabled, so the *same* reference
drives both the bulk wNNLS solver and the spatial NB-CAR model.

### Bulk: pseudobulk with known ground truth
Real single cells are sampled per target composition and their raw counts
summed into pseudobulk mixtures. Three regimes stress different conditions:

- **easy** — one dominant cell type (sparse simplex),
- **medium** — 2–3 types share the mass,
- **hard** — near-uniform across many types.

Ground truth is recorded as **mRNA proportions** (fraction of total counts per
cell type), because that is exactly what the bulk deconvolver estimates.
Comparing estimates to *cell* fractions would be a category error and would
violate the "never report bulk estimates as absolute cell fractions" rule.
Cell-count proportions are still recorded in the metadata for transparency.

Metrics: Pearson, Spearman, RMSE, MAE, signed bias (overall and per cell type).

### Spatial: real Visium, qualitative validation
A real Visium breast-cancer section has no per-spot ground-truth composition.
Validation is therefore qualitative:

- estimates fit the NB model per spot,
- proportions are spatially coherent (Moran's I per cell type),
- per-spot QC and all warnings (non-convergence, separability) are surfaced,
- the smoothing strength `lambda_spatial` is recorded.

## Acceptance signals (not hard thresholds)

- Bulk: high Pearson/Spearman and low RMSE on easy/medium regimes; degradation
  on hard regimes is expected and informative.
- Spatial: positive Moran's I for spatially structured cell types; convergence;
  no hidden warnings.

## Guardrails

- Core algorithms are **not** modified by this harness.
- If validation exposes a real bug: stop, explain it, add a *failing*
  regression test under `tests/`, fix the core, rerun, then continue.
- Default tests are offline and mock all downloads.
- Legacy `chimera_v1/` and `spatcar/` are never touched.

## Limitations

- mRNA-proportion ground truth depends on the reference's count depths; it is
  self-consistent for validation but not an external cell-fraction truth.
- Spatial validation is qualitative (no spot-level ground truth).
- Census/OpenProblems availability and exact dataset contents can change;
  the manifest records what was actually used.
