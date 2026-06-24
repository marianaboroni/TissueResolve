# Tutorial — Spatial (10x Visium) deconvolution

TissueResolve spatial estimates **spot-level RNA-derived composition** (not single-cell
counts) with an NB-CAR model. No real-Visium accuracy is claimed (real Visium has no
ground truth); accuracy statements come from synthetic gold-truth.

## 1. Default run

```bash
tissueresolve run \
  --reference reference.h5ad \
  --query visium_folder_or.h5ad \
  --out out_spatial/ \
  --mode spatial
```

Outputs: `deconv/` (spot × cell-type composition), `qc/` (per-spot QC, Moran's I,
convergence), `resolution/` (adaptive resolution + reference reliability),
`run_metadata.json`, `report.html`.

## 2. Experimental smoothing presets (opt-in)

```bash
# less global smoothing (reduces oversmoothing; experimental)
tissueresolve run ... --mode spatial --spatial-preset weak_smoothing

# edge-weighted graph (better boundary preservation; experimental)
tissueresolve run ... --mode spatial --spatial-preset edge_aware_smoothing
```

Both are **experimental**; the default keeps `lambda_spatial = 0.1`. An explicit
`--lambda-spatial` overrides the preset.

> **Not recommended** (benchmark-only): `combined_weak_edge_smoothing` (redundant with
> `weak_smoothing`) and the `state_regularized_solver_*` presets (negative results).
> They remain available for reproducibility but should not be used as standard modes.

## 3. Interpreting oversmoothing & boundary limitations

On synthetic gold truth (see `docs/PERFORMANCE_BENCHMARK_REPORT.md` §7):

- TissueResolve's default spatial solver has the **best local RMSE** and competitive
  fine accuracy, but **oversmooths more than CARD** (oversmoothing ≈1.95 vs CARD ≈0.59;
  1.0 = truth-like). This is the main spatial weakness.
- `weak_smoothing` reduces oversmoothing (e.g. breast 1.80→1.31) but does **not** reach
  CARD's level. `edge_aware_smoothing` improves boundary F1 but not the oversmoothing
  gate. Neither fully closes the gap.

Use the smoothing presets when boundary/oversmoothing behaviour matters for your
analysis, and report the oversmoothing score from QC alongside estimates.

## 4. Reading spatial QC

`qc/` includes per-spot entropy, NB log-likelihood, spatial residual, Moran's I per
cell type, convergence status, and the recorded smoothing parameters
(`lambda_spatial`, `alpha`, `edge_aware_smoothing_used`). High-entropy / low-loglik
spots are flagged as uncertain; Moran's I quantifies spatial structure per type. The
`resolution/` tables apply here too — interpret collinear fine states at the family
level.
