# Poisson GLM — Recommended Experimental Bulk Solver

**Status: recommended *experimental* solver (opt-in). NOT the default.** The default
bulk solver remains weighted-NNLS (`wNNLS`); switching the package default requires
broader real-data + external validation. On the evidence below, the Poisson GLM is
the recommended choice when a user opts into an experimental bulk solver.

Outputs remain **RNA-derived (mRNA) proportions**, never cell fractions.

## Why Poisson GLM

The default `wNNLS` minimises a Gaussian/L2 loss on L1-normalised profiles — a
noise model mismatched to RNA-seq counts (heteroscedastic, library-size-dependent).
The Poisson GLM (`experimental/nb_bulk_solver.py`) minimises a count likelihood with
an explicit library-size term, solved by a monotone multiplicative MM update:

```
μ_gn = ℓ_n · Σ_k Φ_gk θ_kn ,   min_{θ≥0} Σ_g [ μ_g − y_g log μ_g ] ,  θ renormalised to the simplex
```

## Evidence (benchmarked, gold-truth pseudobulk, donor-held-out)

Reports: `docs/NB_BULK_SOLVER_REPORT.md` (full), plus
`benchmarks/diagnostics/nb_bulk_solver_benchmark.py` and
`benchmarks/diagnostics/external_bulk_comparison.py`.

- **vs wNNLS (breast + lung, 5 seeds × 3 scenarios):** beats wNNLS on fine Pearson
  in **30/30** cells; fine 0.70→0.85 (breast) / 0.56→0.76 (lung); broad 0.87→0.99 /
  0.86→0.96; **conditional within-family RMSE −39% / −33%**; effective-N closer to
  truth; rare recall up; FP-subtype and absent mass not increased; 0 failures;
  comparable runtime. Promotion gates: **lung 8/8, breast 7/8** (the one miss is a
  rare precision/recall tradeoff, addressed below).
- **vs MuSiC / BisqueRNA (identical donor-disjoint breast inputs):** the Poisson GLM
  **leads** on fine (0.81), broad (0.98), conditional within-family RMSE (0.16,
  lowest), and effective-N calibration (16.1 vs truth 16.2). MuSiC is the closest
  competitor; BisqueRNA is conservative. (BayesPrism/DWLS/SCDC deferred — not
  installed.)

## How to enable

CLI:

```bash
tissueresolve run --reference ref.h5ad --query bulk.tsv --out out/ \
    --mode bulk --bulk-solver poisson_glm_experimental
```

Python / config:

```python
cfg = TissueResolveConfig()
cfg.bulk_solver.method = "poisson_glm_experimental"   # default "wNNLS"
result = deconv_bulk(bulk, ref, config=cfg)
```

YAML: `bulk_solver: { method: poisson_glm_experimental }`.

## Caveats and companion layers

- **Still experimental / opt-in.** Not the package default in this release.
- **Rare-subtype FPR tradeoff:** the count likelihood raises rare *recall* but can
  raise rare false positives (notably on breast). Use the opt-in **rare-detection
  layer** (`experimental/rare_detection.py`) as a post-hoc decision step — a
  marker-support + calibrated-probability gate that cuts rare FPR while preserving
  recall, moving dropped mass to `unresolved_<family>` (mass-conserving).
- **Prefer Poisson over NB.** With a Poisson warm start, NB is empirically
  indistinguishable from Poisson on these data — even with a better donor-level
  dispersion estimate (`experimental/nb_dispersion.py`). The gain is the count
  likelihood, not overdispersion; NB adds cost without benefit. Use `loss="poisson"`
  / `poisson_glm_experimental`.
- **Uncertainty:** if bootstrap CIs are requested with a GLM point solver, resampling
  currently uses wNNLS — recorded as `bootstrap_solver="wNNLS"` in metadata
  (no silent mislabel).
- **Output type:** mRNA proportions; apply `MRNAContentCorrector` only with valid
  per-type mRNA-content data.

## Companion identifiability-aware reporting

Every run now also emits (best-effort, reporting only — estimates unchanged), under
`<out>/resolution/`: `adaptive_resolution.tsv` (most strongly supported resolution per family)
and, when the raw reference is available, `reference_uncertainty.tsv` /
`state_reliability.tsv` / `family_reliability.tsv`. Use these to interpret which fine
subtypes are trustworthy vs better reported at the family level.
