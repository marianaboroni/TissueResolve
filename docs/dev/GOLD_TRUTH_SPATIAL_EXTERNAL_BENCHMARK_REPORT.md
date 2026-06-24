# Gold-Truth Spatial External Tool Benchmark Report

Run 2026-06-23. Spatial deconvolution benchmark of TissueResolve against external
spatial tools on a **synthetic** spatial scenario with known ground truth. No
real Visium accuracy is claimed (real Visium has no ground truth).

Reproduce:

```bash
python benchmarks/spatial/run_synthetic_spatial.py --run-real-data       # TissueResolve + NNLS baseline
python benchmarks/spatial/export_spatial_external.py --run-real-data      # export same seed-0 scenario
Rscript benchmarks/spatial/methods/run_card_synthetic.R                   # CARD
# RCTD + cell2location predictions were produced in prior sessions on the same seed-0 scenario
python benchmarks/spatial/score_external_synthetic.py                     # score + merge
```

## 1. Dataset and truth design

- Reference: real breast-cancer single-cell atlas
  (`examples/real_breast_cancer/data/reference/breast_cancer_sc_reference.h5ad`),
  donor-disjoint split (20 reference donors, 21 held-out query donors).
- Spatial pseudo-spots are realized from **held-out query-donor cells** at
  count level (3022 reference cells, 30 cell types, 400 spots, 40 cells/spot,
  20×20 grid, seed 0). Ground truth is the per-spot **RNA-derived (mRNA)
  proportion** — the same estimate type all scored methods produce.
- Truth, coordinates, and domain labels saved under
  `benchmarks/outputs/spatial_synthetic/` and
  `.../external_inputs/spatial/`.

## 2. Spatial-like scenarios

The synthetic generator embeds, in a **single structured scenario**, the core
spatial-difficulty features: a **sharp domain boundary**, a **smooth gradient**,
and a **rare spatial niche**, with mixed-composition spots. True mean Moran's I
≈ 0.30 (genuine spatial structure).

**Limitation (honest):** this is ONE scenario on ONE dataset at ONE seed. The
full multi-scenario × multi-dataset matrix in the task design (10 scenarios ×
breast+lung) was **not** executed — see §5. The benchmark is therefore a
single-scenario probe, not the full matrix.

## 3. Tools attempted

TissueResolve_spatial_flat, TissueResolve_spatial_hierarchical, NNLS-per-spot
(spatial-naive control), RCTD/spacexr, cell2location, CARD, SPOTlight, Tangram,
DestVI.

## 4. Tools executed

| method | status | runtime | notes |
|---|---|---|---|
| TissueResolve_spatial_flat | executed | 10.2 s | NB-CAR spatial, default lambda |
| TissueResolve_spatial_hierarchical | executed | 32.5 s | broad→fine + soft gating |
| NNLS_per_spot | executed | 5.3 s | spatial-naive control |
| RCTD (spacexr) | executed | 2078 s | prior session, same seed-0 scenario |
| cell2location | executed | 1561 s | prior session (250/6000 epochs, CPU/MPS) |
| CARD | executed | 34.6 s | this session (`run_card_synthetic.R`) |

## 5. Tools failed / skipped and why

| method | status | reason |
|---|---|---|
| SPOTlight | skipped | not installed (Bioconductor; Seurat absent) |
| Tangram | skipped | not installed |
| DestVI | skipped | not installed |
| CellTrek / NovoSpaRc | not attempted | not installed |
| Redeconve | not attempted | not installed |

RCTD/spacexr is **not currently importable** in the R environment
(`requireNamespace('spacexr')` → FALSE); its prediction file persists from a
prior session on the identical deterministic seed-0 scenario and is scored as
such (explicitly flagged). cell2location ran in a dedicated `c2l_py39` env.

**Runtime infeasibility at scale:** RCTD (~35 min) and cell2location (~26 min)
per run make the full 10-scenario × 2-dataset external matrix infeasible in a
single session (≈10+ h of external compute). This is the binding constraint on
spatial benchmark completeness, recorded honestly rather than hidden.

## 6. Input harmonization

All methods consume the SAME exported scenario: reference counts (genes×cells) +
cell metadata, spot counts (genes×spots), coordinates, and the shared truth
table. CARD reads via `run_card_synthetic.R`; cell2location via the c2l env
runner; predictions are row-normalized and reindexed to the shared fine-label
space before scoring.

## 7–8. Broad and fine results

| method | fine Pearson | broad Pearson | dominant acc |
|---|---|---|---|
| **TissueResolve_spatial_flat** | **0.799** | 0.933 | **0.708** |
| cell2location | 0.681 | 0.819 | 0.528 |
| NNLS_per_spot | 0.678 | 0.913 | 0.507 |
| CARD | 0.674 | 0.859 | 0.460 |
| RCTD | 0.538 | **0.971** | 0.428 |
| TissueResolve_spatial_hierarchical | 0.164 | −0.186 | 0.020 |

- **TissueResolve_spatial_flat has the highest fine Pearson and dominant
  accuracy of all methods**, including RCTD/cell2location/CARD.
- **RCTD has the highest broad Pearson (0.971)**; TissueResolve_flat is second
  (0.933).
- `hierarchical` abstains heavily (unresolved mass) → low raw scores by
  construction, as in the bulk benchmark.

## 9–11. Conditional / rare / spillover

The cross-method merge reports fine/broad accuracy, dominant accuracy, local
RMSE, oversmoothing, and domain recovery. Per-family conditional, rare-niche
sensitivity/precision, and pairwise-spillover tables were computed for the
TissueResolve modes in `spatial_synthetic_metrics.tsv` but were **not** uniformly
recomputed for every external tool in this merge; they are therefore reported as
**partial** and not used for cross-method claims. The rare spatial niche is
embedded in the scenario and contributes to the fine-accuracy numbers above.

## 12. Spatial structure preservation

| method | local RMSE | oversmoothing score | domain recovery ARI |
|---|---|---|---|
| TissueResolve_spatial_flat | **0.026** | 1.95 | 0.271 |
| cell2location | 0.032 | 0.96 | 0.295 |
| CARD | 0.033 | **0.59** | 0.253 |
| RCTD | 0.041 | 0.87 | **0.313** |
| NNLS_per_spot | 0.040 | 0.42 | 0.270 |

- **TissueResolve_flat has the lowest local RMSE (0.026)** — best local accuracy.
- **TissueResolve_flat oversmooths the most (1.95; >1 ⇒ smoother than truth)** —
  the spatial NB-CAR prior smooths more than RCTD/CARD/cell2location. This is the
  one structural dimension where TissueResolve is clearly weaker.
- RCTD recovers domains best (ARI 0.313); TissueResolve_flat is mid-pack (0.271).

## 13. Reconstruction accuracy

Held-out-gene / spot-level expression reconstruction was **not** computed
uniformly across external tools in this run and is not reported cross-method.

## 14. Sparsity / effective-N

Reported for TissueResolve modes in the internal spatial metrics; not uniformly
available for external tools → not used for cross-method claims.

## 15. Runtime

TissueResolve_flat 10.2 s and CARD 34.6 s are practical; RCTD (~35 min) and
cell2location (~26 min) are 1–2 orders of magnitude slower on this scenario.

## 16. Where TissueResolve performs well

- Highest **fine Pearson** (0.799) and **dominant accuracy** (0.708) of all
  methods.
- Lowest **local RMSE** (0.026).
- Competitive **broad Pearson** (0.933, second to RCTD).
- Fast (10 s).

## 17. Where TissueResolve underperforms

- **Oversmoothing**: TissueResolve_flat (1.95) smooths more than truth and more
  than RCTD/CARD/cell2location. This is the clearest weakness.
- **Domain recovery ARI**: RCTD (0.313) > TissueResolve_flat (0.271).
- **Broad Pearson**: RCTD (0.971) > TissueResolve_flat (0.933).
- `hierarchical` raw spatial accuracy is low (expected — abstention).

## 18. Supported claims

- On this synthetic gold-truth spatial scenario, **TissueResolve (flat spatial
  mode) is competitive with, and on fine-level accuracy leads, established
  spatial tools (RCTD, cell2location, CARD)**.
- TissueResolve achieves the best per-spot fine accuracy and local RMSE while
  remaining fast.

## 19. Unsupported claims

- **NOT** claimed: TissueResolve is best overall (RCTD leads broad Pearson and
  domain recovery; TissueResolve oversmooths).
- **NOT** claimed: real Visium spatial accuracy (no ground truth).
- **NOT** claimed: results generalize across many scenarios/datasets — only ONE
  synthetic scenario on ONE dataset at ONE seed was scored.

## 20. Whether algorithmic improvement is justified (Phase 3 gate)

**Phase 3 (Redeconve-inspired spatial regularization) is NOT triggered.** Two
independent gate conditions both block it:

1. **TissueResolve spatial does not underperform.** The flat spatial mode has
   the best fine Pearson, dominant accuracy, and local RMSE among all tools. The
   trigger ("only if TissueResolve spatial underperforms") is not met.
2. **The spatial benchmark is incomplete.** Only one synthetic scenario, one
   dataset, one seed was scored; the rule explicitly forbids starting Phase 3
   on an incomplete spatial benchmark.

The one genuine weakness is **oversmoothing**. A similar-state / collinearity-
aware regularizer is more about stabilizing collinear fine subtypes than about
reducing oversmoothing, so it is not clearly indicated even for that weakness.
If pursued in future, it should be evaluated against the multi-scenario matrix
(once external-tool runtime is amortized) with the promotion gates from the
task spec, and it must not increase oversmoothing further.
