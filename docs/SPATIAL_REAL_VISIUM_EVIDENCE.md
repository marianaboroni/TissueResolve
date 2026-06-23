# Real Visium Breast Cancer — No-Ground-Truth Evidence (PART 13)

**This is NOT an accuracy benchmark.** The real 10x Visium breast-cancer section
has no per-spot cell-type ground truth. We report *evidence* metrics only.
Reading rules carried from the brief: **concordance is method agreement, not
accuracy**; **marker recovery is a proxy**; **high Moran's I / smoothness is not
automatically better** (it can mean over-smoothing — see the synthetic benchmark,
where TissueResolve over-smooths ~2×).

**Data.** `human_breast_cancer_1.h5ad` (10x Visium V1_Breast_Cancer_Block_A_Section_1),
**500 of 3798 spots** subsampled (seed 0), 36 601 genes; reference = the saved
32-type breast atlas signature. Runner: `benchmarks/spatial/run_real_visium_evidence.py`;
outputs: `benchmarks/outputs/spatial_real_visium/`.

## Evidence metrics (no ground truth)

| Method | marker recovery (proxy) | mean Moran's I | mean entropy | near-zero frac | reconstruction Pearson | runtime |
|---|---|---|---|---|---|---|
| TissueResolve_spatial_flat | **0.285** | 0.319 | 2.92 | 0.024 | 0.400 | 7 s |
| TissueResolve_spatial_hierarchical | 0.200 | **0.360** | 2.09 | 0.681 | 0.328 | 10 s |
| NNLS_per_spot | 0.149 | 0.217 | 0.72 | 0.906 | **0.448** | 15 s |

**Cross-method concordance (agreement, NOT accuracy):**

| Pair | fine Pearson |
|---|---|
| NNLS_per_spot vs TissueResolve_spatial_flat | 0.474 |
| NNLS_per_spot vs TissueResolve_spatial_hierarchical | −0.045 |
| TissueResolve_spatial_flat vs TissueResolve_spatial_hierarchical | 0.067 |

## What the evidence suggests (carefully, no accuracy claims)

1. **Marker-recovery proxy favours TissueResolve_spatial_flat** (0.285 vs NNLS
   0.149): its per-spot abundances track known cell-type marker expression more
   closely than spatial-naive NNLS. A proxy only — markers are imperfect and the
   absolute values are low (real Visium + 500 spots + a 32-type atlas is hard).
2. **NNLS-per-spot is extremely sparse on real data** (near-zero fraction 0.91,
   entropy 0.72) — it explains each spot with very few cell types. This mirrors the
   bulk audit's NNLS sparsity. TissueResolve_flat is diffuse (near-zero 0.02).
   Neither is "right" without truth; the contrast is the same solver-geometry
   effect documented in `BULK_PREDICTION_SPARSITY_AUDIT.md`.
3. **NNLS has the highest reconstruction Pearson** (0.448) — expected, since
   per-spot least squares directly optimises reconstruction; it is a
   *self-consistency* check, not accuracy, and high reconstruction with high
   sparsity can mean fitting the dominant signal only.
4. **Moran's I is highest for hierarchical, lowest for NNLS** — TissueResolve
   produces more spatially-structured maps. Per the brief, this is **not
   automatically better**: on synthetic data TissueResolve over-smooths ~2×, so
   elevated real-data Moran's I is consistent with added structure that may be
   partly over-smoothing.
5. **Hierarchical barely agrees with the other methods** (−0.05, 0.07) — combined
   with its reproducible failure on synthetic ground truth, this flags the
   spatial hierarchical path as the component to investigate, not trust.

## Honest limitations

- No ground truth → **no accuracy, no ranking by correctness**.
- 500/3798 spots (subsampled for tractable external runs); marker proxy uses the
  top-15 specificity genes per type from the reference.
- Duplicate gene symbols in the Visium var were not collapsed (overlap-based
  deconvolution tolerates this; a minor caveat for the marker proxy).
- **External methods (RCTD, cell2location) not yet run on this section** — adding
  them would give cross-method concordance against established tools (the most
  informative no-GT comparison), at ~1 h CPU each on 500 spots. Inputs are
  exported to `external_inputs/` ready to run.

## Addendum — external methods (RCTD, cell2location) + cross-method concordance

RCTD (spacexr 2.2.1) and cell2location 0.1.4 were run on the **same 500 spots**
(sc reference = 3633 atlas cells of the 32 types). Honest status
(`external_real_status.tsv`): both **executed** (RCTD 1279 s, cell2location 1309 s,
CPU). RCTD first **failed** on duplicate Visium gene symbols (`duplicate
'row.names'`) — fixed by collapsing duplicate symbols (36 601→36 591 genes) and
re-running; cell2location tolerated the duplicates (pandas), a minor inconsistency.

**Evidence panel, all methods (NO GROUND TRUTH):**

| Method | marker recovery (proxy) | Moran's I | entropy | near-zero frac | reconstruction r |
|---|---|---|---|---|---|
| cell2location | **0.293** | 0.271 | 3.09 | 0.000 | 0.400 |
| TissueResolve_spatial_flat | 0.285 | 0.319 | 2.92 | 0.024 | 0.400 |
| RCTD | 0.220 | 0.191 | 1.74 | 0.560 | **0.462** |
| TissueResolve_spatial_hierarchical | 0.200 | **0.360** | 2.09 | 0.681 | 0.328 |
| NNLS_per_spot | 0.149 | 0.217 | 0.72 | 0.906 | 0.448 |

**Cross-method concordance (fine Pearson; agreement, NOT accuracy), sorted:**

| Pair | concordance |
|---|---|
| NNLS_per_spot ↔ RCTD | 0.661 |
| TissueResolve_flat ↔ cell2location | 0.544 |
| RCTD ↔ cell2location | 0.506 |
| NNLS_per_spot ↔ TissueResolve_flat | 0.474 |
| TissueResolve_flat ↔ RCTD | 0.453 |
| NNLS_per_spot ↔ cell2location | 0.369 |
| hierarchical ↔ {cell2location, RCTD, flat, NNLS} | 0.15, 0.05, 0.07, −0.05 |

**Concordance findings (no accuracy claims):**

1. **The four non-hierarchical methods form a moderately-concordant cluster**
   (pairwise 0.37–0.66) — but no pair agrees *strongly*. On real data without
   truth, methods diverge substantially; this is a caution against over-trusting
   any single deconvolution.
2. **TissueResolve_spatial_flat sits naturally inside the established-tool cluster:**
   it concords with cell2location (0.544) and RCTD (0.453) at levels comparable to
   how RCTD and cell2location agree with each other (0.506), and it ~ties
   cell2location for the best marker-recovery proxy (0.285 vs 0.293).
3. **The hierarchical spatial path is the clear outlier** — it concords only
   0.05–0.15 with every established tool *and* with its own flat sibling (0.07).
   Triangulated with its reproducible failure on synthetic ground truth and the
   bulk hierarchical issue, this is strong convergent evidence that the
   **hierarchical (spatial) component is the part to fix**, not to trust.
4. **NNLS and RCTD are the sparser methods** (near-zero 0.91 / 0.56) and concord
   most with each other (0.661) and reconstruct expression best — a
   self-consistency pattern, not accuracy.

Figure: `figures/figM_real_concordance` (full pairwise heatmap). Tables:
`spatial_real_metrics_all.tsv`, `method_concordance_{fine,family}_all.tsv`,
`external_real_status.tsv`. Wrappers: `methods/run_rctd_real.R`,
`methods/run_cell2location_real.py`; scorer `score_real_visium_external.py`.

## Output completeness / run success

All three internal methods executed and produced complete spots×types output
(`output_complete=True`). Full panel incl. Geary's C, dominant fraction,
per-method tables: `spatial_real_metrics.tsv`, `method_concordance_{fine,family}.tsv`,
figure `figures/figL_real_evidence`.

---

*Source: `benchmarks/outputs/spatial_real_visium/`. Runner:
`benchmarks/spatial/run_real_visium_evidence.py`. No core algorithm modified;
nothing committed. NOT an accuracy benchmark.*
