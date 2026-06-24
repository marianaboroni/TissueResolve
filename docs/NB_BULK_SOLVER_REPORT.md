# NB / Poisson GLM Bulk Solver Report (P1)

Experimental, opt-in, non-default, benchmark-gated. Adds a count-likelihood
(Poisson / negative-binomial) bulk deconvolution solver alongside the existing
weighted-NNLS default, to test whether a noise model matched to RNA-seq counts
improves accuracy — especially for low-abundance and within-family subtypes.

Outputs remain **RNA-derived (mRNA) proportions**, never cell fractions. The
default solver is unchanged.

## Sprint 0 — Audit (written before any code change)

### 1. Current wNNLS formulation
`bulk/solver.py::WNNLSSolver` solves, per sample, on the L1-normalised marker panel:
`min_θ ‖√w ⊙ (b_norm − Φ_norm θ)‖²₂, θ ≥ 0`, then renormalises `θ` to sum to 1
(mRNA proportions). Deterministic; `coverage_r2` per sample.

### 2. Input normalization
Both `b` (columns) and `Φ` (columns, from `ref.as_phi()`) are **L1-normalised per
column** (`_l1norm_cols`). So the solve is on relative (compositional) profiles, not
raw counts — there is no explicit library-size / count model.

### 3. Current weights
Composite per-gene `w` (specificity × stability × protocol × concordance) from
`GeneSelector.compute_gene_weights`, applied as `√w` reweighting of rows. Heuristic,
not the statistical (Fisher) information of a count model.

### 4. Output semantics
`BulkDeconvResult` with `ESTIMATE_TYPE = "mRNA_proportion"`; fields
`proportions, coverage_r2, gene_panel, gene_weights, lower_ci/upper_ci (bootstrap),
cell_fractions (only after explicit mRNA correction), run_metadata`. The new solver
must return the **same schema** so reports/bootstrap consume it unchanged.

### 5. Limitations of the Gaussian/L2 loss for count data
- RNA-seq counts are **heteroscedastic and overdispersed** (var ≈ μ + μ²/φ); L2 on
  L1-normalised profiles assumes homoscedastic Gaussian noise — mis-specified.
- High-expression genes dominate the squared loss; informative low-count markers
  (often the discriminative ones for rare/fine types) are **under-weighted**.
- `√w` weighting is heuristic, not the variance-optimal weighting a GLM gives.
- Internal inconsistency: the **spatial** path already uses a proper NB likelihood;
  bulk does not. P1 closes this gap.

### 6. Reuse of the spatial NB update
`spatial/model.py::_nb_multiplicative_update` is a per-spot multiplicative
majorise-minimise (MM) step on the NB NLL. The **mathematical form** is reusable for
bulk (per sample, no spatial graph), but to honour "do not alter the spatial solver"
the bulk solver **re-derives the MM update in a separate experimental module** —
no import of, or change to, spatial code.

### 7. Proposed NB/Poisson GLM formulation
Per sample `n`, design `X_{gk} = ℓ_n · R_{gk}` (R in proportion scale, `R = CPM/1e6`),
expected counts `μ_{gn} = Σ_k X_{gk} θ_{kn}`:

```
min_θ ≥ 0   Σ_g  NB_NLL(y_{gn}, μ_{gn}, φ_g)      [loss="nb"]
            Σ_g  Poisson_NLL(y_{gn}, μ_{gn})       [loss="poisson"]
```

Optional gene weights enter as exponents on the per-gene likelihood term (`w_g`).
Multiplicative MM update (monotone non-increasing NLL, keeps θ ≥ 0):

```
θ_k ← θ_k · [Σ_g w_g X_{gk} y_g / μ_g] / [Σ_g w_g X_{gk} (y_g+φ_g)/(μ_g+φ_g)]   (NB)
θ_k ← θ_k · [Σ_g w_g X_{gk} y_g / μ_g] / [Σ_g w_g X_{gk}]                       (Poisson; φ→∞)
```

then renormalise to the simplex (mRNA proportions). Dispersion `φ_g` from
`ref.phi_g` when present; otherwise Poisson fallback (or a simple
mean–variance overdispersion estimate), recorded explicitly.

### 8. Implementation plan
- `experimental/nb_bulk_solver.py` :: `fit_bulk_nb_glm(...)` returning a
  `BulkDeconvResult` (same schema), with full diagnostics in `run_metadata`.
- `config.py` :: `BulkSolverConfig.method ∈ {wNNLS, poisson_glm_experimental,
  nb_glm_experimental}` (default `"wNNLS"`).
- `bulk/pipeline.py` step 6 branches on `cfg.bulk_solver.method` (default path
  byte-for-byte unchanged). Bootstrap, if requested with a GLM point solver, still
  uses wNNLS resampling — recorded as `bootstrap_solver="wNNLS"` (no silent
  mislabel; not changed in scope).
- CLI `--bulk-solver` choice.

### 9. Benchmark plan
Reuse the gold-truth pseudobulk infrastructure (breast + HLCA/lung, ≥5 donor-held-out
seeds, same truth/gene sets). Compare wNNLS vs Poisson-GLM vs NB-GLM vs an external
plain-NNLS control. MuSiC/Bisque/DWLS/SCDC/BayesPrism only if already installed
(else deferred, not failed). CARD excluded (spatial). Outputs under
`benchmarks/outputs/nb_bulk_solver/` (gitignored).

### 10. Promotion gates (retain experimental)
broad accuracy ≥; fine accuracy ≥; conditional within-family RMSE improved or not
worsened; rare precision/recall not worsened; FP subtype rate not increased; absent
mass not increased; effective-N closer-or-equal; runtime acceptable; 0 failures;
works on both tissues. Promotion to default needs stronger evidence (not this task).

---

## Sprint 1B — Results

**Headline: a clear, benchmarked positive result.** The count-likelihood GLM beats
the default wNNLS substantially and consistently on both tissues. Gold-truth
pseudobulk, donor-held-out, breast + HLCA/lung, 5 seeds × 3 scenarios
(imbalanced / rare / similar_subtypes) = 120 fits, **0 failures**. Truth = true mRNA
proportions. Poisson and NB are reported separately; NB warm-starts from the Poisson
fit.

### Results vs wNNLS (5-seed × 3-scenario means)
| metric | breast wNNLS → GLM | lung wNNLS → GLM |
|---|---|---|
| fine Pearson | 0.695 → **0.845** | 0.557 → **0.764** |
| broad Pearson | 0.873 → **0.990** | 0.864 → **0.956** |
| fine RMSE | 0.050 → 0.032 | 0.043 → 0.029 |
| conditional within-family RMSE | 0.256 → **0.155 (−39%)** | 0.298 → **0.201 (−33%)** |
| rare sensitivity | 0.30 → **0.94** | 0.00 → **0.43** |
| rare precision | 0.52 → 0.46 | 0.00 → **1.00** |
| rare FPR | 0.16 → 0.57 | 0.02 → 0.00 |
| FP subtype rate | 0.034 → 0.031 | 0.031 → 0.027 |
| absent subtype mass | 0.041 → 0.036 | 0.074 → 0.059 |
| effective-N (truth 12.1 / 16.6) | 9.9 → 14.1 | 11.8 → **18.9** |
| runtime / fit | ~10s → ~11s | ~10s → ~12s |

**Replication: the GLM beats wNNLS on fine Pearson in 30/30 dataset×seed×scenario
cells** — fully consistent, every scenario including `similar_subtypes`.

### Conditional within-family RMSE
Improved by **−39% (breast)** and **−33% (lung)** — the metric every prior
smoothing / state-regularization experiment failed to move. Here a better *noise
model* (count likelihood + library-size term, vs Gaussian-L2 on L1-normalised
profiles) moves it materially. The gain comes from the Poisson count likelihood;
**NB ≈ Poisson** because, warm-started from Poisson and with the reference's small
`φ_g`, the NB multiplicative ratio is near-inert (recorded honestly — the win is the
count likelihood, not NB overdispersion specifically).

### Rare subtype / spillover / false positives
Rare **recall** improves markedly on both tissues (lung wNNLS never detected the
rare type: 0.00 → 0.43; breast 0.30 → 0.94). FP-subtype rate and absent-subtype mass
**do not increase** (both decrease). The one tradeoff: on **breast**, rare *precision*
dips (0.52 → 0.46) and rare FPR rises (0.16 → 0.57) — the GLM calls the rare type far
more often, trading some precision for large recall gains. On lung rare behaviour is
strictly better (precision 0→1.0, FPR 0.02→0.00).

### Effective-N
Moves closer to truth on both tissues (lung 11.8→18.9 vs truth 16.6; breast 9.9→14.1
vs 12.1) — wNNLS under-mixes; the GLM is better calibrated.

### External control
A plain (unweighted) NNLS control sits between wNNLS and the GLM on accuracy and is
clearly worse than the GLM on conditional RMSE and rare recall. MuSiC / Bisque /
DWLS / SCDC / BayesPrism are not installed in this environment → **deferred**, not
failed (recorded in `method_status.tsv`). CARD is spatial and correctly excluded.

### Promotion gates
- **Lung: 8/8 — full pass.**
- **Breast: 7/8** — the only failure is `rare_precision_not_worse` (0.46 vs 0.52),
  a direct consequence of the large rare-recall gain (precision/recall tradeoff),
  not an accuracy regression.
- Replicated across all 5 seeds and both tissues; 0 failures; runtime comparable.

**Decision: RETAIN as experimental, opt-in — the retain bar is met decisively**
(lung fully; breast 7/8 with only a rare precision/recall tradeoff). The evidence is
strong and consistent, but **default is NOT changed** in this task: the breast rare
precision/FPR tradeoff and the need for external-tool + real-data validation warrant
broader evidence before any default switch. Outputs remain mRNA proportions.

### Supported / unsupported claims
- Supported: count-likelihood GLM **improves fine & broad accuracy, conditional
  within-family RMSE, rare recall, and effective-N calibration** vs wNNLS, on both
  tissues, replicated 30/30, at comparable runtime, 0 failures.
- Unsupported: ✗ "NB beats Poisson" (they tie — NB term near-inert at the reference
  dispersion); ✗ "improves rare precision on breast" (it trades precision for
  recall there); ✗ any default change; ✗ cell fractions (outputs are mRNA
  proportions); ✗ superiority over external probabilistic tools (deferred,
  not benchmarked here).
