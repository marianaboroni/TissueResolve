# Distribution-aware pilot on collinear T/NK states (breast)

Pilot probing whether distribution-aware reference modeling (low-rank covariance /
donor×state metacells) can break the collinear-T/NK identifiability ceiling that mean
profiles (BC > 0.98) hit. Script: `benchmarks/dev/distribution_aware_tnk_pilot.py`
(experimental; outputs gitignored; no defaults changed).

## Probe A — headroom: are collinear states separable beyond the mean? (cell level, donor-held-out AUROC)

| collinear pair | BC (mean) | mean-only centroid | LDA (shared cov) | QDA low-rank cov | logistic |
|---|---|---|---|---|---|
| CD4 helper ∥ CD8 mem | 0.993 | 0.972 | 0.981 | 0.970 | 0.980 |
| CD4 αβ ∥ Treg | 0.989 | 0.930 | 0.967 | 0.935 | 0.953 |
| CD8 mem ∥ eff-mem CD8 | 0.988 | 0.931 | 0.985 | 0.978 | 0.980 |
| CD4 helper ∥ Treg | 0.986 | 0.935 | 0.958 | 0.925 | 0.956 |
| CD4 helper ∥ CD4 αβ | 0.986 | 0.869 | 0.946 | 0.938 | 0.952 |

**Key finding: the collinear states ARE separable at the cell level (AUROC ≈ 0.93–0.98
on held-out donors), despite mean BC ≈ 0.99.** The mean-CPM Bhattacharyya coefficient
**overstates** collinearity — it is dominated by shared pan-T-cell genes, masking the
real discriminating markers (CD4/CD8 etc.). So the ceiling is **not** absolute
biological identity: separating signal exists.

Nuance: the **mean-direction (NearestCentroid) already captures most of it** (0.87–0.97);
covariance/distribution-aware models add only a modest boost (largest on the hardest
pair CD4 helper ∥ CD4 αβ: 0.869 → 0.95, +0.08). So low-rank covariance helps a little
at the cell level, but the mean direction in a discriminative feature space is the bulk
of the signal. *(Caveat: HVGs were selected on all T/NK cells (unsupervised); PCA was
fit train-only and folds are donor-held-out, so the AUROC leak is minor.)*

## Probe B — does it translate to bulk? metacell vs mean NNLS (T/NK-restricted pseudobulk)

| reference | conditional RMSE ↓ | Pearson to truth |
|---|---|---|
| mean profile (1/state) | **0.082** | **0.832** |
| donor×state metacells (29 cols / 6 states) | 0.106 (worse) | 0.735 (worse) |

**Metacells did NOT help — they hurt.** Expanding the dictionary with donor metacells
adds many near-collinear columns, making the NNLS more ill-conditioned and the
aggregated split noisier.

## Interpretation — the reframed problem

1. **Headroom exists, but at the cell level, not the bulk level.** Individual cells of
   collinear T/NK states are ~95% separable; yet bulk is a *sum*, so the few
   discriminating genes are diluted and the **mean profiles remain collinear** — which
   is exactly why bulk deconvolution abstains/struggles on these states.
2. **Metacells are the wrong lever here** (negative: worsens conditioning).
3. **Low-rank covariance adds only a modest cell-level edge** over the mean direction;
   on its own it is unlikely to move bulk recovery much (the summation dilution remains).
4. **The promising lever is discriminative-feature-driven deconvolution:** the genes
   that give the high cell-level AUROC (within-family discriminative markers) must
   *drive* the bulk fit (e.g., a within-family discriminative-marker GLS / reweighting),
   rather than the full mean profile where they are diluted. This is a different,
   targeted direction — not metacells, not raw covariance.

## Conclusion / recommendation

- **Positive:** collinear T/NK states are genuinely distinguishable from data; the
  mean-CPM BC overstates the ceiling. There is signal to exploit.
- **Negative (this pilot):** donor×state metacells do **not** help bulk recovery
  (worse RMSE/Pearson); low-rank covariance gives only a modest cell-level gain.
- **Next step tested below:** a within-family **discriminative-marker reweighted**
  bulk solve for collinear families.

## Follow-up — discriminative-marker reweighting (tested; negative)

Implemented `experimental/discriminative_within_family.py`: a within-family
conditional split (weighted NNLS) where per-gene weights ∝ between-sibling variance
(power sweep), so the genes that separate states drive the split. Tested on the same
T/NK pseudobulk (8 states, 800 HVG, 5 seeds);
`benchmarks/dev/discriminative_within_family_benchmark.py`.

| weight power | conditional within-T/NK RMSE ↓ | Pearson |
|---|---|---|
| **0.0 (uniform baseline)** | **0.080** | **0.654** |
| 1.0 | 0.098 | 0.523 |
| 2.0 | 0.118 | 0.414 |
| 4.0 | 0.133 | 0.346 |

**Reweighting toward discriminative genes does NOT help — it monotonically hurts.**
Why: the discriminating genes are few, noisy, and **donor-variable**; upweighting them
in bulk amplifies donor shift / sampling noise instead of signal, while the shared
genes that stably anchor the NNLS are down-weighted. The per-cell separability does
**not** transfer to a bulk reference reweighting.

## Overall conclusion (distribution-aware pilot)

| lever | collinear T/NK recovery | verdict |
|---|---|---|
| mean profile (baseline) | RMSE 0.080 / Pearson 0.654 (uniform NNLS) | reference point |
| low-rank covariance | modest cell-level AUROC gain only | does not move bulk |
| donor×state metacells | RMSE 0.106 (worse) | negative |
| discriminative reweighting | RMSE 0.098–0.133 (worse) | negative |

- **Headroom is real but cell-level:** collinear T/NK states are ~95% separable per
  cell; the mean-CPM BC overstates the ceiling.
- **None of the bulk-reference levers exploit it:** summation dilutes the discriminating
  genes and donor variability makes them unreliable, so metacells (ill-conditioning),
  covariance (modest), and discriminative reweighting (amplifies noise) all fail to
  improve — or worsen — the bulk collinear split.
- **Honest implication:** in *bulk*, collinear within-family states remain
  fundamentally limited; the right response stays **abstention / family-level
  reporting** (unresolved mass, adaptive resolution). Exploiting the real cell-level
  separability requires data that preserves it — single-cell / single-nucleus or
  spatial context — not a reweighted bulk mean reference. Default behaviour unchanged;
  all of this remains experimental and is retained as a documented negative result.

## Update — supervised conditional estimator (positive)

The mean-geometry levers above all failed, but a **supervised, donor-held-out
conditional family estimator** (different hypothesis) DOES partially recover collinear
families (breast T/NK +14.8%, Myeloid +32%, several lung families +30-58%) without
raising rare FPR, with a reject option for non-learnable families. See
`docs/dev/CONDITIONAL_FAMILY_ESTIMATOR_RESULTS.md` (verdict: promising_experimental).
