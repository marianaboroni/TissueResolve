# Subpopulation Identifiability Report (Stage 1)

Diagnosis only — no algorithm changed; nothing committed. Estimates the ceiling of
within-family subtype discriminability **before** building any new solver/contrast,
using three complementary probes. Donor-held-out throughout; query-detectable =
genes seen in TCGA-TNBC bulk (4971/5000); simple classifiers only (rule: no deep
models). Scripts: `benchmarks/diagnostics/identifiability_ceiling.py`,
`mixture_recovery_ceiling.py`.

## The three probes give an apparently contradictory — then coherent — picture

| Probe | Question | Result (all families) |
|---|---|---|
| **(1) cell classification** | Can a classifier tell single cells of subtype A from B? | **AUROC ≈ 1.0** (T/NK 0.998, Myeloid 0.997, Endothelial 1.0, …) |
| **(2) pairwise mixture recovery** | Can NNLS recover A:B ratio from a 2-signature mixture? | **MAE ≈ 0.003–0.005** (trivial, even at signature Pearson 0.98) |
| **(3) full deconvolution (Phase 2)** | Within-family conditional error in the real benchmark | **RMSE ≈ 0.31** (hard) |

Probes (1) and (2) say subtypes are **easily** distinguishable; probe (3) says the
real within-family deconvolution is **hard**. Reconciling them is the key insight.

## Resolution of the paradox
- Subtypes are **biologically real and cell-level separable** (AUROC ≈ 1.0).
- Their **signatures are highly collinear**: subtype signature Pearson **0.84–0.98**,
  and the **shared-lineage fraction ‖B‖/‖S‖ = 0.96–0.99** — i.e. the subtype-specific
  contrast `D = S − B` is a **tiny residual (1–4%) on a dominant shared family
  program**.
- In *isolation* (2 signatures, matched reference, mild noise) NNLS still exploits
  that small contrast and recovers ratios well (probe 2).
- But the **real** within-family error (0.31) arises from the combination of
  **(a) full multi-subtype ill-conditioning** (≈30 collinear signatures at once, far
  worse than pairwise cond 3–10) and **(b) reference↔query donor mismatch**, which
  makes the small contrast `D` **unreliable** across donors. Phase-1 already showed
  an oracle within-family DE panel did *not* help — consistent with the limit being
  the **donor-stability of the contrast**, not its selection.

## Stop/go gate
The gate "if most pairs have donor-held-out AUROC < 0.60 on query-detectable genes →
do not force with a solver" is **NOT triggered** (0% of 91 pairs below 0.60; the
resolution is *not* fundamentally impossible). **However**, the deconvolution-relevant
evidence (shared-lineage dominance + donor-fragile contrast) means a naive
within-family solver/panel will *not* help — the principled (and only untested)
lever is **donor-stable contrast modules**.

## Recommended resolution per family — IMPORTANT caveat
The cell-classification and pairwise-mixture probes would label *every* family
"fine-resolvable", but **that is the wrong basis for a deconvolution report**: the
real benchmark shows fine resolution is fragile (conditional RMSE 0.31, soft-gate
fine Pearson ~0.62 vs flat ~0.72). The report's **trusted resolution must be driven
by the practical deconvolution reliability** (conditional error + soft-gating
confidence + unresolved mass), NOT by cell-classifiability. Operationally: families
with high shared-lineage fraction (≥0.95: T/NK, Mural, Myeloid, Endothelial) should
be presented as **broad-confident, fine-cautious** (soft-gated, unresolved mass
shown); B/Plasma and Epithelial (lower collinearity) are the most fine-resolvable.

## Outputs
`subtype_identifiability_pairs.tsv` (91 pairs: AUROC/AUPRC/bacc by panel×classifier,
signature similarity, condition number), `subtype_identifiability_family_summary.tsv`,
`subtype_query_detectability.tsv`, `subtype_donor_stability.tsv`,
`mixture_recovery_ceiling_pairs.tsv`, `mixture_recovery_family_summary.tsv`. Figures
(`figures/identifiability/`): AUROC distribution, ref-vs-query, recommended
resolution, similarity-vs-AUROC, mixture-recovery-vs-collinearity (each with
`.data.tsv` + caption).

## Implications for the remaining stages
- **Stage 3 (contrastive signatures):** *conditionally* justified — the shared
  program dominates (96–99%), so explicitly modelling `S=B+D` and using **donor-
  stable contrast modules** is the principled attempt. But Phase-1 (oracle DE panel
  no help) tempers expectations: if the contrast is donor-unstable, contrast
  selection alone won't fix it. **Recommend a focused, gated Stage-3 attempt; do not
  assume success.**
- **Stage 2 (broad panel):** broad families are well-recovered pairwise; the broad
  RMSE gap vs oracle is multi-family conditioning — a panel audit is reasonable but
  lower priority than honest reporting.
- **Report (Stage 7):** the central deliverable — communicate that fine subtypes are
  *real but deconvolution-fragile*; gate fine predictions by practical reliability.

## Limitations
Single tissue; cell-classification AUROC is a *necessary-not-sufficient* ceiling;
the pairwise mixture recovery omits donor mismatch (deliberately, to isolate
collinearity) so it is optimistic; the donor-mismatch contribution is inferred from
the Phase-2 benchmark rather than measured in a single combined experiment.
