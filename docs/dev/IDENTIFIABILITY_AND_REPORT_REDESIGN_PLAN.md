# Identifiability & Report-Redesign Plan (Stage 0)

Written before any algorithm/report change. Baseline frozen + verified
(`git diff src/tissueresolve/ d15a3da..HEAD` empty). 19 Phase-2 tests + 119
benchmark tests pass. sklearn 1.6.1 available. Nothing committed.

## Current methodological conclusions (from Phases 1–2)
- Hard confidence gating was the primary hierarchical failure; **partial
  confidence-weighted gating (soft gate)** fixes it (validated, 9/9 Stage-1 gates;
  experimental — 2nd tissue pending).
- **Soft joint hierarchy is a negative result** (≡ sequential; broad–fine coupling
  already consistent) — not pursued.
- Residual error localised to: (a) broad-signature accuracy, (b) **within-family
  conditional discrimination** (dominant), (c) weak confidence predictors,
  (d) spatial over-smoothing (default λ too strong).

## Hypotheses to test (this stage)
- **H1 (Stage 1, decisive):** the requested fine subtypes are *not* reliably
  distinguishable from donor-stable, query-detectable genes — i.e. there is a low
  **identifiability ceiling** that no solver/contrast can exceed.
- H2 (Stage 2): a donor/protocol/query-aware broad panel can reduce broad RMSE
  toward oracle.
- H3 (Stage 3): subtype-contrast signatures beat conventional panels — **only
  pursued if Stage 1 shows headroom** (rule 17).
- H4 (Stage 4): pipeline-derived features improve soft-gating confidence.
- H5 (Stage 5): level-specific spatial λ (λ_broad>λ_fine) reduces over-smoothing.

## Stop/go gates
- **Stage 1 (gate):** if most within-family pairs have donor-held-out AUROC < 0.60
  on query-detectable genes → the resolution is **not identifiable**; do NOT build
  a more complex solver/contrast (skip Stage 3). Preserve unresolved mass; report
  broad-only; drive the report's "trusted resolution" from these results.
- Stage 2/3/4/5 each promote only if their predefined margins pass (see each
  stage); otherwise the component stays experimental and the simpler path wins.
- Stage 6: no promotion of soft gating until a real second tissue validates.

## Modules changed vs untouched
- **Changed (experimental, isolated):** new `experimental/contrastive_signatures/`
  *only if Stage 1 warrants*; new diagnostics scripts under `benchmarks/diagnostics/`.
- **Untouched (stable core):** bulk/spatial solvers, default gating, gene-panel
  selection, `api.py`. Soft gating stays experimental. No PyTorch/JAX/OT/NB-CAR-VI.

## Datasets & splits
Breast atlas (41 donors, raw counts). **Donor-held-out** classification (GroupKFold
by donor); query-detectable = genes detected in TCGA-TNBC bulk. No test-donor use
for gene selection (rules 10–12). Second tissue: none adequate locally (Stage 6 is
a selection plan).

## Metrics
Identifiability: donor-held-out AUROC/AUPRC/balanced-acc, donor variance, stable &
query-detectable gene counts, signature cosine/Pearson, condition number.
Broad: held-out broad RMSE by panel, conditioning. Contrastive: conditional RMSE,
spillover, richness. Confidence: AUROC, Brier, ECE. Spatial: fine/broad RMSE,
boundary F1, oversmoothing, Moran's-I preservation, local RMSE.

## Report redesign (Stage 7) — central rule
**A figure appears in the MAIN report only if it helps decide whether a prediction
is reliable, or helps interpret the final prediction.** Everything else → technical
appendix / source data. QC-first: reference quality → signature/identifiability →
input compatibility → prediction quality/trusted resolution → predictions →
benchmark → warnings → methods. ≤10 sections, ~15–22 essential figures. Fine
predictions **gated by Stage-1 trusted resolution**. Bulk/spatial benchmarks
separate. No spatial-accuracy claim without truth; RNA-proportion (not cell
fraction) labelling.

### Main-report figures (planned): exec status cards; cells&donors per broad &
fine; donor-coverage heatmap; reference-suitability components; pairwise
identifiability heatmap; top confusable pairs; query-detectable marker support;
recommended-resolution-by-family; gene overlap by level; prediction-quality
dashboard; unresolved-mass by family; effective-N/entropy; trusted-resolution
summary; final broad composition (bulk); trusted fine; spatial dominant-broad map;
unresolved-mass map; bulk & spatial benchmark summaries.
### Appendix: full separability/spillover matrices, all pairwise tables, large
signature heatmaps, spot pie charts, raw benchmark tables, composite-score detail,
install logs, experimental diagnostics.

## Test plan
donor-held-out split integrity; no leakage; identifiability metric correctness;
query-detectable filtering; (contrast centering if built); level-specific λ
behaviour; soft-gating default unchanged; failed joint hierarchy not promoted;
report ordering (QC before predictions), trusted-resolution gating, separate
benchmark summaries, source data per figure, deterministic palette.

## Risks
Single tissue limits generality; identifiability ceiling may itself be panel-
sensitive (mitigated by panel sensitivity analysis); report redesign is large
(staged: Stage-1 results first drive the trusted-resolution layer).

## This-turn scope (evidence-gated)
Execute **Stage 1 (identifiability ceiling)** — the decisive gate — plus Stage-6
selection plan and the report-redesign spec. Sequence Stages 2/3/4/5/7 after the
Stage-1 evidence (Stage 3 likely *not* warranted if the ceiling is low). Stop and
report at the Stage-1 gate per the prompt's discipline.
