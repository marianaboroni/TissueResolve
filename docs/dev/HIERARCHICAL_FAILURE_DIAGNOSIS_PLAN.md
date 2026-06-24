# Hierarchical-Failure Diagnosis Plan (Stage 0 → Phase 1)

Methodology for the diagnosis executed in `HIERARCHICAL_FAILURE_DIAGNOSIS.md`.
Diagnosis only — no algorithm changes. Baseline frozen first (Stage 0).

## Stage 0 — baseline freeze (done)
`benchmarks/baselines/current_baseline/`: `manifest.json` (git commit, package +
external-tool versions, seeds, donor splits, dataset versions, default method
params), `config.yaml` (frozen benchmark configuration), `method_versions.tsv`,
`summary_metrics.tsv` (frozen baseline headline CIs). These are the immutable
comparison point; v2 must beat them on the promotion gates before any default
change (rules 3, 4).

## Phase 1 — four oracle experiments
Goal: localise the hierarchical failure to (broad stage | broad→fine coupling |
fine signatures/panel | gating | spatial smoothing) before committing to any of
the expensive proposed methods.

| Exp | Question | Manipulation | Decision rule |
|---|---|---|---|
| E1 oracle broad mass | Is fine failure propagated broad error? | inject TRUE family mass into `combine_family_and_conditional_estimates` | if it restores accuracy → broad stage/coupling is the cause |
| E2 oracle fine signatures | Is within-family gene selection insufficient? | family panels = reference within-family DE vs global vs default | if oracle panel ≫ default → panel is the cause |
| E3 pre- vs post-gating | Does gating over-abstain / over-sparsify? | all-subtype combine vs gated/unresolved vs `allow_unresolved=False` | if ungated ≫ gated → gating is the cause |
| E4 spatial smoothing | Does the CAR term over-smooth fine? | λ_spatial ∈ {0, 0.02, 0.1, 0.5} | if weak-λ ≥ default on accuracy + calibration → default λ too strong |

Metrics: fine Pearson/RMSE, within-family conditional Pearson, effective-N /
complexity error, unresolved mass fraction, n unresolved families (bulk);
fine Pearson, oversmoothing, local RMSE, Moran's-I MAE, domain ARI (spatial).

Data: held-out-donor synthetic (broad+fine mRNA truth); donor-disjoint;
oracle panels from reference donors only (rule 11). Deterministic seeds.

## Exit criterion
Rank the dominant failure modes and map each proposed v2 change to "justified /
not justified" by the evidence. **Stop. Report. Await sign-off before Phase 2.**
(Outcome: gating = primary; broad coupling = secondary; gene panel = not a
bottleneck; spatial = default λ too strong. See the diagnosis doc.)
