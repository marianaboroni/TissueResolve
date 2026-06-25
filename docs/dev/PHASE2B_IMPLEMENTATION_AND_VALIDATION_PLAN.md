# Phase 2B — Implementation & Validation Plan (Stage 0)

Written **before** any algorithm change. Baseline verified frozen: `git diff
src/tissueresolve/ d15a3da..HEAD` is empty (core unchanged); Phase 2A is opt-in
(`tissueresolve.experimental.soft_hierarchy`, not imported by default); the
production default still uses the original stable hierarchical/flat path. Phase 2A
result tables present (6). Nothing committed.

## Recorded gate history (rule 5 — do NOT rewrite failed gates)
- **Phase 2A false-resolution gate: FAILED as written** (soft 0.097 vs hard 0.069).
- **Reason:** the hard gate achieves low false-resolution mainly through excessive
  abstention (68% unresolved), so "soft must beat hard on false-resolution" is an
  apples-to-oranges comparison.
- **Revised PROSPECTIVE criterion (declared here, before Stage-1 results):** soft
  gating's false-resolution must be **non-inferior to the *ungated* hierarchy**
  (margin ≤ +10% relative), **while** false-abstention is materially below hard
  gating. The original failed result is preserved in
  `PARTIAL_CONFIDENCE_GATING_FINAL_REPORT.md` and is not retroactively changed.

## Critical data limitation (affects the replication requirement)
The only single-cell atlas available is **breast cancer**. The sole non-breast
dataset (scanpy PBMC68k_reduced) has **no donor labels, is not raw counts, and is
700 cells** — it cannot support donor-disjoint, count-level validation. **A second
tissue satisfying the requirements is NOT available in this environment.**
Consequence: the **2-tissue replication requirement cannot be met**, so soft
gating (and any Stage-2 method) **cannot be PROMOTED** regardless of metrics.

## Stop/go logic (integrity-preserving reading)
- Stage 1 runs an **expanded single-tissue (breast) validation** (≥10 test seeds,
  more scenarios) — strengthening evidence within one tissue.
- The **prospective gates** (metric thresholds vs ungated/hard) are the operational
  stop/go for *implementing + experimentally benchmarking* Stage 2 (which stays
  experimental, never promoted).
- The **2-tissue replication requirement is recorded as UNMET** → promotion of any
  method is blocked this stage; everything remains experimental (rules 2, 12).
- Stage 3 (spatial level-specific smoothing) is **independent** of Stage 2 and the
  tissue gate; it may be implemented + tested (no default change without replication).

## Independent Phase-2A validation design (Stage 1)
- Splits (donor- AND seed-disjoint; rules 3,4): calibration seeds {0,1}, validation
  {2}, **test {3..12}** (10 seeds). Each seed = an independent donor split.
- Scenarios: balanced, imbalanced, similar_subtypes (partially/non-separable),
  rare, missing_population, reduced_gene_overlap. (cross-platform/batch via the
  assay-aware donor split where feasible.)
- Modes: flat / hard_gate / ungated / soft_gate.
- Metric panel: broad (Pearson/Spearman/RMSE/MAE/JSD/Aitchison); fine (absolute +
  conditional Pearson/RMSE, macro error, rare sensitivity, false-positive
  detection); resolution (unresolved precision/recall, false-resolution,
  false-abstention, correct-resolution-level, resolved/unresolved family counts);
  complexity (effective-N, entropy, richness, Gini, dominant, top-5); safety (mass
  conservation, negatives, numerical failures, runtime).
- Statistics: bootstrap 95% CI + paired tests over seed×scenario replicates.

### Prospective soft-gating gates (declared pre-results)
vs **ungated**: broad RMSE ≤ +2%; fine RMSE ≤ +2%; false-resolution ≤ +10%
relative (revised criterion above); unresolved precision ≥; effective-N & entropy
closer to truth; mass conserved.
vs **hard**: materially lower false-abstention; higher rare sensitivity;
substantially lower effective-N error; lower fine RMSE.
Stability: benefits consistent across the 10 test seeds (CI excludes null).

## Stage 2 — soft joint broad–fine hierarchy (math)
Constrained objective per sample (no strong L1; rule 6):
`min_{π,θ≥0} ‖W_b(y−Bπ)‖² + α‖W_f(y−Sθ)‖² + λ_h‖π−Aθ‖² + λ₂‖θ‖²`,
`Σπ=1, Σθ=1`. Two solvers: **(A) alternating** (NNLS-init θ → π=Aθ → update π|θ →
update θ|π → reconcile → iterate) and **(B) joint** (scipy SLSQP/`minimize` with
simplex constraints). Prefer the simpler if equivalent (rule 16). Module:
`src/tissueresolve/experimental/soft_hierarchy/joint_solver.py`. Ablations:
init {flat,sequential,uniform,broad-first}; λ_h {0,weak,mod,strong}; λ₂ {0,weak,mod}.
Error decomposition: broad / conditional-fine / consistency / gating / total.

## Stage 3 — level-specific spatial smoothing (math, no NB-CAR VI; rule 8)
Experimental `λ_broad`, `λ_fine` applied at the two resolution levels of the
existing NB-CAR (no new model). Sweep λ_broad∈{0,.005,.01,.02,.05,.1},
λ_fine∈{0,.0025,.005,.01,.02,.05}; modes none/broad-only/fine-only/equal/level-
specific; ≥10 seeds; patterns gradient/sharp/rare-niche/mixed. Metrics: broad/fine
RMSE, JSD, boundary preservation, edge blurring, Moran's-I preservation, domain
recovery, local RMSE, rare-niche detection, runtime.

## Tests & deliverables
Tests: split-leakage, deterministic splits, prospective-gate evaluation, recorded-
failed-gate; joint-solver invariants (non-neg, Σπ=Σθ=1, aggregation consistency,
mass conservation, convergence, λ_h=0 behavior, strong-λ_h ↑ consistency, no
sparsity collapse, backward-compat); spatial (separate λ, both-zero equivalence,
broad/fine-only, deterministic graph, boundary metrics).
Deliverables: the four reports + the six benchmark TSVs named in the brief.

## Stop/go summary
1. Stage 1 prospective gates **fail** → keep soft gating experimental, **stop**.
2. Stage 1 prospective gates **pass** (breast) → implement + experimentally
   validate Stage 2; **promotion blocked** (replication unmet); proceed to Stage 3
   independently.
