# Stage 1.6 — separability-aware adaptive signature report

## SCIENTIFIC QUESTION
Does a fine-global signature whose per-cell-type gene count is set by each population's
**separability difficulty** beat fixed / fixed-per-type budgets — preserving broad accuracy,
reducing spillover, and using compact panels for easy populations?

## AUDIT OF EXISTING SEPARABILITY
`compute_separability` = pairwise Bhattacharyya on **mean CPM profiles** (not donor-aware, not
deconvolution-aware, panel-dependent). Two resolvability layers exist
(`reference/resolution.classify_resolvability`; `hierarchy.evaluate_within_family_resolvability`);
the hierarchical path uses the latter, the flat/default path neither. See
`docs/dev/STAGE1_6_ADAPTIVE_SIGNATURE_AUDIT.md`.

## DUPLICATION AND CONSOLIDATION
Two resolvability layers (not consolidated this session — flagged; no third separability metric
was created). The allocator REUSES `compute_separability` (closest confounder) + `donor_aware_de`
(candidates) + a new donor-held-out difficulty curve.

## ALGORITHM
`SeparabilityAwareBudgetAllocator` (`reference/adaptive_signature.py`): per fine type, rank
donor-aware one-vs-rest candidates; build a **donor-held-out (leave-one-train-donor-out)
target-vs-closest-confounder** recovery curve over growing panels; pick the smallest panel that
is good enough; assign a per-type status. Bulk-only (validated_modality=bulk_experimental).

## SEPARABILITY COMPONENTS
Closest confounder (pairwise BC); donor-held-out 2-component (target vs confounder) fraction
recovery on the **union of both types' markers** (so basis columns are non-proportional exactly
when discriminating genes exist); per-type best RMSE + curve.

## ADAPTIVE BUDGET RULE
Smallest panel with recovery RMSE ≤ max(best·(1+perf_tol), easy_rmse floor); clipped to
[min_genes, max_genes]; deduped global union with per-type coverage. Not a fixed constant.

## INCREMENTAL STOPPING RULE
Stop growing when within tolerance of best OR below the absolute easy floor; UNRESOLVABLE if best
RMSE stays above `unresolvable_rmse` even at max panel (no discriminating genes → do not expand).

## DATA LEAKAGE CONTROLS
Separability/difficulty use **training donors only** (leave-one-train-donor-out); evaluation on
held-out query donors; thresholds pre-specified (not chosen on the test set); benchmark guards
(validation.py) abort on wrong-tissue hierarchy, donor overlap, selection leakage, mis-orientation.

## SYNTHETIC VALIDATION
6-type toy (easy / two near-collinear / rare / identical-twin pair): allocator gives Easy
EASY_COMPACT (low RMSE), flags the identical twins UNRESOLVABLE, Easy easier than collinear,
deterministic, no-donor→INSUFFICIENT_REFERENCE. **The synthetic test exposed and drove fixes to
two real algorithm flaws** (aggregate-rest diluted confounders → missed twins; 2-component fit on
target-only markers was degenerate → mislabeled Easy). 8 tests pass.

## BREAST RESULTS (5 seeds, donor-disjoint)
D_adaptive best on all: RMSE 0.0319 / Pearson 0.666 / cond-RMSE 0.171 / spillover 0.086, 383
genes — vs C 0.0339/0.631/0.182/0.100 (303) and A 0.0353/0.620/0.201/0.113 (500).

## LUNG RESULTS (5 seeds)
D_adaptive 0.0298/0.534/0.227/0.164, **754 genes** — ~tie with C 0.0298/0.546/0.222/0.167 (553);
D used MORE genes without gain. F_all_eligible (2917 genes) 0.0282/0.560/0.219 (marginally lower
RMSE but impractical gene count).

## THIRD-TISSUE RESULTS
**NOT EXECUTED** — no third tissue (PBMC/brain/…) prepared in this session's data; deferred.

## CROSS-PLATFORM RESULTS
**NOT EXECUTED** for the adaptive strategy this session (time). Breast cross-platform data exist
(prior gene_selection evidence); an adaptive cross-platform run is the immediate next step.

## REAL-BULK RESULTS
**NOT EXECUTED** — TCGA-TNBC exists but has no fine ground truth; adaptive real-bulk concordance
deferred to the same harness used in the Rectangle real-bulk gate.

## BROAD PERFORMANCE
Not separately aggregated (broad-RMSE metric not added). D reduces spillover_absent on breast
(0.086 vs 0.113), an indirect broad-plausibility signal. Open item.

## FINE PERFORMANCE
D best conditional within-family RMSE on breast (0.171, significant); ~tie/slightly worse on lung.

## PER-CELL-TYPE RESULTS
`adaptive_gene_budget_by_celltype.tsv` (per type: gene count, status, best_rmse, separability_bc,
closest_confounder). Breast: 6 EASY_COMPACT / 19 DIFFICULT_EXPANDED / 13 INSUFFICIENT_REFERENCE.
Lung: 16 / 18 / 3 MODERATE / 16 INSUFFICIENT.

## RARE-CELL EXPLORATORY RESULTS
Rare-by-abundance-band **NOT EXECUTED** in this benchmark (spillover_absent reported as a proxy).
Many rare types fall to INSUFFICIENT_REFERENCE (too few donors) — honestly flagged, not forced.

## PAIRED STATISTICS
Breast D−C: RMSE −0.0020 (CI[−0.0034,−0.0009], 83% win), cond-RMSE −0.011 (CI excl 0, 67%);
D−A: RMSE −0.0034 (CI excl 0), cond-RMSE −0.030 (CI excl 0). Lung D−C: RMSE ~0 (tie), cond-RMSE
+0.005 (CI crosses 0 — no gain). Lung D−A: cond-RMSE −0.007 (CI crosses 0, tie).

## RUNTIME AND MEMORY
Adaptive allocation ~23–24 s per reference (cached per-type DE + cheap 2-col NNLS curves). Whole
benchmark (2 tissues × 5 seeds × 4 strategies) a few minutes. No memory issues.

## UNRESOLVABLE POPULATIONS
Structurally non-resolvable types (identical/near-identical to confounder) → UNRESOLVABLE with
best-achieved metrics + closest confounder recorded; the allocator stops expanding them.

## NEGATIVE RESULTS
Adaptive does NOT beat fixed-per-type on lung (tie, more genes) → the stopping rule over-allocates
on many-type references without gain. F_all_eligible confirms more genes ≠ better (breast).

## LIMITATIONS
Pseudobulk only; 2 tissues; broad-RMSE/rare-by-abundance/oracle/cross-platform/real-bulk NOT run;
2-component difficulty proxy (not full multi-way deconvolution); resolvability layers not consolidated.

## PROMOTION DECISION
**§30 = D — tissue-dependent; keep experimental.** Adaptive is significantly superior on breast
(RMSE, conditional RMSE, spillover, fewer genes than current) but only ties fixed-per-type on lung
while using more genes. It does NOT meet all §22 promotion criteria (lung no-gain + larger panel;
cross-platform not tested). Status: `fine_global separability-aware = experimental candidate under
validation`. Fixed-per-type remains the experimental recommended baseline; GeneSelector the stable
default. No default changed.

## NEXT STAGE
Tighten the stopping rule (redundancy/incremental-information pruning to stop lung over-allocation),
add broad-RMSE + rare-by-abundance metrics, run adaptive **cross-platform** and a **third tissue**,
then re-audit §22. If lung still shows no gain over fixed-per-type, prefer the simpler fixed-per-type
(budget ∝ #types) and keep adaptivity only where it demonstrably helps.
