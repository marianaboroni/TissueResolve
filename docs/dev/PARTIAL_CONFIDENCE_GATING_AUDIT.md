# Partial-Confidence Gating Audit (Phase 2A, step 2)

Audit of the current hierarchical gating **before** any code change. Baseline
verified frozen: `git diff src/tissueresolve/ d15a3da..HEAD` is **empty** (the new
commit `9133194` only added Phase-1 docs/baseline/diagnostics), so the frozen
metrics remain valid. Oracle table md5 unchanged; `hier_standard` fine = 0.128.

## Where gating lives
`src/tissueresolve/reference/hierarchy.py`:
- `assemble_hierarchical_estimates()` (L921) — orchestrates gating + mass allocation.
- `evaluate_within_family_resolvability()` (L675) — **family-level** gate.
- `compute_within_family_subtype_confidence()` (L795) — **per-subtype** confidence.
- `estimate_partial_subtype_resolution()` (L860) — applies the subtype gate (default path).
- `decide_unresolved_families()` (L883) / `add_unresolved_family_mass()` (L613) — family-level path.
- `combine_family_and_conditional_estimates()` (L~595) — `family × conditional` → absolute fine.

## Thresholds (current defaults)
| Threshold | Default | Role |
|---|---|---|
| `unresolved_threshold` | 0.10 | family resolvable if mean within-family separability (1−BC) ≥ this |
| `min_discriminating_genes` | 10 | family resolvable if worst pair has ≥ this many \|log2FC\|>1 genes |
| `within_family_spillover_threshold` | 0.30 | family resolvable if mean sibling Pearson ≤ this |
| `subtype_confidence_threshold` | 0.10 | **subtype kept iff confidence ≥ this** (the binary subtype gate) |

## Two code paths
1. **Default** (`allow_partial_resolution=True AND allow_unresolved=True`, L969):
   per-subtype "partial" branch via `estimate_partial_subtype_resolution`.
2. **Fallback** (`allow_unresolved=False` or partial off, L992): family-level
   `decide_unresolved_families` + `add_unresolved_family_mass` (zeros whole
   families).

## The core defect — confidence is computed continuous, then collapsed to BINARY
`compute_within_family_subtype_confidence` returns a **continuous** score per
subtype (`confidence = (1 − BC) to closest sibling`, halved if discriminating
genes < `min_discriminating_genes`). But `estimate_partial_subtype_resolution`
uses it through a **hard threshold**:
```
confident = [m for m in members if confidence[m] >= subtype_confidence_threshold]  # 0.10
sub[m]   = family_mass * conditional[m]        for m in confident   # FULL mass
residual = family_mass * (1 - conditional[confident].sum())          # everything else → unresolved
```
So a subtype is **all-or-nothing**: confidence ≥ 0.10 → keeps its *entire*
conditional share; confidence < 0.10 → **zeroed**, its share dumped to
`unresolved_<family>`. The continuous evidence is discarded.

## Why this produces the measured collapse
The collinear families (T/NK, Myeloid, Epithelial, Endothelial, Mural) have
within-family separability (1−BC) ≈ 0.03–0.12 — at or below the 0.10 threshold —
so most of their subtypes fall under it and are zeroed. Phase-1 measured the
consequence: **70–86% of mass parked unresolved, effective-N 10.2 → 2.5, fine
Pearson 0.128** (vs 0.609 ungated).

## Confidence inputs today
- **Continuous? Internally yes; used as binary** (the defect).
- **Sample-specific? No** — confidence is reference-derived only (separability +
  discriminating-gene count), identical across samples/spots.
- **Missing signals:** bootstrap stability, gene-masking stability, query gene
  overlap, panel-perturbation stability, cross-solver agreement, optimisation
  uncertainty — none feed confidence today.

## Mass handling
- Allocation is **mass-conserving by construction**: for each family,
  `resolved_f + unresolved_f = family_mass_f`, and `Σ_f family_mass_f = 1`, so the
  combined row sums to 1 *before* any normalisation. No post-hoc renormalisation
  alters richness.
- `combined = concat(resolved_fine, unresolved_mass)` — unresolved columns are
  explicit, never silently dropped.

## Report / display
Phase-1 confirmed the bulk report applies **no top-N population filter by
default** (only warnings/table-row caps). The richness collapse is a *data-level
gating* effect, not display filtering.

## Implication for the Phase-2A fix
Keep the exact mass-conservation arithmetic; **replace the binary subtype
membership with a continuous, calibrated weight** `c_{k,s} ∈ [0,1]`:
```
resolved_theta[k,s] = (family_mass[f,s] · conditional[k,s]) · c_{k,s}
unresolved[f,s]     = family_mass[f,s] − Σ_{k∈f} resolved_theta[k,s]   (≥ 0, since c≤1, Σ conditional =1)
```
This is a strict generalisation: `c∈{0,1}` reproduces today's hard gate; the
calibrated `c` recovers the mass the hard threshold over-discards. It reuses the
existing confidence signal and adds richer, *sample-aware* features (step 4),
calibrated on held-out mixtures (step 5). No gene-panel or solver change.
