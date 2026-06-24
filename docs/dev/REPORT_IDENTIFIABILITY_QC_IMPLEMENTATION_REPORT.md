# QC-First Report — Implementation Report (Task A)

**Status: complete + validated.** Self-contained builder
`benchmarks/build_qc_first_report.py` (does NOT touch stable report code; rule "no
default change"). Outputs under `benchmarks/outputs/qc_report/`: `report.html`,
`technical_appendix.html`, `figures/figure_manifest.tsv`, `source_data/` (9 TSVs),
`warnings.json`, `methods.txt`, `figures/`. Nothing committed.

## Evidence chain (QC before predictions)
reference support → donor coverage → 3-level signature identifiability → query
compatibility → **full-panel deconvolution reliability** → trusted resolution →
final predictions. 10 main sections, 10 main figures (≤ the 15–22 target; only
decision-relevant figures included per the central rule).

## Central inclusion rule — enforced
A figure is in the main report only if it helps decide reliability or interpret a
prediction. Heavy material (full 91-pair identifiability, mixture-recovery per pair,
donor-stability per pair, TCGA per-sample audit, raw multisplit/multiseed tables,
soft-gating comparison, joint-hierarchy negative result) is referenced in
`technical_appendix.html`, not embedded as main figures.

## Three identifiability levels kept DISTINCT (not merged)
The report shows, per family, three separate bars (figure `figI1_three_levels`):
1. cell-level donor-held-out AUROC (~1.0 — labelled "NOT deconv proof");
2. pairwise mixture recovery (1−MAE, ~1.0);
3. **full-panel conditional reliability (1−RMSE)** — the basis for trust.
A banner repeats: "cell-level separability does not guarantee accurate
deconvolution in a complex mixture."

## Trusted resolution — derived from full-panel metrics (NOT cell AUROC)
`trusted_resolution()` combines shared-lineage fraction, condition number, and the
per-family full-panel conditional RMSE (computed from one held-out seed not used
for tuning). Result (metric-derived, `source_data/trusted_resolution.tsv`):

| family | shared-lineage | trusted resolution |
|---|---|---|
| B/Plasma | 0.86 | **selected_fine** |
| Epithelial | 0.94 | broad_only |
| Myeloid | 0.96 | broad_only |
| Endothelial | 0.99 | broad_only |
| T/NK | 0.99 | broad_only |
| Mural | 1.00 | broad_only |

**No family is full-fine.** This is *more conservative* than the provisional
classification: Epithelial fell to broad-only because its full-panel conditional
RMSE exceeded the cross-family median — a derived, not hardcoded, decision (as the
brief requires). Fine predictions in §6 are **gated**: only selected-fine families
may be shown at subtype level; broad-only families are shown at family level with
unresolved mass.

## Status cards (§1) — metric-derived
reference suitability (donor range/family), donor coverage (41 donors), **fine
full-panel reliability (median conditional RMSE)**, query compatibility (% markers
detectable), recommended resolution (n selected-fine / n broad-only), soft-gating
status (experimental, no 2nd tissue). Each card carries status · metric ·
explanation · affected result (no bare PASS/FAIL).

## Predictions
- **Bulk (§6):** TCGA-TNBC broad-family composition (aggregated, clustered),
  labelled "RNA-derived proportions, NOT cell fractions"; fine gated out (no
  full-fine family).
- **Spatial (§7):** real-Visium cross-method concordance embedded with explicit
  "concordance is NOT accuracy / no ground truth" + "Moran's I is not accuracy".
- **Benchmark (§8):** bulk (held-out CI) and spatial (synthetic CI) **kept
  separate**.

## Tests (`tests/test_qc_report.py`, 9 pass)
trusted-resolution uses full-panel not cell-AUROC; no family marked full-fine;
status cards carry metric+threshold+affected; QC sections precede predictions;
3 identifiability levels distinguished; fine gated + no spatial-accuracy claim;
every main figure has caption + source data; donor support present; technical
material in appendix not main.

## Limitations
Single tissue; per-family conditional RMSE from one held-out seed (directional);
trusted-resolution thresholds (shared<0.95, conditional≤median) are defensible but
sensitivity-worthy; builder is self-contained (not yet merged into the production
report module — a deliberate isolation pending the second tissue).

## Gate
Task A complete + validated. **Per the strict order, stopping before Task C**
(donor-stable contrastive signatures) as instructed ("do not automatically
continue"). Task B (second-tissue candidates) delivered in parallel.
