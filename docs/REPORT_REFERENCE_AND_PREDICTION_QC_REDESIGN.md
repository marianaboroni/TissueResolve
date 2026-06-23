# Report Redesign — Reference & Prediction QC (Stage 7 blueprint)

Specification for the QC-first report redesign. **Blueprint only this turn** — the
full `report.html` rebuild is a large front-end effort sequenced after the Stage-1
evidence (which supplies the trusted-resolution layer) and is implemented next.
Nothing committed.

## Central inclusion rule (the most important rule)
A figure appears in the **main report** only if it helps the user **decide whether
a prediction is reliable** or **interpret the final prediction**. Everything else →
`technical_appendix.html` / `source_data/`.

## Section order (≤10; QC BEFORE predictions)
1. **Executive decision summary** — status cards (each: metric, threshold,
   explanation, affected result; no bare PASS/FAIL).
2. **Reference quality** — cells & donors per broad/fine label (distinguish "many
   cells, one donor" from "many donors"); donor-coverage heatmap; suitability
   components; imbalance; low-support populations; batch/donor confounding.
3. **Signature & identifiability quality** — pairwise identifiability heatmap; top
   confusable pairs; query-detectable marker support; **recommended resolution per
   family driven by deconvolution reliability (shared-lineage fraction + conditional
   error + soft-gating confidence), NOT cell-classifiability** (per Stage-1 caveat);
   signature conditioning.
4. **Input compatibility** — reference/query gene overlap by level; signature-gene
   detection by family; count/normalisation status; library-size; protocol/sc-vs-sn
   mismatch; (spatial) H&E/spot alignment.
5. **Prediction quality & trusted resolution** (precedes predictions) —
   reconstruction/residual QC; solver mode; soft-gating confidence; unresolved mass
   by family; effective-N/entropy/dominant; trusted-resolution summary; fine-caution
   table.
6. **Final bulk predictions** — broad composition (stacked bar + heatmap, clustered);
   **trusted fine only** (gated by §5); unresolved gray; RNA-proportion labelling.
7. **Final spatial predictions** — H&E + dominant-broad map; top broad maps; trusted
   fine maps; entropy/dominant/unresolved maps; Moran's-I ranked; boundary/over-
   smoothing QC where truth/histology exists. Broad before fine; no accuracy claim
   without truth.
8. **Benchmark evidence** — bulk and spatial **kept separate**; show *why* a method
   wins, not just rank.
9. **Warnings & limitations.**
10. **Methods & source data.**

## Trusted-resolution gating (from Stage 1 + soft gating)
Per family, classify `broad_only | selected_fine | fine` from: shared-lineage
fraction, full-panel conditional error, soft-gating confidence, unresolved mass.
**Families flagged broad_only must NOT be shown as confident fine subtypes in the
main report** (fine figures gated; unresolved mass displayed instead).

## Deterministic hierarchical palette
Broad families → distinct hues; subtypes → subtones of the parent hue; unresolved →
gray. Consistent across reference QC / signature QC / bulk / spatial / benchmark /
unresolved.

## Per-figure metadata (all main figures)
title · one-sentence purpose · full caption · axes def · color def · data/method
source · "what to check" · limitation · `source_data/` link. Save PNG + SVG/PDF +
`.data.tsv`. Plotly where interactive helps.

## Appendix (moved out of main)
full separability/spillover matrices · complete pairwise tables · large signature
heatmaps · all spot pie charts · raw benchmark tables · composite-score detail ·
install logs · model diagnostics · experimental figures.

## Targets
main ≤10 sections, ~15–22 essential figures, no raw tables/debug dumps; appendix
holds the rest. Reuse existing `src/tissueresolve/report/` machinery; add a
QC-first ordering layer + the trusted-resolution gate; do not change stable
prediction code.

## Implementation status
Blueprint complete; build sequenced next (after Stage-1 trusted-resolution layer is
wired). Existing report tables that already feed this: reference suitability
(`outputs/bulk/reference_suitability_*`), identifiability
(`subtype_identifiability_*`, `mixture_recovery_*`), prediction QC
(`tcga_prediction_complexity_audit`, soft-gating outputs), benchmark
(`benchmark_report.html`).
