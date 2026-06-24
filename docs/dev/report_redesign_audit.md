# TissueResolve report redesign — audit

**Scope:** read-only audit of the current generated unified report
(`examples/real_breast_cancer/outputs/report.html`) and the generator
(`examples/real_breast_cancer/scripts/07_generate_reports.py`,
`src/tissueresolve/report/*`, `src/tissueresolve/plotting/palette.py`), written
**before** any redesign code change (per the redesign brief, Part 1).

It records the state of the report as the basis for the redesign. Where a
Stage-0 stabilization pass (warnings population, figure embedding, composite
de-leak, license, tuning quarantine) has already addressed a problem, that is
noted as **[fixed in Stage 0]** so the redesign does not re-do it.

---

## 0. Current structure (as generated)

Section order in `generate_unified_report`:

| # | anchor | title |
|---|--------|-------|
| 1 | `summary` | Executive summary |
| 2 | `reference` | Reference quality |
| 3 | `input` | Input data quality |
| 4 | `bulk` | Bulk deconvolution |
| 5 | `spatial` | Spatial deconvolution |
| 6 | `hierarchical` | Hierarchical broad→fine deconvolution |
| 7 | `resolution` | Resolution, separability & spillover |
| 8 | `benchmark` | Benchmark comparison |
| 9 | `warnings` | Warnings & limitations |
| 10 | `methods` | Methods |
| 11 | `outputs` | Output files & source data |

Report size: ~77 KB single HTML, sticky sidebar nav, component-based
(`metric-card`, `figure-card`, `methods-card`, `warning-card`, …).

---

## 1. Which sections are too complex?

* **Executive summary (1):** a 12-cell `metric_grid` of raw counts plus two long
  prose boxes. It lists numbers but does **not** give a go/no-go decision view
  (no PASS/CAUTION/WARNING/FAIL status cards, no pre-interpretation checklist).
  It is informative but not *decision-oriented*.
* **Resolution / separability / spillover (7):** mixes three distinct concepts
  (pairwise separability, spillover risk, unresolved family mass) in one section
  with large tables; no compact "what resolution can I trust" summary.
* **Benchmark (8):** raw `to_html` tables of method status + composite scores;
  hard to scan; no separation of executed vs imported vs skipped at a glance.

## 2. Which figures are missing, empty, or only linked?

* **[fixed in Stage 0]** Empty figure bodies: every `figure_card` was built with
  `body_html=""` — the prior report had **22 empty `fig-body` divs** (figures
  were *linked* only). Stage 0 now embeds a PNG `<img>` when present, else an
  `<iframe>` to the interactive HTML. Regenerated report: **0 empty bodies**,
  16 PNG + 6 iframe embeds.
* **Still missing (redesign work):** there is no dedicated signature heatmap,
  hierarchy tree/sankey, within-family separability lollipop, marker-support
  plot, reference-suitability component bar, normalization/library-size
  diagnostic, spatial entropy/dominant-fraction/unresolved-mass maps, solver-CV
  plot, or benchmark status/leaderboard plots. The figure set is generic
  (whatever PNG/HTML happens to be in each `figures/` dir) rather than the
  decision-workflow panels the brief specifies.

## 3. Which captions are generic and need rewriting?

* The `_CAPTIONS["_default"]` entry is generic and is applied to **any figure
  not explicitly listed**:
  * subtitle: `"TissueResolve figure"`
  * caption: `"Visual summary of a TissueResolve output. Axes, colors and …"`
* The regenerated report still contains **2 generic-caption hits**
  (`"Visual summary"` / `"TissueResolve output"`). Any figure whose stem is not a
  key in `_CAPTIONS` inherits these. The redesign must (a) give every emitted
  figure a specific caption keyed by stem, and (b) make the `_default` itself at
  least name the figure (title-derived) rather than say "a TissueResolve output".

## 4. Which tables are too large and should be replaced by plots?

* Pairwise separability (497 rows) — shown/linked as a table; should become a
  "top-20 most confusable pairs" lollipop/severity plot, full table as source
  data.
* Spillover risk by cell type — table; should become a high-risk spillover
  network (top edges only).
* Unresolved families / recommended merges — tables; should become an
  "unresolved mass by family" barplot + a compact "recommended interpretation
  level by family" panel.
* Reference suitability components — currently a collapsible HTML table; should
  become a traffic-light component bar.
* Benchmark status + composite — `to_html`; should become a status plot + an
  executed-only leaderboard.

## 5. Are colors consistent across figures?

* A hierarchical palette **exists** (`plotting/palette.py`:
  `build_hierarchical_color_map`, `save_hierarchical_color_map`,
  `load_hierarchical_color_map`) and a run-level `color_map.json` /
  `cell_type_color_map.tsv` are written under `outputs/hierarchical/`.
* **Gap:** the palette is not demonstrably *threaded through every figure* in the
  report — figures are collected as pre-rendered PNG/HTML from each `figures/`
  dir, so color consistency depends on whichever script produced them, not on a
  single enforced color map. The redesign must guarantee one saved color map is
  reused by all plotting and on regeneration (reuse existing, only assign new
  labels).

## 6. Are colors hierarchical?

* Yes by design: `build_hierarchical_color_map` assigns a distinct base color per
  broad family and related shades per fine subtype, with neutral grey for
  `Other` / `Unresolved` / low-confidence (`_blend` / `_family_shades`). The
  mechanism is correct; the brief's requirement (Part 6) is mostly about
  *enforcing reuse and testing the relationships*, not building it from scratch.

## 7. Are broad cell types visually distinguishable?

* At the palette level yes (distinct base hues per family). Not independently
  verified inside the rendered report because figures are pre-rendered upstream.

## 8. Are fine subpopulations represented as subtones of their broad family?

* At the palette level yes (`_family_shades` produces a ramp from the base
  color). Needs a test asserting fine colors are *near* their broad base in color
  space, and a figure (signature heatmap column annotation / hierarchy tree) that
  visibly uses it.

## 9. Are warnings visible before results?

* **[partially fixed in Stage 0]** Warnings are now *populated honestly*
  (Stage 0 fixed the false "No warnings": the warnings section now aggregates
  bulk QC recommendations, spillover risk, low gene overlap, separability
  ("N of M pairs not separable"), and unresolved families — 4 warning cards in
  the regenerated report).
* **Gap (redesign):** the Warnings section is still **#9 — after all results**.
  The brief requires the *decision summary and QC to precede results*; warnings
  and suitability status must be surfaced at the top (executive decision summary)
  and QC sections (reference, signature, input) must come before bulk/spatial.

## 10. Are reference QC and signature QC shown before deconvolution?

* **Reference QC:** yes — Reference quality is section 2 (before bulk/spatial),
  and it already shows a suitability badge + components.
* **Signature / separability / hierarchy QC:** **no** — these live in sections 6
  (hierarchical) and 7 (resolution/separability/spillover), **after** the bulk
  (4) and spatial (5) results. This is the **main structural defect** the
  redesign targets: signature quality, separability, and hierarchy usability must
  be promoted to a QC dashboard *before* results.

## 11. Does the report help decide whether predictions are trustworthy?

* Partially. The pieces exist (suitability badge, separability/spillover tables,
  unresolved mass, populated warnings) but they are **scattered and mostly after
  the results**, and there is no single executive decision view (status cards +
  "before interpreting" checklist + per-family recommended interpretation level).
  A reader cannot, from the top of the report, see "reference = CAUTION, fine
  separability = WARNING → trust family-level, not subtype-level" without
  scrolling through results first.

---

## Summary of required redesign changes (derived from the audit)

**Already addressed in Stage 0 (do not redo):** empty figure bodies (now
embedded), false "No warnings" (now honest), composite-score spatial leak (now
per-modality).

**Structural (highest priority):**
1. Add an **executive decision summary** at the top: status cards
   (reference suitability / hierarchy / separability / input compatibility =
   PASS/CAUTION/WARNING/FAIL; bulk/spatial/benchmark availability = yes/no) and a
   "before interpreting results" checklist, with neutral, non-biological framing.
2. **Reorder** so QC precedes results: Reference QC → Signature/hierarchy QC →
   Input compatibility **before** Bulk/Spatial results. Promote separability,
   spillover, and hierarchy-usability content into a signature-QC dashboard ahead
   of results (keep a results-adjacent resolution/uncertainty section too).

**Content / quality:**
3. Replace the generic `_default` caption; give every emitted figure a specific
   caption (title, axes, colors, method, source, interpretation, limitation).
4. Replace large tables with summary plots (top-20 confusable pairs, spillover
   network, unresolved-mass barplot, suitability component bar, benchmark
   status + executed-only leaderboard); keep full tables as collapsible/source
   data only.
5. Enforce one saved hierarchical color map reused across all figures and on
   regeneration; add palette tests (distinct broad bases, related fine subtones,
   determinism, reuse).

**Tests (Part 17):** QC-before-results ordering; signature QC before
bulk/spatial; no empty figure cards; no generic captions; suitability
WARNING/FAIL surfaced before results; spatial section leads with H&E alignment;
broad-before-fine; no biological-causality claims; no spatial-accuracy claims
without ground truth; deterministic color map.

**Out-of-scope reminders (project rules):** do not modify core deconvolution
algorithms, do not change prediction values, do not download data, do not commit
outputs, do not overinterpret biology, do not claim spatial accuracy without
ground truth.
