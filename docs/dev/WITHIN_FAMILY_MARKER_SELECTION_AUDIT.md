# Within-family marker selection — audit

Read-only audit (Part 1) of how genes/markers are currently selected, written
**before** any code change. Citations are to the current source.

## Summary of findings

| # | Question | Answer |
|---|----------|--------|
| 1 | Genes selected globally, within-family, or both? | **Global only** |
| 2 | HVGs recalculated within each broad family? | **No — HVGs are never computed at all** |
| 3 | Pairwise discriminative genes used in the solver or only diagnostics? | **Diagnostics + global panel augmentation only; never family-specific in the solver** |
| 4 | Does the fine-level step use family-specific gene panels? | **No — it reuses the global panel** |
| 5 | Why only 2 of 8 families resolved? | **Resolvability is judged on the global gene set, where within-family signal is diluted** |
| 6 | Which families lack within-family markers? | T/NK, B/Plasma, Myeloid, Epithelial, Endothelial, Mural (6 of 8) flagged unresolved |
| 7 | What must change? | Build family-specific panels from raw cells; thread them into resolvability + conditional fine estimation |

## Details

### 1. Global vs within-family
Marker/signature selection is **global**: `GeneSelector.select()`
(`reference/markers.py:138-371`) scores all genes against all cell types
(Phase 1 composite scoring) and adds pairwise-discriminative genes across **all**
pairs (Phase 2), then takes a single global top-N
(`markers.py:305`). No family stratification.

### 2. Highly variable genes
**HVGs are never computed** anywhere (no `scanpy.pp.highly_variable_genes`, no
per-cell variance ranking). The reference is built from **aggregated
per-cell-type mean profiles** (`reference/build.py` → CPM/log means); a global
cross-donor CV (`donor_cv`) is the only variability signal, and it is global,
not family-stratified.

### 3. Pairwise discriminative genes
`pairwise_markers.py` provides `score_pairwise_markers`,
`select_within_family_discriminative_genes`, `build_family_specific_gene_panels`
— but these operate on the **aggregated `ReferenceSignature`** and are used only
(a) to augment the **global** panel (`markers.py:326-330`) and (b) in
**diagnostics** (separability / resolution / report). They are **never** passed
to the bulk or spatial hierarchical solver.

### 4. Fine-level gene panel
`run_hierarchical_bulk` (`bulk/hierarchical.py`) runs the family-level and
fine-level deconvolutions with the **same** `run_kwargs`/global panel, then
computes conditional within-family proportions **post-hoc**
(`hierarchy.compute_conditional_subtype_proportions`). The within-family NNLS
(`bulk/solver.py`) subsets `phi` to `gene_panel ∩ bulk ∩ ref`, but `gene_panel`
is **global**. So fine subtypes are separated using genes chosen to separate
**broad families**, not subtypes.

### 5. Why only 2/8 resolved
`evaluate_within_family_resolvability` (`hierarchy.py:675-772`) gates each family
on mean separability (1 − Bhattacharyya), worst-pair discriminating-gene count
(|log2FC|>1), and within-family spillover (Pearson r) — **all computed on the
GLOBAL gene set** (the family view keeps `gene_names = all genes`). Because the
global panel is optimised for between-family contrast, within-family pairs look
inseparable and 6/8 families are gated to unresolved.

### 6. Families lacking within-family markers (current verdict)
From the latest run: resolved = 2 of 8; unresolved = **T/NK, B/Plasma, Myeloid,
Epithelial, Endothelial, Mural**. This is a *consequence of measuring
separability on global genes*, not necessarily a true biological ceiling.

### 7. Data-flow constraint (critical) and what must change
- **Raw per-cell AnnData is available only at reference-build time**
  (`build.py` reads the h5ad, aggregates, and **discards** the cells). The
  solver only ever sees the aggregated `ReferenceSignature`
  (`gene_names, cell_types, R_cpm, R_log, phi, donor_cv`). There is **no**
  per-family gene-panel field on `ReferenceSignature` or the hierarchy bundle.
- Therefore **within-family HVG/DE selection must happen at build time** (from
  AnnData) and the resulting **family-specific gene panels must be stored** and
  threaded into the solver.

**Code that must change (additive, backward-compatible):**
1. **New** `reference/within_family_markers.py` — AnnData-based within-family HVG
   + DE + pairwise selection and a `build_family_specific_gene_panels(adata,…)`
   that returns panels + per-gene scores + excluded genes + marker-support
   tables.
2. **`reference/separability.py`** (or a new helper) — compute separability
   **restricted to each family's gene panel** (within-family separability),
   alongside the existing global separability.
3. **`reference/hierarchy.py`** — `evaluate_within_family_resolvability` and
   `compute_conditional_subtype_proportions` gain an **optional**
   `family_gene_panels` argument; when provided, separability/resolvability and
   the conditional fine split are computed on the family-specific panel. Default
   (no panels) keeps current behaviour, so existing tests are unaffected.
4. **`bulk/hierarchical.py`** (and later `spatial/hierarchical.py`) — accept and
   pass `family_gene_panels` through to the conditional fine step.
5. **Report** (`07_generate_reports.py`) — a "Within-family marker selection"
   subsection with before/after separability and per-family panel sizes.

## Constraints honored
Do not rewrite the package; **keep** existing global marker selection; additive
and backward-compatible (defaults unchanged); do not weaken tests; do not hide
warnings; do not overclaim fine-subtype accuracy (only report improvement where
metrics improve); nothing committed.
