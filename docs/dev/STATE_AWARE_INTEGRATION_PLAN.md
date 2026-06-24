# State-aware integration plan (Part 1)

Read-only plan, written before any code change. Goal: **integrate** the already-
built granular-signature + state-aware components behind one explicit
experimental flag, **measure** before/after, and **report honestly** — without
expanding the package or implementing the deferred features.

## Already implemented (committed)
- `reference/three_level_hierarchy.py` — broad→cell_type→state model (+ 2-level
  fallback) — `60efe97`.
- `reference/granular_signatures.py` — broad / cell-type / state gene panels via
  reuse of `build_family_specific_gene_panels` — `409e9f2`.
- `bulk/state_aware_hierarchical.py` — `run_state_aware_hierarchical_bulk`
  (broad→cell_type→state, correct-level unresolved mass, optional panels) —
  `409e9f2`. Tested (57 affected tests green).

## Not integrated (the gap this step closes)
- The state-aware solver is **not reachable** from `api.deconv_bulk` / the CLI.
- No `outputs/signatures/*` are generated in the real workflow.
- No benchmark row for the state-aware strategy.
- No report subsection; no feature-status doc.

## What will be wired now (minimal, reuse-first)
1. **api.py** — add `state_aware: bool = False` (+ optional `reference_adata`)
   to `deconv_bulk`. When enabled in hierarchical mode: build a
   `ThreeLevelHierarchy` (2-level fallback from the `{fine→broad}` mapping when no
   state labels), optionally build granular panels from `reference_adata` (else
   `panels=None` → global genes), call `run_state_aware_hierarchical_bulk`, and
   record metadata (`hierarchy_mode`, `state_aware_enabled`, `broad_col`,
   `cell_type_col`, `state_col`, `fallback_reason`, `feature_status`). Default
   path unchanged.
2. **cli.py** — add `--state-aware` flag (default off) → passes through.
3. **Real workflow signatures** — a small helper (reusing
   `granular_signatures.build_multigranularity_panels` +
   `write_granular_signature_outputs`) that runs **only when the reference
   AnnData (h5ad) is present**; writes `outputs/signatures/*` +
   `family_solver_gene_panels.tsv` + `granularity_signature_metadata.json`.
   Degrades gracefully (records "no AnnData / no state labels") otherwise.
4. **Benchmark** — register `TissueResolve_state_aware` (and keep
   `TissueResolve_hierarchical` as standard) as a distinct bulk method row;
   reuse existing metric functions; write
   `state_aware_integration_benchmark.tsv` + family summary + solver metadata.
5. **Report** — a "State-aware / multi-granularity deconvolution" subsection that
   reads the signature outputs + solver metadata if present and **states plainly
   when no third-level state labels were available** (the breast-cancer reference
   has only broad + fine, so this runs in two-level fallback).
6. **docs/FEATURE_STATUS.md** — classify every feature stable/experimental/
   partial/planned/not-implemented.

## Deferred (explicitly NOT in this step)
Reference adaptation; cell-type-specific expression reconstruction. README/docs
must not claim these exist.

## Tests to add
CLI/config (state-aware enable; default unchanged; 2-level fallback +
`fallback_reason` recorded); api routing to the state-aware solver; benchmark row
present and distinct; report subsection states whether state labels exist;
guardrails (deferred features remain not-implemented / not claimed).

## Risks
- The real breast-cancer reference has **no state labels** → state level is a
  2-level fallback; "state-level accuracy" cannot be measured on it. Report must
  say so (no overclaim).
- Granular panels need per-cell AnnData (only at build time); when absent the
  solver runs with global genes — must be labelled, not silently equated with
  "granular".
- Real-data benchmark numbers require the (git-ignored, downloaded) reference;
  if absent, document that real-data measurement was not run rather than faking
  it.
- Keep behind one flag; default behavior and all prior tests must stay green.
