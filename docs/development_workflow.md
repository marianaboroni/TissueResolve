# Development workflow

## Environment

Use the project virtual environment and install the editable package with the
optional extras needed for spatial work, plotting, and reports:

```bash
source .venv/bin/activate
python -m pip install -e ".[all]"   # spatial + report + plots + dev
```

`scikit-learn` and `scanpy` come with the `[spatial]` extra; `matplotlib` with
`[plots]`; `jinja2` is optional (`[report]`) — the HTML reporter uses f-strings
and does not require it.

## Running tests

```bash
python -m pytest -v          # full suite
python -m pytest tests/spatial -v
python -m pytest -k compliance -v
```

### Testing rules (from `CLAUDE.md`)

- The default suite is **fast, deterministic, and offline** — no network, no
  dataset downloads.
- Unit tests use toy synthetic data; integration tests use tiny local fixtures.
- Real-data tests run only behind an explicit flag and are **not** part of the
  default suite.
- Do not weaken or skip a test to make the suite pass.
- If a bug is found, add a failing regression test first, then fix the bug.

## Staged build

The package was built in stages (see `AUDIT_AND_MIGRATION_PLAN.md`):

| Stage | Scope |
|---|---|
| 0 | Scaffold + pre-migration bug fixes |
| 1 | Shared reference layer (`ReferenceSignature`, markers, separability) |
| 2 | Protocol layer (metadata, risk, mismatch) |
| 3 | Bulk workflow (wNNLS, bootstrap, QC, pipeline, CLI) |
| 4 | Spatial workflow (graph, NB-CAR model, QC, neighbourhood, benchmark, CLI) |
| 5 | Unification: plotting, reports, methods text, docs, compliance |

Make only the changes for the current stage; do not auto-advance.

## Non-negotiable rules

Enforced by `tests/test_compliance.py` and by the result containers:

1. No silent gene removal.
2. No hidden warnings.
3. Bulk outputs are mRNA proportions, not cell fractions (unless explicit
   mRNA-content correction was applied).
4. Spatial smoothing parameters are always recorded.
5. No confident estimates for non-separable types without warnings.
6. No plots without saved underlying data.
7. No mixing of bulk and spatial model assumptions.

## Legacy code

`chimera_v1/` and `spatcar/` are the read-only legacy packages. Do not modify
them; port logic into `src/tissueresolve/` instead.
