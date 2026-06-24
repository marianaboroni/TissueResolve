# Tutorial — Adaptive resolution & reference uncertainty

TissueResolve writes identifiability-aware diagnostics on every run (reporting only —
they do **not** change the estimates). They tell you which fine cell states you can
trust and which should be interpreted at a coarser level. Outputs live under
`<out>/resolution/`.

## The files

| file | what it tells you |
|---|---|
| `adaptive_resolution.tsv` | per broad family: the best *supported* resolution |
| `reference_uncertainty.tsv` | per state: cells/donors, donor variability, marker stability, reliability |
| `state_reliability.tsv` | per state: reliability score in [0,1] |
| `family_reliability.tsv` | per family: mean/min reliability, total cells |

## adaptive_resolution.tsv — the resolution categories

`family | recommended_resolution | supported_states | grouped_states | unresolved_mass | reason | confidence`

- **resolved_fine** — subtypes are separable; trust the fine estimates.
- **partially_resolved_group** — some subtypes separable, others collinear: trust the
  separable ones; interpret the listed `grouped_states` together, not individually.
- **broad_only** — all subtypes mutually collinear; interpret at the family level only.
- **unresolved_family** — most family mass is `unresolved_<family>` (gating abstained).
- **diagnostic_only** — borderline; inspect before trusting any fine split.

## What drives the category

- **Collinearity** — pairwise Bhattacharyya separability + within-family signature
  condition number. Collinear fine states (e.g. closely related T-cell subtypes) cannot
  be split reliably from a mean-profile reference and are grouped/broad.
- **Unresolved mass** — `unresolved_<family>` columns in the estimates: a family
  dominated by unresolved mass is `unresolved_family`.
- **Low support** — few cells / few donors / high donor-to-donor variability lower a
  state's reliability (`reference_uncertainty.tsv`), explaining why a rare state may be
  `diagnostic_only`.

## Why some fine estimates should not be trusted

A fine proportion can be numerically produced yet not be identifiable: if two states
have near-identical reference signatures, the split between them is noise. The adaptive
resolution table flags exactly these cases so you report the family (or the group)
rather than an overconfident subtype. This is intentional — TissueResolve declines to
promise fine resolution the reference cannot support (see the negative results in
`docs/PERFORMANCE_BENCHMARK_REPORT.md` §10).

## Example reading

```
family   recommended_resolution     supported_states   grouped_states         unresolved_mass  confidence
T/NK     partially_resolved_group   NK cell            CD4 T | CD8 T | Treg   0.18             0.42
Myeloid  broad_only                                    Macro | cDC | Mono     0.05             0.20
Epithelial resolved_fine            Luminal; Basal                            0.01             0.91
```

Here NK is trustworthy; the CD4/CD8/Treg split is not (report them as a group); the
myeloid subtypes are all collinear (report at family level); epithelial subtypes are
resolvable.
