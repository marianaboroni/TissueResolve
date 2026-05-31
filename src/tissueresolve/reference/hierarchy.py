"""
Hierarchical (broad → subtype) deconvolution support.

When subtypes within a family are not reliably separable, a hierarchical view
is more honest than forcing confident subtype splits:

1. estimate **broad group** proportions from an aggregated reference,
2. estimate **conditional subtype** proportions *within* each group,
3. combine into **absolute subtype** proportions = broad × conditional,
4. for families that are not resolvable, keep the mass at the family level as
   **unresolved family mass** rather than splitting it.

This module builds the aggregated reference and performs the arithmetic; it
does not change the core solver — the same bulk/spatial pipelines are run on
the aggregated reference.
"""
from __future__ import annotations

import re
from typing import Optional

import numpy as np
import pandas as pd

from tissueresolve.results import ReferenceSignature

__all__ = [
    "infer_broad_groups_from_labels",
    "apply_user_label_mapping",
    "aggregate_reference_by_group",
    "merge_reference_cell_types",
    "build_hierarchical_reference",
    "decompose_within_family",
    "aggregate_predictions_by_family",
    "compare_fine_vs_merged_predictions",
]

# Coarse keyword → broad group, checked in order (first match wins).
_BROAD_KEYWORDS = [
    ("T cell", "T cell"), ("CD8", "T cell"), ("CD4", "T cell"),
    ("regulatory T", "T cell"), ("B cell", "B cell"), ("plasma", "Plasma cell"),
    ("NK", "NK cell"), ("natural killer", "NK cell"),
    ("macrophage", "Myeloid"), ("monocyte", "Myeloid"), ("dendritic", "Myeloid"),
    ("mast", "Myeloid"), ("myeloid", "Myeloid"), ("neutrophil", "Myeloid"),
    ("endothelial", "Endothelial"), ("fibroblast", "Fibroblast"),
    ("pericyte", "Perivascular"), ("smooth muscle", "Perivascular"),
    ("epithelial", "Epithelial"), ("luminal", "Epithelial"),
    ("basal", "Epithelial"), ("myoepithelial", "Epithelial"),
    ("lymphocyte", "Lymphoid"),
]


def infer_broad_groups_from_labels(cell_types: list[str]) -> dict[str, str]:
    """Heuristically map fine cell-type labels to broad groups by keyword.

    This is a convenience default; for publication use an explicit mapping via
    :func:`apply_user_label_mapping`.  Unmatched labels map to themselves.
    """
    mapping: dict[str, str] = {}
    for ct in cell_types:
        low = ct.lower()
        group = None
        for kw, grp in _BROAD_KEYWORDS:
            if kw.lower() in low:
                group = grp
                break
        mapping[ct] = group if group is not None else ct
    return mapping


def apply_user_label_mapping(cell_types: list[str],
                             mapping: dict[str, str]) -> dict[str, str]:
    """Validate a user-provided fine→broad mapping; unmapped types map to self."""
    out: dict[str, str] = {}
    missing = []
    for ct in cell_types:
        if ct in mapping:
            out[ct] = str(mapping[ct])
        else:
            out[ct] = ct
            missing.append(ct)
    if missing:
        import warnings

        warnings.warn(
            f"{len(missing)} cell type(s) had no broad-group mapping and were "
            f"kept at fine level: {missing[:5]}"
            + ("…" if len(missing) > 5 else ""),
            stacklevel=2,
        )
    return out


def aggregate_reference_by_group(
    ref: ReferenceSignature, mapping: dict[str, str],
) -> ReferenceSignature:
    """Aggregate a reference to broad groups (n_cells-weighted mean of CPM).

    Returns a new :class:`ReferenceSignature` whose ``cell_types`` are the
    broad groups.  ``phi_g`` and gene names are preserved; ``donor_cv`` is
    dropped (cross-donor CV is not meaningfully aggregable across subtypes).
    """
    cts = list(ref.cell_types)
    groups = sorted({mapping.get(ct, ct) for ct in cts})
    R = ref.as_R_cpm().astype(np.float64)              # (K, G)
    n_cells = ref.n_cells_per_type or {}

    rows = []
    new_counts: dict[str, int] = {}
    for g in groups:
        members = [i for i, ct in enumerate(cts) if mapping.get(ct, ct) == g]
        weights = np.array([max(1, int(n_cells.get(cts[i], 1))) for i in members],
                           dtype=float)
        weights /= weights.sum()
        rows.append((R[members] * weights[:, None]).sum(axis=0))
        new_counts[g] = int(sum(int(n_cells.get(cts[i], 0)) for i in members))

    R_group = np.vstack(rows).astype(np.float32)       # (n_groups, G)
    broad = ReferenceSignature(
        gene_names=list(ref.gene_names),
        cell_types=groups,
        R_cpm=R_group,
        R_log=np.log1p(R_group).astype(np.float32),
        phi_g=ref.phi_g,
        n_cells_per_type=new_counts,
        genome=ref.genome,
    )
    broad.validate()
    return broad


def build_hierarchical_reference(
    ref: ReferenceSignature, mapping: Optional[dict[str, str]] = None,
) -> dict:
    """Build a broad reference plus the group→members structure.

    Returns ``{"broad": ReferenceSignature, "mapping": {...},
    "groups": {group: [members]}, "subtype_reference": ref}``.
    """
    mapping = mapping or infer_broad_groups_from_labels(list(ref.cell_types))
    groups: dict[str, list[str]] = {}
    for ct in ref.cell_types:
        groups.setdefault(mapping.get(ct, ct), []).append(ct)
    broad = aggregate_reference_by_group(ref, mapping)
    return {"broad": broad, "mapping": mapping, "groups": groups,
            "subtype_reference": ref}


def decompose_within_family(
    broad_proportions: pd.DataFrame,
    conditional_subtype_proportions: pd.DataFrame,
    mapping: dict[str, str],
    *,
    unresolved_groups: Optional[list[str]] = None,
) -> dict[str, pd.DataFrame]:
    """Combine broad and within-group conditional estimates.

    Parameters
    ----------
    broad_proportions:
        samples × groups (rows sum to 1).
    conditional_subtype_proportions:
        samples × subtypes; within each group the subtype values should sum to
        1 (conditional on the group).  Renormalised per group defensively.
    mapping:
        subtype → group.
    unresolved_groups:
        groups to keep at the broad level (their subtype mass is reported as
        ``unresolved_<group>`` rather than split).

    Returns
    -------
    dict with ``broad_proportions``, ``conditional_subtype_proportions``,
    ``absolute_subtype_proportions`` and ``unresolved_family_mass`` — total
    mass per sample is preserved.
    """
    unresolved = set(unresolved_groups or [])
    samples = broad_proportions.index
    subtypes = list(conditional_subtype_proportions.columns)
    group_of = {st: mapping.get(st, st) for st in subtypes}

    # Renormalise conditional proportions within each group.
    cond = conditional_subtype_proportions.copy().astype(float)
    for g in broad_proportions.columns:
        members = [st for st in subtypes if group_of[st] == g]
        if not members:
            continue
        sub = cond[members]
        totals = sub.sum(axis=1).replace(0, np.nan)
        cond[members] = sub.div(totals, axis=0).fillna(0.0)

    abs_rows = pd.DataFrame(0.0, index=samples, columns=subtypes)
    unresolved_mass = pd.DataFrame(0.0, index=samples,
                                   columns=[f"unresolved_{g}" for g in sorted(unresolved)])
    for g in broad_proportions.columns:
        members = [st for st in subtypes if group_of[st] == g]
        if g in unresolved:
            if f"unresolved_{g}" in unresolved_mass.columns:
                unresolved_mass[f"unresolved_{g}"] = broad_proportions[g]
            continue
        for st in members:
            abs_rows[st] = broad_proportions[g] * cond[st]

    return {
        "broad_proportions": broad_proportions,
        "conditional_subtype_proportions": cond,
        "absolute_subtype_proportions": abs_rows,
        "unresolved_family_mass": unresolved_mass,
    }


# ---------------------------------------------------------------------------
# Merged-reference / post-hoc family aggregation workflow
# ---------------------------------------------------------------------------


def merge_reference_cell_types(reference, mapping: dict[str, str]):
    """Pre-deconvolution merge: aggregate a reference to family cell types.

    Thin, explicit wrapper over :func:`aggregate_reference_by_group` — the
    mapping (fine → family) is the recorded, reproducible merge plan.  Returns
    a new :class:`~tissueresolve.results.ReferenceSignature`; the input is not
    modified.
    """
    return aggregate_reference_by_group(reference, mapping)


def aggregate_predictions_by_family(
    predictions: "pd.DataFrame", mapping: dict[str, str],
) -> "pd.DataFrame":
    """Post-hoc: sum subtype proportions into family-level proportions.

    *predictions* is samples/spots × cell types.  Columns are grouped by
    ``mapping`` (unmapped columns keep their own name) and summed, so **total
    mass per row is preserved**.  The input is never modified.
    """
    cols = list(predictions.columns)
    fam_of = {c: mapping.get(c, c) for c in cols}
    families = []
    for c in cols:                      # stable, first-seen order
        f = fam_of[c]
        if f not in families:
            families.append(f)
    out = pd.DataFrame(index=predictions.index)
    for f in families:
        members = [c for c in cols if fam_of[c] == f]
        out[f] = predictions[members].sum(axis=1)
    return out


def compare_fine_vs_merged_predictions(
    predictions: "pd.DataFrame", mapping: dict[str, str],
) -> dict:
    """Return both fine and family-level predictions plus the mapping.

    ``{"fine": <unchanged>, "family": <aggregated>, "mapping": {...},
    "merge_stage": "post_hoc_aggregation"}`` — fine estimates are never
    overwritten.
    """
    family = aggregate_predictions_by_family(predictions, mapping)
    return {
        "fine": predictions,
        "family": family,
        "mapping": dict(mapping),
        "merge_stage": "post_hoc_aggregation",
    }
