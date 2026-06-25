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
import warnings
from dataclasses import dataclass, field
from pathlib import Path
from typing import Any, Optional

import numpy as np
import pandas as pd

from tissueresolve.results import ReferenceSignature

__all__ = [
    # legacy / internal helpers (kept for backward compatibility)
    "infer_broad_groups_from_labels",
    "apply_user_label_mapping",
    "aggregate_reference_by_group",
    "merge_reference_cell_types",
    "build_hierarchical_reference",
    "decompose_within_family",
    "compare_fine_vs_merged_predictions",
    # public broad-to-fine hierarchy API
    "infer_broad_cell_type_family",
    "build_cell_type_hierarchy",
    "validate_hierarchy",
    "save_hierarchy_mapping",
    "load_hierarchy_mapping",
    "aggregate_reference_by_family",
    "aggregate_predictions_by_family",
    "compute_conditional_subtype_proportions",
    "combine_family_and_conditional_estimates",
    "add_unresolved_family_mass",
    "hierarchy_to_frame",
    "evaluate_within_family_resolvability",
    "decide_unresolved_families",
    "assemble_hierarchical_estimates",
    "HierarchicalEstimates",
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

# ---------------------------------------------------------------------------
# Canonical broad cell-type families (spec heuristics).  Used by
# ``infer_broad_cell_type_family``; checked in order, first substring match
# wins.  These are a *convenience fallback* only — for hierarchical mode an
# explicit broad/fine annotation or mapping file is strongly preferred.
# ---------------------------------------------------------------------------
FAMILY_OTHER = "Other"

_FAMILY_KEYWORDS: list[tuple[str, str]] = [
    # B/plasma (check before generic "t cell" etc.)
    ("plasma", "B/plasma"), ("iga", "B/plasma"), ("igg", "B/plasma"),
    ("memory b", "B/plasma"), ("b cell", "B/plasma"), ("b-cell", "B/plasma"),
    # T/NK
    ("regulatory t", "T/NK"), ("treg", "T/NK"), ("cd4", "T/NK"), ("cd8", "T/NK"),
    ("t cell", "T/NK"), ("t-cell", "T/NK"), ("nkt", "T/NK"),
    ("natural killer", "T/NK"), ("nk cell", "T/NK"), ("nk-cell", "T/NK"),
    ("lymphocyte", "T/NK"),
    # Myeloid
    ("macrophage", "Myeloid"), ("monocyte", "Myeloid"), ("dendritic", "Myeloid"),
    ("dc", "Myeloid"), ("neutrophil", "Myeloid"), ("myeloid", "Myeloid"),
    ("mast", "Myeloid"), ("granulocyte", "Myeloid"),
    # Endothelial
    ("endothel", "Endothelial"), ("capillary", "Endothelial"),
    ("vein", "Endothelial"), ("artery", "Endothelial"),
    ("vascular tree", "Endothelial"), ("lymphatic vessel", "Endothelial"),
    # Mural (check "smooth muscle" / "pericyte" before generic stromal)
    ("pericyte", "Mural"),
    ("vascular associated smooth muscle", "Mural"), ("smooth muscle", "Mural"),
    # Stromal / fibroblast
    ("fibroblast", "Stromal/fibroblast"), ("caf", "Stromal/fibroblast"),
    ("stromal", "Stromal/fibroblast"),
    # Epithelial
    ("mammary gland epithelial", "Epithelial"), ("luminal", "Epithelial"),
    ("basal", "Epithelial"), ("myoepithelial", "Epithelial"),
    ("epithelial", "Epithelial"),
    # Adipocyte
    ("adipocyte", "Adipocyte"), ("adipose", "Adipocyte"),
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


# ===========================================================================
# Public broad-to-fine hierarchy API
# ===========================================================================
#
# These functions implement the documented, reproducible hierarchical
# deconvolution contract used by the bulk/spatial hierarchical workflows and
# the CLI.  They never modify their inputs and always preserve total mass.
# The convention is::
#
#     mapping: dict[fine_cell_type -> broad_cell_type_family]
#
# ``fine`` labels are the reference's own ``cell_types``; ``broad`` labels are
# the families they belong to.
# ---------------------------------------------------------------------------


def infer_broad_cell_type_family(label: str) -> str:
    """Map a single fine cell-type *label* to a broad family by keyword.

    Returns one of the canonical families
    (``"T/NK"``, ``"B/plasma"``, ``"Myeloid"``, ``"Endothelial"``,
    ``"Epithelial"``, ``"Stromal/fibroblast"``, ``"Mural"``, ``"Adipocyte"``)
    or ``"Other"`` when no keyword matches.

    This is a *convenience fallback only*.  For hierarchical deconvolution an
    explicit broad/fine annotation (or a mapping file) is strongly preferred;
    callers should never silently rely on inference without warning the user.
    """
    low = str(label).strip().lower()
    if not low or low.startswith("unresolved") or low.startswith("family:"):
        return FAMILY_OTHER
    for kw, fam in _FAMILY_KEYWORDS:
        if kw in low:
            return fam
    return FAMILY_OTHER


def build_cell_type_hierarchy(
    cell_types: list[str],
    mapping: Optional[dict[str, str]] = None,
) -> dict[str, str]:
    """Build a validated fine→broad mapping for *cell_types*.

    Parameters
    ----------
    cell_types:
        Fine cell-type labels (the reference's ``cell_types``).
    mapping:
        Optional user-provided fine→broad mapping.  When ``None``, families
        are inferred per label via :func:`infer_broad_cell_type_family` (a
        warning is emitted, since inference is a fallback, not a guarantee).
        When provided, every cell type must be present; missing keys raise.

    Returns
    -------
    dict[str, str]
        fine cell type → broad family, covering every entry in *cell_types*.
    """
    cell_types = [str(c) for c in cell_types]
    if mapping is None:
        inferred = {ct: infer_broad_cell_type_family(ct) for ct in cell_types}
        n_other = sum(1 for v in inferred.values() if v == FAMILY_OTHER)
        warnings.warn(
            "build_cell_type_hierarchy: no explicit broad/fine mapping was "
            "provided; broad families were inferred heuristically from labels "
            f"({n_other}/{len(cell_types)} mapped to '{FAMILY_OTHER}').  For "
            "hierarchical deconvolution, provide explicit broad/fine "
            "annotation columns or a --cell-type-hierarchy mapping file.",
            stacklevel=2,
        )
        out = inferred
    else:
        out = {ct: str(mapping[ct]) for ct in cell_types if ct in mapping}
        missing = [ct for ct in cell_types if ct not in mapping]
        if missing:
            raise ValueError(
                "build_cell_type_hierarchy: the provided mapping does not "
                f"cover {len(missing)} cell type(s): {missing[:8]}"
                + ("…" if len(missing) > 8 else "")
                + ".  Every fine cell type must map to exactly one broad "
                "family (add them to the mapping file or use auto-inference)."
            )
    validate_hierarchy(out, cell_types)
    return out


def validate_hierarchy(mapping: dict[str, str], cell_types: list[str]) -> None:
    """Validate a fine→broad *mapping* against *cell_types*.

    Raises ``ValueError`` when:

    * a cell type has no mapping,
    * a broad family name is empty/blank.

    Emits a warning when a broad family contains only one fine subtype (such
    families cannot be split and are reported at the broad level).

    A ``dict`` mapping is inherently one-broad-per-fine, so the
    "one fine label maps to multiple broad labels" condition is impossible to
    represent here; loaders (:func:`load_hierarchy_mapping`) enforce it at
    parse time.
    """
    cell_types = [str(c) for c in cell_types]
    missing = [ct for ct in cell_types if ct not in mapping]
    if missing:
        raise ValueError(
            f"validate_hierarchy: {len(missing)} cell type(s) are not mapped "
            f"to a broad family: {missing[:8]}"
            + ("…" if len(missing) > 8 else "")
        )
    blank = [ct for ct in cell_types if not str(mapping[ct]).strip()]
    if blank:
        raise ValueError(
            f"validate_hierarchy: blank/empty broad family for: {blank[:8]}"
        )
    families: dict[str, list[str]] = {}
    for ct in cell_types:
        families.setdefault(str(mapping[ct]), []).append(ct)
    singletons = [f for f, members in families.items() if len(members) == 1]
    if singletons:
        warnings.warn(
            f"validate_hierarchy: {len(singletons)} broad family/ies contain a "
            f"single fine subtype and cannot be sub-resolved: {singletons[:8]}"
            + ("…" if len(singletons) > 8 else "")
            + ".  These will be reported at the broad level only.",
            stacklevel=2,
        )


def hierarchy_to_frame(mapping: dict[str, str]) -> pd.DataFrame:
    """Return a tidy ``cell_type, family`` DataFrame for the *mapping*."""
    rows = [{"cell_type": ct, "family": fam} for ct, fam in mapping.items()]
    return pd.DataFrame(rows, columns=["cell_type", "family"])


# accepted header aliases when loading a mapping file
_FINE_COL_ALIASES = (
    "cell_type", "fine_cell_type", "sub_cell_type", "subtype",
    "cell_state", "fine", "fine_label",
)
_BROAD_COL_ALIASES = (
    "family", "broad_cell_type", "broad", "major_cell_type",
    "compartment", "lineage", "broad_label",
)


def save_hierarchy_mapping(mapping: dict[str, str], path: Path | str) -> Path:
    """Persist a fine→broad *mapping* as a two-column TSV (``cell_type``/``family``)."""
    path = Path(path)
    path.parent.mkdir(parents=True, exist_ok=True)
    hierarchy_to_frame(mapping).to_csv(path, sep="\t", index=False)
    return path


def load_hierarchy_mapping(path: Path | str) -> dict[str, str]:
    """Load a fine→broad mapping TSV.

    The file must have a header with one fine-label column and one broad-label
    column.  Accepted fine headers: ``cell_type``/``fine_cell_type``/
    ``sub_cell_type``/…; broad headers: ``family``/``broad_cell_type``/
    ``compartment``/….  Raises ``ValueError`` if the columns cannot be
    identified or if a fine label maps to more than one broad family.
    """
    path = Path(path)
    if not path.exists():
        raise FileNotFoundError(f"hierarchy mapping file not found: {path}")
    sep = "\t" if path.suffix.lower() in (".tsv", ".txt") else ","
    df = pd.read_csv(path, sep=sep, dtype=str, comment="#").fillna("")
    cols_lower = {c.lower().strip(): c for c in df.columns}
    fine_col = next((cols_lower[a] for a in _FINE_COL_ALIASES if a in cols_lower), None)
    broad_col = next((cols_lower[a] for a in _BROAD_COL_ALIASES if a in cols_lower), None)
    if fine_col is None or broad_col is None:
        raise ValueError(
            f"load_hierarchy_mapping: could not identify fine/broad columns in "
            f"{path} (columns: {list(df.columns)}).  Expected a fine column "
            f"(one of {_FINE_COL_ALIASES}) and a broad column "
            f"(one of {_BROAD_COL_ALIASES})."
        )
    mapping: dict[str, str] = {}
    conflicts: dict[str, set] = {}
    for fine, broad in zip(df[fine_col].astype(str), df[broad_col].astype(str)):
        fine, broad = fine.strip(), broad.strip()
        if not fine or not broad:
            continue
        if fine in mapping and mapping[fine] != broad:
            conflicts.setdefault(fine, {mapping[fine]}).add(broad)
        mapping[fine] = broad
    if conflicts:
        detail = "; ".join(f"{k!r}→{sorted(v)}" for k, v in list(conflicts.items())[:5])
        raise ValueError(
            "load_hierarchy_mapping: every fine cell type must map to exactly "
            f"one broad family, but found conflicting mappings: {detail}."
        )
    if not mapping:
        raise ValueError(f"load_hierarchy_mapping: no usable rows in {path}.")
    return mapping


def aggregate_reference_by_family(
    reference: ReferenceSignature, mapping: dict[str, str],
) -> ReferenceSignature:
    """Aggregate a *reference* to broad families (n_cells-weighted CPM mean).

    Spec-named alias of :func:`aggregate_reference_by_group`.  The input is
    never modified; gene order is preserved; ``donor_cv`` is dropped (not
    meaningfully aggregable across heterogeneous subtypes).
    """
    return aggregate_reference_by_group(reference, mapping)


def compute_conditional_subtype_proportions(
    fine_props: pd.DataFrame,
    family_props: pd.DataFrame,
    mapping: dict[str, str],
) -> pd.DataFrame:
    """Conditional within-family subtype proportions ``P(subtype | family)``.

    For each broad family, the fine-level estimates of its member subtypes are
    renormalised to sum to 1 within that family (per sample/spot).  When the
    fine estimates for a family sum to zero for a given sample, the family's
    subtypes are split *uniformly* (rather than left undefined) so that
    ``combine_family_and_conditional_estimates`` preserves the family's mass.

    Parameters
    ----------
    fine_props:
        ``(n_obs × n_fine)`` flat fine-level estimates.
    family_props:
        ``(n_obs × n_family)`` broad-level estimates.  Only its column set is
        used (to know which families exist); the magnitudes are applied later
        in :func:`combine_family_and_conditional_estimates`.
    mapping:
        fine → broad mapping.

    Returns
    -------
    pd.DataFrame
        ``(n_obs × n_fine)`` conditional proportions; within each family the
        member columns sum to 1 (per row).
    """
    subtypes = [str(c) for c in fine_props.columns]
    group_of = {st: mapping.get(st, st) for st in subtypes}
    families = list(family_props.columns)
    cond = fine_props.copy().astype(float)
    for fam in families:
        members = [st for st in subtypes if group_of[st] == fam]
        if not members:
            continue
        sub = cond[members]
        totals = sub.sum(axis=1)
        # uniform split where the family has no fine signal
        zero = totals <= 0
        if zero.any():
            cond.loc[zero, members] = 1.0 / len(members)
            totals = cond[members].sum(axis=1)
        cond[members] = cond[members].div(totals, axis=0).fillna(0.0)
    return cond


def combine_family_and_conditional_estimates(
    family_props: pd.DataFrame,
    conditional_props: pd.DataFrame,
    mapping: dict[str, str],
) -> pd.DataFrame:
    """Final absolute subtype estimates = family proportion × conditional.

    ``fine(subtype) = family(broad) × P(subtype | broad)``.  Families present
    in *family_props* but absent from the mapping's value set are passed
    through as their own single column (so flat family-only labels survive).

    Returns
    -------
    pd.DataFrame
        ``(n_obs × n_fine)`` absolute subtype proportions.  Together with any
        unresolved-mass columns the per-row total equals the family-level total.
    """
    subtypes = [str(c) for c in conditional_props.columns]
    group_of = {st: mapping.get(st, st) for st in subtypes}
    out = pd.DataFrame(0.0, index=family_props.index, columns=subtypes)
    for fam in family_props.columns:
        members = [st for st in subtypes if group_of[st] == fam]
        for st in members:
            out[st] = family_props[fam] * conditional_props[st]
    return out


def add_unresolved_family_mass(
    fine_props: pd.DataFrame,
    family_props: pd.DataFrame,
    mapping: dict[str, str],
    unresolved_families: list[str],
) -> tuple[pd.DataFrame, pd.DataFrame]:
    """Move unresolved families' mass out of fine subtypes into ``unresolved_*``.

    For every family in *unresolved_families*, its member subtype estimates in
    *fine_props* are zeroed and the family's broad-level mass is reported in a
    dedicated ``unresolved_<family>`` column instead — never split into
    subtypes for which there is no evidence.

    Parameters
    ----------
    fine_props:
        ``(n_obs × n_fine)`` absolute subtype proportions (e.g. the output of
        :func:`combine_family_and_conditional_estimates`).
    family_props:
        ``(n_obs × n_family)`` broad-level proportions.
    mapping:
        fine → broad mapping.
    unresolved_families:
        Broad families to keep at the broad level.

    Returns
    -------
    (resolved_fine_props, unresolved_mass)
        *resolved_fine_props* has the unresolved families' subtype columns set
        to 0; *unresolved_mass* is ``(n_obs × n_unresolved)`` with one
        ``unresolved_<family>`` column per unresolved family.  Total mass per
        row (resolved fine + unresolved) is preserved.
    """
    unresolved = [f for f in unresolved_families]
    subtypes = [str(c) for c in fine_props.columns]
    group_of = {st: mapping.get(st, st) for st in subtypes}
    resolved = fine_props.copy().astype(float)
    cols = [f"unresolved_{f}" for f in unresolved]
    unresolved_mass = pd.DataFrame(0.0, index=fine_props.index, columns=cols)
    for fam in unresolved:
        members = [st for st in subtypes if group_of[st] == fam]
        for st in members:
            resolved[st] = 0.0
        col = f"unresolved_{fam}"
        if fam in family_props.columns:
            unresolved_mass[col] = family_props[fam].astype(float)
        else:  # fall back to the (pre-zeroing) summed subtype mass
            unresolved_mass[col] = fine_props[members].sum(axis=1).astype(float)
    return resolved, unresolved_mass


# ===========================================================================
# Within-family resolvability gating and result assembly
# ===========================================================================
#
# The scientific heart of hierarchical mode: a family's mass is only split
# into fine subtypes when there is *evidence* the subtypes are distinguishable
# within that family.  The gate combines three independent signals so a single
# noisy metric cannot force (or block) a split.
# ---------------------------------------------------------------------------


def evaluate_within_family_resolvability(
    fine_ref: ReferenceSignature,
    mapping: dict[str, str],
    *,
    unresolved_threshold: float = 0.10,
    min_discriminating_genes: int = 10,
    within_family_spillover_threshold: float = 0.30,
    family_gene_panels: Optional[dict[str, list[str]]] = None,
) -> pd.DataFrame:
    """Score how separable the fine subtypes are *within* each broad family.

    For every family with ≥ 2 fine members, the pairwise separability of its
    members is computed on the fine reference (Bhattacharyya coefficient,
    Pearson correlation, number of discriminating genes).  A family is flagged
    *resolvable* only when **all** of the following hold:

    * mean within-family separability score (``1 − BC``) ≥
      ``unresolved_threshold``;
    * the worst within-family pair has ≥ ``min_discriminating_genes``
      discriminating genes (|log2FC| > 1);
    * mean within-family spillover (max Pearson r to a sibling) ≤
      ``within_family_spillover_threshold``.

    Single-member families are trivially "resolvable" (fine == family).

    Returns
    -------
    pd.DataFrame
        Indexed by broad family with columns ``n_subtypes``,
        ``mean_separability``, ``min_discriminating_genes``,
        ``mean_spillover``, ``resolvable`` (bool) and ``reason`` (str).
    """
    from tissueresolve.reference.separability import compute_separability

    cell_types = [str(c) for c in fine_ref.cell_types]
    fam_members: dict[str, list[str]] = {}
    for ct in cell_types:
        fam_members.setdefault(str(mapping.get(ct, ct)), []).append(ct)

    rows: dict[str, dict[str, Any]] = {}
    for fam, members in fam_members.items():
        if len(members) <= 1:
            rows[fam] = {
                "n_subtypes": len(members),
                "n_panel_genes": len(fine_ref.gene_names),
                "mean_separability": float("nan"),
                "min_discriminating_genes": -1,
                "mean_spillover": float("nan"),
                "resolvable": True,
                "reason": "single subtype (fine == family)",
            }
            continue
        # When a family-specific gene panel is supplied, restrict the reference
        # to that panel BEFORE computing separability — so within-family
        # resolvability is judged on genes selected to separate the subtypes,
        # not the global (broad-family) panel.  Default (no panel) keeps the
        # original global-gene behaviour.
        ref_g = fine_ref
        n_panel = len(fine_ref.gene_names)
        if family_gene_panels:
            panel = [g for g in (family_gene_panels.get(fam) or [])
                     if g in set(fine_ref.gene_names)]
            if len(panel) >= 2:
                ref_g = fine_ref.subset_genes(panel)
                n_panel = len(panel)
        g_cell_types = [str(c) for c in ref_g.cell_types]
        idx = [g_cell_types.index(m) for m in members]
        sub_ref = ReferenceSignature(
            gene_names=list(ref_g.gene_names),
            cell_types=members,
            R_cpm=(ref_g.as_R_cpm()[idx]).astype(np.float32),
            R_log=(ref_g.as_R_log()[idx]).astype(np.float32),
            phi=(ref_g.as_phi()[:, idx]) if ref_g.phi is not None else None,
            phi_g=ref_g.phi_g,
            n_cells_per_type={m: fine_ref.n_cells_per_type.get(m, 0) for m in members},
            genome=fine_ref.genome,
        )
        with warnings.catch_warnings():
            warnings.simplefilter("ignore")
            sep = compute_separability(sub_ref, warn_threshold=2.0,
                                       raise_on_critical=False)
        seps = [1.0 - p.bhattacharyya_coeff for p in sep.pairs]
        disc = [p.n_discriminating_genes for p in sep.pairs]
        spill = [max(p.pearson_r, 0.0) for p in sep.pairs]
        mean_sep = float(np.mean(seps)) if seps else float("nan")
        min_disc = int(min(disc)) if disc else 0
        mean_spill = float(np.mean(spill)) if spill else float("nan")

        reasons = []
        if not (mean_sep >= unresolved_threshold):
            reasons.append(
                f"mean separability {mean_sep:.3f} < {unresolved_threshold}")
        if not (min_disc >= min_discriminating_genes):
            reasons.append(
                f"min discriminating genes {min_disc} < {min_discriminating_genes}")
        if not (mean_spill <= within_family_spillover_threshold):
            reasons.append(
                f"within-family spillover {mean_spill:.3f} > "
                f"{within_family_spillover_threshold}")
        resolvable = not reasons
        rows[fam] = {
            "n_subtypes": len(members),
            "n_panel_genes": n_panel,
            "mean_separability": round(mean_sep, 4),
            "min_discriminating_genes": min_disc,
            "mean_spillover": round(mean_spill, 4),
            "resolvable": resolvable,
            "reason": "resolvable" if resolvable else "; ".join(reasons),
        }
    out = pd.DataFrame(rows).T
    # single-subtype rows lack n_panel_genes; fill for a stable schema
    if "n_panel_genes" not in out.columns:
        out["n_panel_genes"] = len(fine_ref.gene_names)
    return out.reindex(
        columns=["n_subtypes", "n_panel_genes", "mean_separability",
                 "min_discriminating_genes", "mean_spillover", "resolvable",
                 "reason"]
    )


def compute_within_family_subtype_confidence(
    fine_ref: ReferenceSignature, mapping: dict[str, str],
    *, min_discriminating_genes: int = 10,
    family_gene_panels: Optional[dict[str, list[str]]] = None,
) -> dict[str, dict[str, Any]]:
    """Per-subtype confidence *within* its family.

    A subtype is confident when it is well separated from its **closest** family
    sibling (high ``1 − BC`` and enough discriminating genes).  Returns
    ``{subtype: {"confidence": float, "min_disc_genes": int, "family": str}}``.
    Singletons get confidence 1.0 (fine == family).

    When *family_gene_panels* is supplied, separability is computed on each
    family's own gene panel (default ``None`` = global genes, unchanged).
    """
    from tissueresolve.reference.separability import compute_separability

    cell_types = [str(c) for c in fine_ref.cell_types]
    fam_members: dict[str, list[str]] = {}
    for ct in cell_types:
        fam_members.setdefault(str(mapping.get(ct, ct)), []).append(ct)

    out: dict[str, dict[str, Any]] = {}
    for fam, members in fam_members.items():
        if len(members) == 1:
            out[members[0]] = {"confidence": 1.0, "min_disc_genes": -1, "family": fam}
            continue
        ref_g = fine_ref
        if family_gene_panels:
            panel = [g for g in (family_gene_panels.get(fam) or [])
                     if g in set(fine_ref.gene_names)]
            if len(panel) >= 2:
                ref_g = fine_ref.subset_genes(panel)
        g_cell_types = [str(c) for c in ref_g.cell_types]
        idx = [g_cell_types.index(m) for m in members]
        sub_ref = ReferenceSignature(
            gene_names=list(ref_g.gene_names), cell_types=members,
            R_cpm=(ref_g.as_R_cpm()[idx]).astype(np.float32),
            R_log=(ref_g.as_R_log()[idx]).astype(np.float32),
            phi=(ref_g.as_phi()[:, idx]) if ref_g.phi is not None else None,
            phi_g=ref_g.phi_g, genome=fine_ref.genome,
            n_cells_per_type={m: fine_ref.n_cells_per_type.get(m, 0) for m in members})
        with warnings.catch_warnings():
            warnings.simplefilter("ignore")
            sep = compute_separability(sub_ref, warn_threshold=2.0,
                                       raise_on_critical=False)
        # closest-sibling separability + discriminating genes per subtype
        best_sep = {m: 1.0 for m in members}
        best_disc = {m: 10 ** 9 for m in members}
        for p in sep.pairs:
            s = 1.0 - p.bhattacharyya_coeff
            for a in (p.type_a, p.type_b):
                if s < best_sep[a]:
                    best_sep[a] = s
                if p.n_discriminating_genes < best_disc[a]:
                    best_disc[a] = p.n_discriminating_genes
        for m in members:
            disc = 0 if best_disc[m] == 10 ** 9 else int(best_disc[m])
            conf = float(best_sep[m]) if disc >= min_discriminating_genes else \
                float(best_sep[m]) * 0.5
            out[m] = {"confidence": round(conf, 4), "min_disc_genes": disc,
                      "family": fam}
    return out


def estimate_partial_subtype_resolution(
    family: str, members: list[str], family_mass: pd.Series,
    conditional: pd.DataFrame, confidence: dict[str, dict],
    *, subtype_confidence_threshold: float = 0.10,
) -> tuple[pd.DataFrame, pd.Series]:
    """Split a family into confident subtype mass + an unresolved residual.

    Confident subtypes (confidence ≥ threshold) receive
    ``family_mass × P(subtype | family)``; the remaining (uncertain) conditional
    share becomes the family's unresolved residual.  Mass is preserved:
    ``confident subtype mass + residual = family_mass``.
    """
    confident = [m for m in members
                 if confidence.get(m, {}).get("confidence", 0.0)
                 >= subtype_confidence_threshold]
    sub = pd.DataFrame(0.0, index=family_mass.index, columns=members)
    for m in confident:
        sub[m] = family_mass * conditional[m]
    residual = family_mass * (1.0 - conditional[confident].sum(axis=1)) \
        if confident else family_mass.copy()
    return sub, residual


def estimate_soft_subtype_resolution(
    family: str, members: list[str], family_mass: pd.Series,
    conditional: pd.DataFrame, confidence: dict[str, dict],
    *, min_confidence: float = 0.0, max_confidence: float = 1.0,
) -> tuple[pd.DataFrame, pd.Series]:
    """Partial confidence-WEIGHTED resolution (the validated default gate).

    Unlike :func:`estimate_partial_subtype_resolution` (a binary keep/drop
    threshold), each subtype's mass is scaled *continuously* by its calibrated
    confidence ``c_k ∈ [0,1]``:

        resolved_k = family_mass · P(k|family) · c_k
        residual   = family_mass − Σ_k resolved_k        (≥ 0 since c_k ≤ 1)

    ``c ∈ {0,1}`` reproduces the legacy hard gate; intermediate ``c`` recovers the
    mass the threshold over-discards.  Mass is conserved exactly
    (``Σ_k resolved_k + residual = family_mass``).  Validated on breast and lung
    benchmarks (8–9/9 prospective gates; see SECOND_TISSUE_LUNG_VALIDATION_REPORT).
    """
    import numpy as _np
    sub = pd.DataFrame(0.0, index=family_mass.index, columns=members)
    for m in members:
        c = float(confidence.get(m, {}).get("confidence", 0.0))
        c = float(_np.clip(c, min_confidence, max_confidence))
        sub[m] = family_mass * conditional[m] * c
    residual = (family_mass - sub.sum(axis=1)).clip(lower=0.0)
    return sub, residual


def decide_unresolved_families(
    resolvability: pd.DataFrame, *, allow_unresolved: bool = True,
) -> list[str]:
    """Return the families to keep at the broad level (unresolved).

    A family is unresolved when ``resolvable`` is ``False`` **and**
    ``allow_unresolved`` is True.  When ``allow_unresolved`` is False, an empty
    list is returned (fine splits are forced) — callers should warn the user.
    """
    if not allow_unresolved:
        return []
    if resolvability is None or resolvability.empty:
        return []
    mask = ~resolvability["resolvable"].astype(bool)
    return [str(f) for f in resolvability.index[mask]]


@dataclass
class HierarchicalEstimates:
    """Container for the full set of broad-to-fine estimates (one modality).

    All proportion frames are indexed by sample (bulk) or spot (spatial).
    Total mass per row is preserved: ``fine_proportions`` (resolved subtypes)
    plus ``unresolved_mass`` sum to 1 per row.
    """

    family_proportions: pd.DataFrame
    conditional_proportions: pd.DataFrame
    fine_proportions: pd.DataFrame
    unresolved_mass: pd.DataFrame
    combined_fine: pd.DataFrame  # fine_proportions + unresolved_mass columns
    resolvability: pd.DataFrame
    mapping: dict[str, str]
    unresolved_families: list[str]
    qc: pd.DataFrame
    metadata: dict[str, Any] = field(default_factory=dict)


def assemble_hierarchical_estimates(
    family_props: pd.DataFrame,
    fine_props: pd.DataFrame,
    fine_ref: ReferenceSignature,
    mapping: dict[str, str],
    *,
    allow_unresolved: bool = True,
    unresolved_threshold: float = 0.10,
    min_discriminating_genes: int = 10,
    within_family_spillover_threshold: float = 0.30,
    allow_partial_resolution: bool = True,
    subtype_confidence_threshold: float = 0.10,
    gating: str = "soft",
    gating_version: str = "soft_gating-1.0",
    family_gene_panels: Optional[dict[str, list[str]]] = None,
    extra_metadata: Optional[dict] = None,
) -> HierarchicalEstimates:
    """Combine broad + fine estimates into the final hierarchical result.

    Given flat *family_props* (broad-level estimates) and flat *fine_props*
    (subtype-level estimates from the same query), this:

    1. evaluates within-family resolvability from *fine_ref*,
    2. decides which families to keep unresolved,
    3. computes conditional ``P(subtype | family)`` from the fine estimates,
    4. combines into absolute subtype proportions ``family × conditional``,
    5. moves unresolved families' mass into ``unresolved_<family>`` columns.

    When *family_gene_panels* (``{family: [genes]}``) is supplied, within-family
    resolvability and subtype confidence are judged on each family's own gene
    panel rather than the global panel — so families that the global signature
    could not sub-resolve may become resolvable.  Default (``None``) preserves
    the original global-gene behaviour exactly.

    The deconvolution itself is done by the caller (bulk or spatial pipeline);
    this function performs only the (modality-agnostic) hierarchy arithmetic
    and gating, so the core solvers are never modified.
    """
    resolvability = evaluate_within_family_resolvability(
        fine_ref, mapping,
        unresolved_threshold=unresolved_threshold,
        min_discriminating_genes=min_discriminating_genes,
        within_family_spillover_threshold=within_family_spillover_threshold,
        family_gene_panels=family_gene_panels,
    )
    conditional = compute_conditional_subtype_proportions(
        fine_props, family_props, mapping)
    subtypes = [str(c) for c in fine_props.columns]
    group_of = {st: mapping.get(st, st) for st in subtypes}

    # ----- Gating-mode dispatch -------------------------------------------
    # Resolve the effective mode. ``allow_unresolved=False`` forces the
    # diagnostic ungated behaviour (resolve every family, no unresolved mass),
    # regardless of the nominal gating mode, preserving the historic semantics
    # of that flag.
    gating_mode = str(gating).lower()
    if gating_mode not in ("soft", "hard", "ungated"):
        raise ValueError(
            f"hierarchical_gating must be 'soft', 'hard', or 'ungated' "
            f"(got {gating!r})")
    effective_mode = gating_mode if allow_unresolved else "ungated"

    if effective_mode in ("soft", "hard") and allow_partial_resolution:
        # Per-subtype gating: keep (a fraction of) subtype mass, residual →
        # unresolved_<family>.  soft = confidence-WEIGHTED (default, validated
        # on breast + lung); hard = legacy binary threshold gate.
        confidence = compute_within_family_subtype_confidence(
            fine_ref, mapping, min_discriminating_genes=min_discriminating_genes,
            family_gene_panels=family_gene_panels)
        resolved_fine = pd.DataFrame(0.0, index=family_props.index, columns=subtypes)
        residuals = {}
        for fam in family_props.columns:
            members = [st for st in subtypes if group_of[st] == fam]
            if not members:
                continue
            if effective_mode == "soft":
                sub, residual = estimate_soft_subtype_resolution(
                    fam, members, family_props[fam], conditional, confidence)
            else:
                sub, residual = estimate_partial_subtype_resolution(
                    fam, members, family_props[fam], conditional, confidence,
                    subtype_confidence_threshold=subtype_confidence_threshold)
            for m in members:
                resolved_fine[m] = sub[m]
            if float(residual.abs().mean()) > 1e-9:
                residuals[f"unresolved_{fam}"] = residual
        unresolved_mass = (pd.DataFrame(residuals)
                           if residuals else
                           pd.DataFrame(index=family_props.index))
        unresolved_families = sorted(c[len("unresolved_"):] for c in unresolved_mass.columns)
        subtype_confidence = confidence
    elif effective_mode == "ungated":
        # Diagnostic only: resolve every family into subtypes, no abstention.
        unresolved_families = []
        resolved_fine = combine_family_and_conditional_estimates(
            family_props, conditional, mapping)
        unresolved_mass = pd.DataFrame(index=family_props.index)
        subtype_confidence = {}
    else:
        # Legacy family-level hard gate (hard mode with partial resolution off).
        unresolved_families = decide_unresolved_families(
            resolvability, allow_unresolved=allow_unresolved)
        absolute = combine_family_and_conditional_estimates(
            family_props, conditional, mapping)
        resolved_fine, unresolved_mass = add_unresolved_family_mass(
            absolute, family_props, mapping, unresolved_families)
        subtype_confidence = {}

    combined = pd.concat([resolved_fine, unresolved_mass], axis=1)

    # per-family QC table (counts of subtypes, gating metrics, decision)
    qc = resolvability.copy()
    qc.insert(0, "broad_family", qc.index)
    qc["decision"] = [
        "unresolved (report at family level)" if f in unresolved_families
        else "resolved into subtypes"
        for f in qc.index
    ]
    qc = qc.reset_index(drop=True)

    meta = {
        "n_families": int(family_props.shape[1]),
        "n_fine_subtypes": int(fine_props.shape[1]),
        "n_unresolved_families": len(unresolved_families),
        "unresolved_families": list(unresolved_families),
        "allow_unresolved": bool(allow_unresolved),
        "allow_partial_resolution": bool(allow_partial_resolution and allow_unresolved),
        "subtype_confidence_threshold": subtype_confidence_threshold,
        "unresolved_threshold": unresolved_threshold,
        "min_discriminating_genes": min_discriminating_genes,
        "within_family_spillover_threshold": within_family_spillover_threshold,
        "n_confident_subtypes": int(
            sum(1 for v in subtype_confidence.values()
                if v.get("confidence", 0) >= subtype_confidence_threshold))
        if subtype_confidence else None,
    }

    # ----- Gating provenance + mass-conservation diagnostics --------------
    _mass_total = combined.sum(axis=1)
    _fam_total = family_props.sum(axis=1)
    _mass_err = (float((_mass_total - _fam_total).abs().max())
                 if len(_mass_total) else 0.0)
    _grand = float(combined.to_numpy().sum()) or 1.0
    _unres_per_family = {
        c[len("unresolved_"):]: float(unresolved_mass[c].sum())
        for c in unresolved_mass.columns
    }
    meta.update({
        "hierarchical_gating": gating_mode,
        "gating_effective_mode": effective_mode,
        "gating_version": gating_version,
        "gating_default": "soft",
        "gating_is_default": gating_mode == "soft",
        "gating_confidence_model": (
            "within_family_discriminating_gene_confidence"
            if subtype_confidence else None),
        "gating_confidence_features": (
            ["n_discriminating_genes", "within_family_spillover",
             "conditional_separation"] if subtype_confidence else None),
        "mass_conservation_max_error": _mass_err,
        "unresolved_mass_fraction": float(
            unresolved_mass.to_numpy().sum() / _grand)
        if unresolved_mass.shape[1] else 0.0,
        "unresolved_mass_per_family": _unres_per_family,
        "gating_validation_status": (
            "validated on breast and lung benchmarks" if gating_mode == "soft"
            else "legacy" if gating_mode == "hard" else "diagnostic"),
        "hard_gating_status": "legacy",
    })
    # expose per-subtype confidence so the resolution decision layer (and report)
    # can select supported subtypes without recomputation
    if subtype_confidence:
        meta["subtype_confidence"] = {
            str(k): {"confidence": float(v.get("confidence", 0.0))}
            for k, v in subtype_confidence.items()}

    # ----- Resolution Decision Layer: trusted resolution per family --------
    # Decide broad_only / selected_fine / full_fine from the (already computed)
    # within-family signature evidence + per-subtype confidence, BEFORE fine
    # predictions are interpreted.  Recorded in metadata + the qc table so it
    # affects the outputs/report, not only visualisation.  Full-panel reliability
    # dominates; cell-level AUROC is never used here.
    try:
        from tissueresolve.resolution import (
            decide_trusted_resolution, decisions_metadata)
        _decisions = decide_trusted_resolution(
            identifiability_metrics=resolvability, family_map=mapping,
            subtype_confidence=subtype_confidence or None)
        meta["resolution_decision"] = decisions_metadata(_decisions)
        meta["trusted_resolution"] = {f: d.status for f, d in _decisions.items()}
        qc["trusted_resolution"] = [
            _decisions[str(f)].status if str(f) in _decisions else "full_fine"
            for f in qc["broad_family"]]
    except Exception as _exc:  # never let the decision layer break deconvolution
        meta["resolution_decision_error"] = str(_exc)

    if extra_metadata:
        meta.update(extra_metadata)

    return HierarchicalEstimates(
        family_proportions=family_props,
        conditional_proportions=conditional,
        fine_proportions=resolved_fine,
        unresolved_mass=unresolved_mass,
        combined_fine=combined,
        resolvability=resolvability,
        mapping=dict(mapping),
        unresolved_families=list(unresolved_families),
        qc=qc,
        metadata=meta,
    )
