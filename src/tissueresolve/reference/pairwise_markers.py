"""
Pairwise marker refinement for confusable cell-type pairs.

For each high-risk (poorly separable) pair, select genes that best distinguish
exactly those two types, scoring them on:

* **pairwise log2 fold-change** between the two types,
* **donor stability** (low cross-donor CV, when available),
* **spatial detectability** (mean CPM — faint genes are unreliable in Visium),
* **low leakage** into the *other* cell types (a discriminative gene should be
  quiet outside the pair).

Then augment a base marker panel with the top genes for each confusable pair.
This sharpens estimates where they are weakest without altering the core
solver.
"""
from __future__ import annotations

from typing import Optional

import numpy as np
import pandas as pd

from tissueresolve.results import ReferenceSignature

__all__ = [
    "select_pairwise_discriminative_genes",
    "score_pairwise_markers",
    "augment_marker_panel_for_confusable_pairs",
    "score_genes_for_subtype_resolution",
    "select_within_family_discriminative_genes",
    "build_family_specific_gene_panels",
]


def _type_index(ref: ReferenceSignature, name: str) -> int:
    try:
        return list(ref.cell_types).index(name)
    except ValueError as exc:
        raise KeyError(f"cell type {name!r} not in reference: {list(ref.cell_types)}") \
            from exc


def score_pairwise_markers(
    ref: ReferenceSignature, type_a: str, type_b: str,
    genes: Optional[list[str]] = None,
) -> pd.DataFrame:
    """Score genes for discriminating *type_a* vs *type_b*.

    Returns a DataFrame indexed by gene with columns ``log2fc`` (signed,
    a−b), ``abs_log2fc``, ``donor_stability`` (∈[0,1]; 1 if no donor CV),
    ``detectability`` (max CPM in the pair), ``leakage`` (max CPM in the
    *other* types, normalised), and a composite ``score``.
    """
    ia, ib = _type_index(ref, type_a), _type_index(ref, type_b)
    R_log = ref.as_R_log().astype(np.float64)     # (K, G)
    R_cpm = ref.as_R_cpm().astype(np.float64)     # (K, G)
    gene_names = [str(g) for g in ref.gene_names]
    gidx = (np.arange(len(gene_names)) if genes is None
            else np.array([gene_names.index(str(g)) for g in genes
                           if str(g) in gene_names], dtype=int))

    log2fc = R_log[ia, gidx] - R_log[ib, gidx]
    abs_fc = np.abs(log2fc)
    detect = np.maximum(R_cpm[ia, gidx], R_cpm[ib, gidx])

    # leakage: highest expression among all OTHER cell types
    other = [k for k in range(ref.n_cell_types) if k not in (ia, ib)]
    if other:
        other_max = R_cpm[np.ix_(other, gidx)].max(axis=0)
    else:
        other_max = np.zeros(len(gidx))
    pair_max = np.maximum(detect, 1e-9)
    leakage = np.clip(other_max / pair_max, 0.0, 1.0)

    # donor stability from cross-donor CV (lower CV → more stable)
    if ref.donor_cv is not None:
        cv = ref.donor_cv[gidx][:, [ia, ib]].mean(axis=1)
        donor_stability = 1.0 / (1.0 + np.asarray(cv, dtype=float))
    else:
        donor_stability = np.ones(len(gidx))

    detect_norm = detect / (detect.max() + 1e-9)
    score = abs_fc * donor_stability * (1.0 - leakage) * (0.5 + 0.5 * detect_norm)

    return pd.DataFrame({
        "log2fc": log2fc,
        "abs_log2fc": abs_fc,
        "donor_stability": donor_stability,
        "detectability": detect,
        "leakage": leakage,
        "score": score,
    }, index=[gene_names[i] for i in gidx]).sort_values("score", ascending=False)


def select_pairwise_discriminative_genes(
    ref: ReferenceSignature, type_a: str, type_b: str,
    *, top_n: int = 20, min_log2fc: float = 1.0,
    candidate_genes: Optional[list[str]] = None,
) -> list[str]:
    """Top *top_n* genes (by composite score) with ``|log2FC| ≥ min_log2fc``."""
    scored = score_pairwise_markers(ref, type_a, type_b, genes=candidate_genes)
    eligible = scored[scored["abs_log2fc"] >= min_log2fc]
    if eligible.empty:                      # never silently return nothing useful
        eligible = scored
    return eligible.head(top_n).index.tolist()


def augment_marker_panel_for_confusable_pairs(
    ref: ReferenceSignature,
    base_panel: list[str],
    pairs: list[tuple[str, str]],
    *, top_n: int = 10, min_log2fc: float = 1.0,
    candidate_genes: Optional[list[str]] = None,
) -> tuple[list[str], dict[tuple[str, str], list[str]]]:
    """Extend *base_panel* with discriminative genes for each confusable pair.

    Returns ``(augmented_panel, added)`` where *added* maps each pair to the
    genes it contributed.  Genes already in the panel are not duplicated; the
    base panel order is preserved and new genes are appended.
    """
    panel = list(dict.fromkeys(str(g) for g in base_panel))
    panel_set = set(panel)
    added: dict[tuple[str, str], list[str]] = {}
    for a, b in pairs:
        genes = select_pairwise_discriminative_genes(
            ref, a, b, top_n=top_n, min_log2fc=min_log2fc,
            candidate_genes=candidate_genes)
        new = [g for g in genes if g not in panel_set]
        added[(a, b)] = new
        for g in new:
            panel.append(g)
            panel_set.add(g)
    return panel, added


# ---------------------------------------------------------------------------
# Within-family (hierarchical) marker refinement
# ---------------------------------------------------------------------------


def _family_members(mapping: dict[str, str]) -> dict[str, list[str]]:
    fams: dict[str, list[str]] = {}
    for fine, broad in mapping.items():
        fams.setdefault(str(broad), []).append(str(fine))
    return fams


def score_genes_for_subtype_resolution(
    ref: ReferenceSignature,
    family: str,
    family_mapping: dict[str, str],
    *,
    candidate_genes: Optional[list[str]] = None,
) -> pd.DataFrame:
    """Score genes for resolving the fine subtypes *within* one broad *family*.

    For each ordered pair of family members the discriminative score from
    :func:`score_pairwise_markers` is computed; the per-gene resolution score is
    the **maximum** pairwise score across all within-family pairs (a gene that
    cleanly splits any sibling pair is valuable).  Leakage is measured against
    *all* other cell types in the reference, so genes that also fire outside the
    family are penalised — exactly the within-family specificity we want.

    Returns a DataFrame indexed by gene with the best-pair components and a
    composite ``resolution_score`` (sorted descending).  Returns an empty frame
    for families with fewer than two members.
    """
    members = [m for m in _family_members(family_mapping).get(family, [])
               if m in ref.cell_types]
    if len(members) < 2:
        return pd.DataFrame(
            columns=["resolution_score", "best_pair", "log2fc", "abs_log2fc",
                     "donor_stability", "detectability", "leakage"]
        )
    best: dict[str, dict] = {}
    for i in range(len(members)):
        for j in range(i + 1, len(members)):
            a, b = members[i], members[j]
            scored = score_pairwise_markers(ref, a, b, genes=candidate_genes)
            for gene, row in scored.iterrows():
                cur = best.get(gene)
                if cur is None or row["score"] > cur["resolution_score"]:
                    best[gene] = {
                        "resolution_score": float(row["score"]),
                        "best_pair": f"{a} vs {b}",
                        "log2fc": float(row["log2fc"]),
                        "abs_log2fc": float(row["abs_log2fc"]),
                        "donor_stability": float(row["donor_stability"]),
                        "detectability": float(row["detectability"]),
                        "leakage": float(row["leakage"]),
                    }
    df = pd.DataFrame(best).T
    if df.empty:
        return df
    return df.sort_values("resolution_score", ascending=False)


def select_within_family_discriminative_genes(
    ref: ReferenceSignature,
    family_mapping: dict[str, str],
    *,
    top_n: int = 20,
    min_log2fc: float = 1.0,
    candidate_genes: Optional[list[str]] = None,
) -> dict[str, list[str]]:
    """Top discriminative genes for each broad family's fine subtypes.

    Returns ``{family: [genes]}`` for every multi-member family.  Genes are
    ranked by within-family resolution score (see
    :func:`score_genes_for_subtype_resolution`); only genes with
    ``|log2FC| ≥ min_log2fc`` on their best pair are kept, falling back to the
    top-scoring genes if none clear the threshold (never silently empty).
    """
    out: dict[str, list[str]] = {}
    for fam, members in _family_members(family_mapping).items():
        if len([m for m in members if m in ref.cell_types]) < 2:
            continue
        scored = score_genes_for_subtype_resolution(
            ref, fam, family_mapping, candidate_genes=candidate_genes)
        if scored.empty:
            continue
        eligible = scored[scored["abs_log2fc"] >= min_log2fc]
        if eligible.empty:
            eligible = scored
        out[fam] = eligible.head(top_n).index.tolist()
    return out


def build_family_specific_gene_panels(
    ref: ReferenceSignature,
    family_mapping: dict[str, str],
    *,
    base_panel: Optional[list[str]] = None,
    top_n: int = 20,
    min_log2fc: float = 1.0,
    candidate_genes: Optional[list[str]] = None,
) -> tuple[dict[str, list[str]], pd.DataFrame]:
    """Build a discriminative gene panel for each broad family.

    Each family's panel is the (optional) *base_panel* augmented with the top
    within-family discriminative genes for that family's subtypes.  Genes are
    restricted to those present in the reference.

    Returns
    -------
    (panels, discriminability)
        *panels* maps family → ordered gene list (base genes first, then new
        within-family genes).  *discriminability* is a tidy DataFrame with
        columns ``family, gene, resolution_score, best_pair, abs_log2fc,
        leakage`` describing the selected within-family genes (for the
        ``within_family_discriminability.tsv`` output).
    """
    ref_genes = set(str(g) for g in ref.gene_names)
    base = [g for g in (base_panel or []) if str(g) in ref_genes]
    per_family = select_within_family_discriminative_genes(
        ref, family_mapping, top_n=top_n, min_log2fc=min_log2fc,
        candidate_genes=candidate_genes)

    panels: dict[str, list[str]] = {}
    disc_rows: list[dict] = []
    for fam, members in _family_members(family_mapping).items():
        present = [m for m in members if str(m) in ref_genes or m in ref.cell_types]
        if len([m for m in members if m in ref.cell_types]) < 2:
            continue  # singleton family: no within-family discrimination
        fam_genes = per_family.get(fam, [])
        panel = list(dict.fromkeys(base + fam_genes))
        panels[fam] = panel
        if fam_genes:
            scored = score_genes_for_subtype_resolution(ref, fam, family_mapping,
                                                        candidate_genes=candidate_genes)
            for g in fam_genes:
                if g in scored.index:
                    r = scored.loc[g]
                    disc_rows.append({
                        "family": fam, "gene": g,
                        "resolution_score": round(float(r["resolution_score"]), 4),
                        "best_pair": r["best_pair"],
                        "abs_log2fc": round(float(r["abs_log2fc"]), 4),
                        "leakage": round(float(r["leakage"]), 4),
                    })
    discriminability = pd.DataFrame(
        disc_rows,
        columns=["family", "gene", "resolution_score", "best_pair",
                 "abs_log2fc", "leakage"],
    )
    return panels, discriminability
