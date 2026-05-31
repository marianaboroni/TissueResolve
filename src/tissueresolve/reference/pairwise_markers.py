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
