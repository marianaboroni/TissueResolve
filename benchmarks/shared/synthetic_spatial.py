"""Structured synthetic spatial generator with ground truth (PART 2B.1).

Builds a square grid of spots with real spatial structure — a **sharp vertical
domain border**, a **smooth top→bottom gradient**, and a **rare-cell niche** —
then realises count-level spots by sampling **held-out (query) donor** cells per
spot (reusing the bulk generator's `realize_pseudobulk`), so the cells forming
the spots are disjoint from the reference donors.  Produces per-spot broad+fine
mRNA-proportion truth, spot coordinates, library sizes, and domain labels.

Returns a scenario dict compatible with `tissueresolve.deconv_spatial`:
``Y`` (spots×genes), ``array_row``, ``array_col``, ``lib_sizes``, ``gene_names``,
``spot_ids``, ``reference`` (filled by caller) plus ``truth`` (spots×types) and
``domain_labels``.
"""
from __future__ import annotations

import numpy as np
import pandas as pd

from benchmarks.shared import synthetic_holdout as SH


def _family_composition(ref_types, mapping, families, rng, *, fg_mass=0.6):
    """Composition vector concentrating *fg_mass* on members of *families*."""
    types = list(ref_types)
    k = len(types)
    members = [t for t in types if mapping.get(t, t) in set(families)]
    if not members:
        members = [rng.choice(types)]
    comp = np.full(k, (1.0 - fg_mass) / k)            # uniform background
    w = rng.dirichlet(np.ones(len(members)))
    for t, wt in zip(members, w):
        comp[types.index(t)] += fg_mass * wt
    return comp / comp.sum()


# Default domain→broad-family assignments (breast); overridable for other tissues.
BREAST_DOMAIN_FAMILIES = {
    "left": ["Epithelial"],
    "right": ["T/NK", "Myeloid"],
    "gradient": ["Stromal/Fibroblast", "Endothelial", "Mural"],
}


def build_spatial_targets(ref_types, mapping, *, n_side=20, seed=0,
                          rare_type=None, niche_size=3, noise=0.04,
                          domain_families=None):
    """Per-spot target proportions with sharp border + gradient + rare niche.

    Domains: left half = L (``domain_families['left']``-rich), right half = R
    (``domain_families['right']``-rich), sharp vertical border at the midline.
    A smooth top→bottom gradient blends each domain toward
    ``domain_families['gradient']``.  A *niche_size* square block (top-right) is
    enriched for *rare_type*.  ``domain_families`` defaults to the breast layout
    (:data:`BREAST_DOMAIN_FAMILIES`); pass tissue-specific broad-family names for
    other datasets.
    """
    rng = np.random.default_rng(seed)
    types = list(ref_types)
    k = len(types)
    fams = domain_families or BREAST_DOMAIN_FAMILIES
    cL = _family_composition(types, mapping, fams["left"], rng)
    cR = _family_composition(types, mapping, fams["right"], rng)
    cG = _family_composition(types, mapping, fams["gradient"], rng)
    mid = n_side // 2
    rare_type = rare_type or types[-1]
    ri = types.index(rare_type)

    rows, cols, P, domains = [], [], [], []
    niche_r0, niche_c0 = 1, n_side - niche_size - 1     # top-right block
    for r in range(n_side):
        for c in range(n_side):
            base = cL if c < mid else cR
            dom = "L" if c < mid else "R"
            g = r / max(n_side - 1, 1)                  # gradient weight along rows
            comp = (1.0 - g) * base + g * cG
            in_niche = (niche_r0 <= r < niche_r0 + niche_size and
                        niche_c0 <= c < niche_c0 + niche_size)
            if in_niche:
                comp = comp.copy(); comp[ri] += 0.10; dom = "niche"
            comp = comp + rng.dirichlet(np.ones(k)) * noise
            comp = np.clip(comp, 0, None); comp /= comp.sum()
            rows.append(r); cols.append(c); P.append(comp); domains.append(dom)
    spot_ids = [f"spot{i}" for i in range(len(rows))]
    targets = pd.DataFrame(np.array(P), index=spot_ids, columns=types)
    coords = np.c_[rows, cols].astype(float)
    return targets, coords, np.array(domains), {"rare_type": rare_type,
                                                "n_side": n_side, "niche_size": niche_size}


def generate_spatial_scenario(adata, ref_types, query_donors, mapping, *,
                              celltype_col="cell_type", donor_col="donor_id",
                              n_side=20, cells_per_spot=40, seed=0, rare_type=None,
                              domain_families=None):
    """Full structured spatial scenario realised from held-out query-donor cells."""
    targets, coords, domains, meta = build_spatial_targets(
        ref_types, mapping, n_side=n_side, seed=seed, rare_type=rare_type,
        domain_families=domain_families)
    ds = SH.realize_pseudobulk(adata, targets, celltype_col=celltype_col,
                               donor_col=donor_col, query_donors=query_donors,
                               seed=seed * 13 + 3, cells_per_sample=cells_per_spot)
    counts = ds.counts                                  # genes × spots
    truth = ds.true_mrna_proportions                    # spots × types
    Y = counts.to_numpy(float).T                        # spots × genes
    return {
        "Y": Y, "array_row": coords[:, 0].astype(int), "array_col": coords[:, 1].astype(int),
        "lib_sizes": counts.sum(axis=0).to_numpy(float),
        "gene_names": list(map(str, counts.index)), "spot_ids": list(counts.columns),
        "truth": truth, "coords": coords, "domain_labels": domains,
        "hierarchy_mapping": mapping, "meta": meta,
    }
