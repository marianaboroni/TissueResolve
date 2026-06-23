"""Held-out-donor synthetic bulk pseudobulk generator (PART 2A.1).

Generates **count-level** pseudobulk mixtures with known ground truth for bulk
deconvolution benchmarking, with the critical property that the donors whose
cells form the mixtures are **disjoint** from the donors used to build the
reference signature — so accuracy is not inflated by reusing the same cells on
both sides ("Do not ... use different donors" / cross-donor requirement).

Design
------
1. ``split_donors`` partitions donors into a *reference* set and a disjoint
   *query* set (optionally honouring an assay/platform split for cross-platform
   scenarios).
2. ``build_target_proportions`` produces per-sample target cell-type proportions
   for a named scenario (balanced / imbalanced / rare / similar / missing / extra).
3. ``realize_pseudobulk`` samples cells from the *query* donors according to each
   target, sums their **raw integer counts**, and records two ground truths:

   * ``true_cell_fractions``  — # cells of a type / total cells (what FACS sees);
   * ``true_mrna_proportions`` — total counts from a type / total counts (what an
     RNA-based bulk deconvolver actually estimates; see CLAUDE.md bulk rule 1).

Everything is deterministic given ``seed``.  No network, no global RNG.
"""
from __future__ import annotations

from dataclasses import dataclass, field
from typing import Optional

import numpy as np
import pandas as pd

SCENARIOS = ("balanced", "imbalanced", "rare", "similar_subtypes",
             "missing_population", "extra_population")


@dataclass
class PseudobulkDataset:
    """One generated scenario."""
    counts: pd.DataFrame                 # genes × samples, integer counts
    true_cell_fractions: pd.DataFrame    # samples × cell_types, rows sum to 1
    true_mrna_proportions: pd.DataFrame  # samples × cell_types, rows sum to 1
    metadata: pd.DataFrame               # samples × params
    reference_donors: list
    query_donors: list
    scenario: str
    held_out_types: list = field(default_factory=list)  # types absent from reference
    params: dict = field(default_factory=dict)


# --------------------------------------------------------------------------- donors
def split_donors(adata, donor_col: str, *, ref_frac: float = 0.5,
                 seed: int = 0, assay_col: Optional[str] = None,
                 cross_platform: bool = False) -> tuple[list, list]:
    """Partition donors into disjoint (reference, query) lists.

    When *cross_platform* and *assay_col* are given, reference donors are drawn
    from the majority assay and query donors from the other assay (a
    cross-platform reference/query split), falling back to a random split if a
    clean assay separation is not possible.
    """
    rng = np.random.default_rng(seed)
    obs = adata.obs
    donors = sorted(map(str, obs[donor_col].unique()))

    if cross_platform and assay_col is not None and assay_col in obs.columns:
        # donor -> dominant assay
        dom = (obs.groupby(donor_col)[assay_col]
               .agg(lambda s: s.value_counts().index[0]).astype(str))
        assays = dom.value_counts().index.tolist()
        if len(assays) >= 2:
            ref = sorted(str(d) for d in dom.index[dom == assays[0]])
            qry = sorted(str(d) for d in dom.index[dom == assays[1]])
            if ref and qry:
                return ref, qry

    perm = rng.permutation(donors)
    n_ref = max(1, int(round(len(donors) * ref_frac)))
    ref = sorted(perm[:n_ref].tolist())
    qry = sorted(perm[n_ref:].tolist())
    if not qry:  # degenerate tiny inputs
        qry = [ref.pop()]
    return ref, qry


# ----------------------------------------------------------------- target proportions
def build_target_proportions(cell_types, n_samples: int, scenario: str, *,
                             seed: int = 0, rare_type: Optional[str] = None,
                             rare_level: float = 0.01,
                             similar_pair: Optional[tuple] = None,
                             missing_type: Optional[str] = None,
                             extra_type: Optional[str] = None,
                             imbalance_alpha: float = 0.3) -> pd.DataFrame:
    """Per-sample target cell-type proportions for *scenario* (rows sum to 1).

    ``missing_population`` / ``extra_population`` shape the *proportions* only;
    which types are withheld from the **reference** is recorded separately by the
    caller (``held_out_types``) so the mixture genuinely contains a population the
    reference lacks.
    """
    rng = np.random.default_rng(seed)
    types = [str(c) for c in cell_types]
    k = len(types)
    idx = {t: i for i, t in enumerate(types)}
    P = np.zeros((n_samples, k))

    if scenario == "balanced":
        P[:] = 1.0 / k
    elif scenario == "imbalanced":
        P[:] = rng.dirichlet(np.full(k, imbalance_alpha), size=n_samples)
    elif scenario == "rare":
        rt = rare_type or types[-1]
        others = [t for t in types if t != rt]
        base = rng.dirichlet(np.full(len(others), 1.0), size=n_samples)
        for s in range(n_samples):
            P[s, idx[rt]] = rare_level
            P[s, [idx[t] for t in others]] = base[s] * (1.0 - rare_level)
    elif scenario == "similar_subtypes":
        a, b = similar_pair if similar_pair else (types[0], types[1])
        for s in range(n_samples):
            share = rng.uniform(0.5, 0.8)  # most mass on the similar pair
            split = rng.uniform(0.3, 0.7)
            P[s, idx[a]] = share * split
            P[s, idx[b]] = share * (1.0 - split)
            rest = [t for t in types if t not in (a, b)]
            if rest:
                r = rng.dirichlet(np.ones(len(rest)))
                P[s, [idx[t] for t in rest]] = r * (1.0 - share)
    elif scenario in ("missing_population", "extra_population"):
        # full Dirichlet over all present types; the held-out one is present in
        # the mixture (and recorded by the caller as absent from the reference).
        special = missing_type or extra_type or types[-1]
        P[:] = rng.dirichlet(np.full(k, 1.0), size=n_samples)
        # ensure the special type carries non-trivial mass so the scenario bites
        for s in range(n_samples):
            P[s, idx[special]] = max(P[s, idx[special]], rng.uniform(0.05, 0.15))
        P = P / P.sum(axis=1, keepdims=True)
    else:
        raise ValueError(f"unknown scenario {scenario!r}; choose from {SCENARIOS}")

    P = np.clip(P, 0.0, None)
    P = P / P.sum(axis=1, keepdims=True)
    samples = [f"{scenario}_s{i:03d}" for i in range(n_samples)]
    return pd.DataFrame(P, index=samples, columns=types)


# ------------------------------------------------------------------- realize counts
def _row_counts(adata, rows) -> np.ndarray:
    """Sum raw count rows (sparse or dense) → 1-D gene vector (float)."""
    X = adata.X[rows]
    if hasattr(X, "toarray"):
        return np.asarray(X.sum(axis=0)).ravel().astype(float)
    return np.asarray(X).sum(axis=0).astype(float)


def realize_pseudobulk(adata, targets: pd.DataFrame, *, celltype_col: str,
                       donor_col: str, query_donors, seed: int = 0,
                       cells_per_sample: int = 500,
                       depth_factor: float = 1.0,
                       noise_cv: float = 0.0) -> PseudobulkDataset:
    """Build count-level pseudobulk from *query_donors* per *targets*.

    For each sample, draw ``round(prop_k · cells_per_sample)`` cells of type *k*
    (with replacement) from the query donors, sum their raw counts, optionally
    rescale total depth and add per-gene multiplicative noise (lognormal, given
    *noise_cv*).  Records realised cell-fraction and mRNA-proportion truths.
    """
    rng = np.random.default_rng(seed)
    obs = adata.obs
    genes = list(map(str, adata.var_names))
    types = [str(c) for c in targets.columns]
    qset = set(map(str, query_donors))

    # candidate row indices per cell type, restricted to query donors
    donor_str = obs[donor_col].astype(str).to_numpy()
    ct_str = obs[celltype_col].astype(str).to_numpy()
    in_query = np.isin(donor_str, list(qset))
    rows_of_type = {t: np.where((ct_str == t) & in_query)[0] for t in types}

    counts = np.zeros((len(genes), targets.shape[0]))
    cell_frac = np.zeros_like(targets.to_numpy(float))
    mrna_prop = np.zeros_like(targets.to_numpy(float))
    meta = []

    for s, samp in enumerate(targets.index):
        tgt = targets.iloc[s].to_numpy(float)
        n_k = np.floor(tgt * cells_per_sample + 0.5).astype(int)
        # guarantee at least one cell where target is clearly non-zero
        n_k[(tgt > 0) & (n_k == 0) & (tgt >= 0.5 / cells_per_sample)] = 1
        sample_vec = np.zeros(len(genes))
        type_counts = np.zeros(len(types))
        realized_cells = np.zeros(len(types))
        for ti, t in enumerate(types):
            if n_k[ti] <= 0:
                continue
            pool = rows_of_type[t]
            if pool.size == 0:
                continue  # type absent among query donors — recorded as 0 below
            pick = rng.choice(pool, size=int(n_k[ti]), replace=True)
            cvec = _row_counts(adata, pick)
            sample_vec += cvec
            type_counts[ti] = cvec.sum()
            realized_cells[ti] = n_k[ti]
        if noise_cv > 0:
            sample_vec = sample_vec * rng.lognormal(0.0, noise_cv, size=len(genes))
        sample_vec = np.floor(sample_vec * depth_factor + 0.5)
        counts[:, s] = sample_vec
        cf = realized_cells / realized_cells.sum() if realized_cells.sum() > 0 else realized_cells
        mp = type_counts / type_counts.sum() if type_counts.sum() > 0 else type_counts
        cell_frac[s] = cf
        mrna_prop[s] = mp
        meta.append({"sample": samp, "n_cells": int(realized_cells.sum()),
                     "total_counts": float(sample_vec.sum()),
                     "depth_factor": depth_factor, "noise_cv": noise_cv})

    counts_df = pd.DataFrame(counts.astype(int), index=genes, columns=list(targets.index))
    return PseudobulkDataset(
        counts=counts_df,
        true_cell_fractions=pd.DataFrame(cell_frac, index=targets.index, columns=types),
        true_mrna_proportions=pd.DataFrame(mrna_prop, index=targets.index, columns=types),
        metadata=pd.DataFrame(meta).set_index("sample"),
        reference_donors=[], query_donors=sorted(qset),
        scenario="", params={"cells_per_sample": cells_per_sample,
                             "depth_factor": depth_factor, "noise_cv": noise_cv},
    )


def generate_scenario(adata, *, celltype_col: str, donor_col: str, scenario: str,
                      n_samples: int = 12, seed: int = 0, ref_frac: float = 0.5,
                      cells_per_sample: int = 500, assay_col: Optional[str] = None,
                      cross_platform: bool = False, held_out_types=None,
                      **target_kwargs) -> PseudobulkDataset:
    """End-to-end: donor split → targets → count-level pseudobulk for *scenario*.

    *held_out_types* (for missing/extra scenarios) are types present in the
    mixture but to be excluded from the reference signature downstream; recorded
    on the returned dataset, not removed from the counts.
    """
    ref_donors, query_donors = split_donors(
        adata, donor_col, ref_frac=ref_frac, seed=seed,
        assay_col=assay_col, cross_platform=cross_platform)
    all_types = sorted(map(str, adata.obs[celltype_col].unique()))
    targets = build_target_proportions(all_types, n_samples, scenario,
                                        seed=seed, **target_kwargs)
    ds = realize_pseudobulk(adata, targets, celltype_col=celltype_col,
                            donor_col=donor_col, query_donors=query_donors,
                            seed=seed, cells_per_sample=cells_per_sample)
    ds.reference_donors = ref_donors
    ds.scenario = scenario
    ds.held_out_types = list(held_out_types or [])
    ds.params.update({"n_samples": n_samples, "ref_frac": ref_frac,
                      "cross_platform": cross_platform, "seed": seed,
                      **target_kwargs})
    return ds
