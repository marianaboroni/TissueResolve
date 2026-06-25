"""Experimental donor-stable gene weighting for bulk deconvolution (opt-in).

Computes a per-gene weight emphasising genes that are cell-type informative,
stable across donors, detected in the query, protocol-compatible, and
non-redundant — inspired by multi-subject tools (e.g. MuSiC) but implemented
independently. It is a **pure scoring function**: it does not replace the bulk
solver, change default weighting, or remove genes (every selected gene receives
a finite, non-negative weight).

    weight = specificity × donor_stability × query_detection
           × protocol_compatibility × non_redundancy

Exposed for experimentation via ``bulk_gene_weighting = donor_stable_experimental``
(default remains the existing weighting). No subtype-contrast amplification.
"""
from __future__ import annotations

from typing import Any, Optional

import numpy as np
import pandas as pd

FEATURE_STATUS = "experimental"
ALGORITHM_VERSION = "bulk_donor_stable_weighting-0.1.0"


def _dense(X):
    return np.asarray(X.toarray() if hasattr(X, "toarray") else X, dtype=float)


def _per_type_means(X, cell_types):
    types = sorted(set(cell_types))
    ct = np.asarray(cell_types)
    M = np.vstack([X[ct == t].mean(0) if (ct == t).any() else np.zeros(X.shape[1])
                   for t in types])     # (T, G)
    return types, M


def compute_donor_stable_gene_weights(
    reference_adata,
    cell_type_col: str,
    donor_col: Optional[str],
    selected_genes,
    query_counts=None,
    family_map: Optional[dict] = None,
    protocol_blacklist: Optional[set] = None,
    eps: float = 1e-3,
) -> pd.Series:
    """Return a pandas Series of per-gene weights in (0, 1], indexed by
    ``selected_genes`` (order preserved; no gene removed).

    Falls back safely when ``donor_col`` is absent (donor_stability = 1).
    """
    genes = [str(g) for g in selected_genes]
    var_names = [str(v) for v in reference_adata.var_names]
    gidx = {g: i for i, g in enumerate(var_names)}
    present = [g for g in genes if g in gidx]
    cols = [gidx[g] for g in present]

    X = _dense(reference_adata.X[:, cols]) if cols else np.zeros((reference_adata.n_obs, 0))
    cell_types = reference_adata.obs[cell_type_col].astype(str).to_numpy()
    types, M = _per_type_means(X, cell_types)            # (T, len(present))

    # --- specificity: concentration of a gene's per-type mean profile (1 - normalized entropy)
    Mpos = np.clip(M, 0, None)
    colsum = Mpos.sum(0, keepdims=True)
    p = np.divide(Mpos, colsum, out=np.full_like(Mpos, 1.0 / max(len(types), 1)),
                  where=colsum > 0)
    ent = -(p * np.log(p + 1e-12)).sum(0)
    max_ent = np.log(max(len(types), 2))
    specificity = 1.0 - ent / max_ent                    # 1 = type-specific, 0 = uniform

    # --- donor stability: 1 / (1 + mean within-type donor CV); fallback 1 if no donor col
    if donor_col and donor_col in reference_adata.obs.columns:
        donors = reference_adata.obs[donor_col].astype(str).to_numpy()
        stability = np.zeros(len(present))
        for gi in range(len(present)):
            cvs = []
            for t in types:
                mask_t = cell_types == t
                if not mask_t.any():
                    continue
                dvals = []
                for d in set(donors[mask_t]):
                    v = X[mask_t & (donors == d), gi]
                    if v.size:
                        dvals.append(v.mean())
                if len(dvals) >= 2:
                    mu = np.mean(dvals)
                    cvs.append((np.std(dvals) / mu) if mu > 1e-9 else 0.0)
            stability[gi] = 1.0 / (1.0 + (np.mean(cvs) if cvs else 0.0))
        donor_fallback = False
    else:
        stability = np.ones(len(present))
        donor_fallback = True

    # --- query detection: fraction of query samples where gene > 0 (else 1)
    if query_counts is not None:
        q = query_counts
        if isinstance(q, pd.DataFrame):
            qd = {str(g): float((q.loc[g] > 0).mean()) if g in q.index else 0.0 for g in present}
        else:
            qd = {g: 1.0 for g in present}
        query_detection = np.array([max(qd.get(g, 0.0), eps) for g in present])
    else:
        query_detection = np.ones(len(present))

    # --- protocol compatibility: downweight blacklisted (stress/dissociation/length-biased)
    bl = {str(x) for x in (protocol_blacklist or set())}
    protocol = np.array([0.2 if g in bl else 1.0 for g in present])

    # --- non-redundancy: 1 - max |corr| of per-type profile with another gene
    if M.shape[1] >= 2 and M.shape[0] >= 2:
        C = np.corrcoef(M.T)
        np.fill_diagonal(C, 0.0)
        max_abs = np.nanmax(np.abs(C), axis=1)
        max_abs = np.where(np.isfinite(max_abs), max_abs, 0.0)
        non_redundancy = 1.0 - 0.5 * max_abs            # mild penalty; never zero
    else:
        non_redundancy = np.ones(len(present))

    w_present = (np.clip(specificity, eps, 1.0)
                 * np.clip(stability, eps, 1.0)
                 * np.clip(query_detection, eps, 1.0)
                 * np.clip(protocol, eps, 1.0)
                 * np.clip(non_redundancy, eps, 1.0))
    w_present = np.clip(w_present, eps, None)
    w_present = w_present / w_present.max() if w_present.size and w_present.max() > 0 else w_present

    # assemble full series (missing genes get neutral weight; never removed)
    weights = pd.Series(1.0, index=genes, dtype=float)
    for g, w in zip(present, w_present):
        weights[g] = float(w)
    weights.attrs["donor_fallback"] = donor_fallback
    weights.attrs["n_present"] = len(present)
    weights.attrs["n_missing"] = len(genes) - len(present)
    weights.attrs["algorithm_version"] = ALGORITHM_VERSION
    return weights
