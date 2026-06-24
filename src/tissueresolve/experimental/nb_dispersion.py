"""Donor-level gene-wise NB dispersion estimation (P3; opt-in).

The NB bulk GLM (P1) reduced to Poisson because the reference's stored ``phi_g`` made
the NB term near-inert. This module estimates a **gene-wise NB dispersion from
donor-level pseudobulk** (biological replicate variability), which is the dispersion
relevant to between-sample bulk variation, so the NB likelihood becomes informative
and can be tested against Poisson.

Parameterisation matches the solver: ``Var(y) = μ + μ²/φ_g`` (so ``φ → ∞`` is
Poisson; small ``φ`` is highly overdispersed). Estimated by method of moments on
per-donor expression *rates* (counts / library size), removing the Poisson sampling
contribution so only biological over-dispersion remains:

    α_g = max( CV²_rate,g − mean_d 1/(ℓ_d · m_g), 0 ),   φ_g = 1/α_g

with ``φ_g`` clipped to ``[phi_min, phi_max]``. Falls back to Poisson (``φ = phi_max``)
for genes with < 2 donors or non-estimable variance. Diagnostic / opt-in; does not
change defaults.
"""
from __future__ import annotations

from typing import Optional, Sequence

import numpy as np
import pandas as pd

__all__ = ["estimate_gene_dispersion_from_donors", "FEATURE_STATUS"]

FEATURE_STATUS = "experimental"
_EPS = 1e-12


def estimate_gene_dispersion_from_donors(
    reference_adata,
    cell_type_col: Optional[str] = None,
    donor_col: str = "donor",
    genes: Optional[Sequence[str]] = None,
    *,
    cell_type: Optional[str] = None,
    phi_min: float = 0.1,
    phi_max: float = 1e4,
) -> pd.Series:
    """Gene-wise NB dispersion ``φ_g`` from donor-level pseudobulk.

    Parameters
    ----------
    reference_adata:
        AnnData with raw counts in ``X`` and ``obs[donor_col]``.
    cell_type_col, cell_type:
        If both given, restrict to cells of ``cell_type`` (per-type dispersion);
        otherwise pool all cells (global gene dispersion).
    genes:
        Subset / order of genes to return (default: all ``var_names``).

    Returns
    -------
    pd.Series indexed by gene → ``φ_g`` (NB dispersion). High = ~Poisson.
    """
    import scipy.sparse as sp
    obs = reference_adata.obs
    if donor_col not in obs.columns:
        # no donors → cannot estimate biological dispersion → Poisson-equivalent
        gl = list(map(str, genes)) if genes is not None else list(map(str, reference_adata.var_names))
        return pd.Series(np.full(len(gl), phi_max), index=gl, name="phi_g")

    ad = reference_adata
    if cell_type_col and cell_type is not None and cell_type_col in obs.columns:
        ad = ad[obs[cell_type_col].astype(str) == str(cell_type)]
        obs = ad.obs
    if genes is not None:
        keep = [g for g in genes if g in set(map(str, ad.var_names))]
        ad = ad[:, keep]
    gene_list = list(map(str, ad.var_names))

    X = ad.X
    X = X.toarray() if sp.issparse(X) else np.asarray(X)
    X = np.asarray(X, dtype=np.float64)                  # (cells, genes)
    donors = obs[donor_col].astype(str).to_numpy()
    uniq = sorted(set(donors))
    if len(uniq) < 2:
        return pd.Series(np.full(len(gene_list), phi_max), index=gene_list, name="phi_g")

    # per-donor pseudobulk total counts + library size → expression rate
    rates, libs = [], []
    for d in uniq:
        m = donors == d
        tot = X[m].sum(axis=0)                           # (genes,) summed counts
        lib = float(tot.sum())
        if lib <= 0:
            continue
        rates.append(tot / lib)                          # rate per gene (sums to 1)
        libs.append(lib)
    if len(rates) < 2:
        return pd.Series(np.full(len(gene_list), phi_max), index=gene_list, name="phi_g")
    R = np.vstack(rates)                                 # (D, genes)
    libs = np.asarray(libs)
    m_g = R.mean(axis=0)                                 # mean rate per gene
    v_g = R.var(axis=0, ddof=1)                          # between-donor variance of rate
    # Poisson sampling contribution to rate variance ≈ mean_d m_g/ℓ_d  (Var(count)=μ)
    pois_cv2 = np.mean(1.0 / libs) / np.maximum(m_g, _EPS)
    cv2 = v_g / np.maximum(m_g ** 2, _EPS)
    alpha = np.maximum(cv2 - pois_cv2, 0.0)              # biological over-dispersion
    phi = np.where(alpha > _EPS, 1.0 / alpha, phi_max)
    phi = np.clip(phi, phi_min, phi_max)
    # genes with ~zero expression → uninformative → Poisson
    phi[m_g <= _EPS] = phi_max
    return pd.Series(phi, index=gene_list, name="phi_g")
