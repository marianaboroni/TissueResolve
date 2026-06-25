"""Tests for donor-level NB dispersion estimation (P3)."""
from __future__ import annotations

import sys
from pathlib import Path

import numpy as np
import pandas as pd
import pytest

REPO = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(REPO / "src"))

from tissueresolve.experimental.nb_dispersion import estimate_gene_dispersion_from_donors


def _adata(n_donors=8, cells_per_donor=40, seed=0):
    """Gene 0 = donor-stable; gene 1 = strongly donor-variable; rest random."""
    import anndata as ad
    rng = np.random.default_rng(seed)
    G = 12
    rows, donors = [], []
    stable_rate = 0.2
    for d in range(n_donors):
        var_rate = 0.2 * (1.0 + rng.normal(0, 0.8))      # big between-donor swing
        var_rate = max(var_rate, 0.01)
        for _ in range(cells_per_donor):
            x = rng.poisson(np.r_[stable_rate * 100, var_rate * 100,
                                  rng.random(G - 2) * 20])
            rows.append(x)
            donors.append(f"D{d}")
    X = np.array(rows, float)
    obs = pd.DataFrame({"donor": donors, "cell_type": ["T"] * len(donors)})
    var = pd.DataFrame(index=[f"g{i}" for i in range(G)])
    return ad.AnnData(X=X, obs=obs, var=var)


def test_returns_phi_per_gene_finite_in_range():
    a = _adata()
    phi = estimate_gene_dispersion_from_donors(a, donor_col="donor")
    assert len(phi) == a.n_vars
    assert np.all(np.isfinite(phi.to_numpy()))
    assert (phi >= 0.1 - 1e-9).all() and (phi <= 1e4 + 1).all()


def test_donor_variable_gene_has_smaller_phi_than_stable():
    a = _adata()
    phi = estimate_gene_dispersion_from_donors(a, donor_col="donor")
    # g1 is engineered donor-variable (overdispersed) → smaller φ than stable g0
    assert phi["g1"] < phi["g0"]


def test_single_donor_falls_back_to_poisson():
    import anndata as ad
    rng = np.random.default_rng(1)
    X = rng.poisson(5, size=(50, 8)).astype(float)
    obs = pd.DataFrame({"donor": ["D0"] * 50})
    a = ad.AnnData(X=X, obs=obs, var=pd.DataFrame(index=[f"g{i}" for i in range(8)]))
    phi = estimate_gene_dispersion_from_donors(a, donor_col="donor", phi_max=1e4)
    assert np.allclose(phi.to_numpy(), 1e4)


def test_missing_donor_col_falls_back_to_poisson():
    a = _adata()
    phi = estimate_gene_dispersion_from_donors(a, donor_col="nope", phi_max=1e4)
    assert np.allclose(phi.to_numpy(), 1e4)


def test_gene_subset_and_order_respected():
    a = _adata()
    phi = estimate_gene_dispersion_from_donors(a, donor_col="donor", genes=["g2", "g0"])
    assert list(phi.index) == ["g2", "g0"]
