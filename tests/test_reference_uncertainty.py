"""Tests for reference uncertainty / donor-aware reliability (P2c). Diagnostics only."""
from __future__ import annotations

import sys
from pathlib import Path

import numpy as np
import pandas as pd
import pytest

REPO = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(REPO / "src"))

from tissueresolve.experimental.reference_uncertainty import (
    estimate_reference_uncertainty, ReferenceUncertainty)


def _adata():
    """Build a reference AnnData with controlled reliability per state:
       - 'reliable'      : many cells, many donors, stable
       - 'few_cells'     : few cells, several donors
       - 'few_donors'    : many cells, 1 donor
       - 'donor_unstable': many cells/donors but large donor-to-donor shift
    """
    import anndata as ad
    rng = np.random.default_rng(0)
    G = 80
    base = rng.random(G) * 5 + 1
    rows, labels, donors = [], [], []

    def add(state, n, donor, shift=0.0):
        prof = base + shift
        X = rng.poisson(np.clip(prof, 0.1, None) * 20, size=(n, G)).astype(float)
        rows.append(X)
        labels.extend([state] * n)
        donors.extend([donor] * n)

    for d in range(6):                         # reliable: 6 donors × 60 cells, stable
        add("reliable", 60, f"D{d}")
    for d in range(4):                         # few_cells: 4 donors × 3 cells
        add("few_cells", 3, f"D{d}")
    add("few_donors", 300, "D0")               # many cells, 1 donor
    for d in range(6):                         # donor_unstable: big per-donor shift
        add("donor_unstable", 60, f"D{d}", shift=rng.normal(0, 4, G))

    X = np.vstack(rows)
    obs = pd.DataFrame({"cell_type": labels, "donor": donors})
    var = pd.DataFrame(index=[f"g{i}" for i in range(G)])
    return ad.AnnData(X=X, obs=obs, var=var)


def test_returns_per_state_table():
    a = _adata()
    ru = estimate_reference_uncertainty(a, "cell_type", donor_col="donor")
    assert isinstance(ru, ReferenceUncertainty)
    for col in ("n_cells", "n_donors", "donor_variability", "centroid_se",
                "marker_stability", "outlier_donor_score", "reliability"):
        assert col in ru.per_state.columns


def test_few_cells_have_higher_centroid_uncertainty():
    a = _adata()
    ru = estimate_reference_uncertainty(a, "cell_type", donor_col="donor").per_state
    assert ru.loc["few_cells", "centroid_se"] > ru.loc["reliable", "centroid_se"]


def test_few_donors_lower_reliability_than_reliable():
    a = _adata()
    ru = estimate_reference_uncertainty(a, "cell_type", donor_col="donor").per_state
    assert ru.loc["few_donors", "n_donors"] == 1
    assert ru.loc["few_donors", "reliability"] < ru.loc["reliable", "reliability"]


def test_donor_unstable_penalized():
    a = _adata()
    ru = estimate_reference_uncertainty(a, "cell_type", donor_col="donor").per_state
    assert ru.loc["donor_unstable", "donor_variability"] > ru.loc["reliable", "donor_variability"]
    assert ru.loc["donor_unstable", "reliability"] < ru.loc["reliable", "reliability"]


def test_reliable_state_high_score():
    a = _adata()
    ru = estimate_reference_uncertainty(a, "cell_type", donor_col="donor").per_state
    assert ru.loc["reliable", "reliability"] == ru["reliability"].max()


def test_missing_donor_column_falls_back_safely():
    a = _adata()
    ru = estimate_reference_uncertainty(a, "cell_type", donor_col=None)
    assert ru.metadata["donor_fallback"] is True
    ps = ru.per_state
    assert ps["reliability"].between(0, 1).all()
    assert (ps["n_donors"] == 0).all()
    assert ps["reliability"].notna().all()


def test_family_aggregation_and_write(tmp_path):
    a = _adata()
    mapping = {"reliable": "famA", "few_cells": "famA",
               "few_donors": "famB", "donor_unstable": "famB"}
    ru = estimate_reference_uncertainty(a, "cell_type", donor_col="donor", mapping=mapping)
    assert not ru.family.empty
    assert set(ru.family.index) == {"famA", "famB"}
    paths = ru.write(tmp_path)
    assert paths["reference_uncertainty"].exists()
    assert paths["state_reliability"].exists()
    assert paths["family_reliability"].exists()


def test_unknown_celltype_col_raises():
    a = _adata()
    with pytest.raises(ValueError):
        estimate_reference_uncertainty(a, "nope")
