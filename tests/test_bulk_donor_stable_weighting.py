"""Tests for experimental donor-stable bulk gene weighting (pure scoring)."""
from __future__ import annotations

import sys
from pathlib import Path

import numpy as np
import pandas as pd
import pytest

REPO = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(REPO / "src"))

from tissueresolve.experimental.bulk_donor_stable_weighting import (
    compute_donor_stable_gene_weights,
)


def _toy_adata(seed=0):
    import anndata as ad
    rng = np.random.default_rng(seed)
    # 4 cell types? keep 2 types, 3 donors. Genes:
    #  STABLE_MARKER: type-specific + consistent across donors
    #  UNSTABLE: type-specific mean but huge donor-to-donor variation
    #  HOUSEKEEPING: uniform across types
    types, donors, rows = [], [], []
    for t in ["T1", "T2"]:
        for d in ["d1", "d2", "d3"]:
            donor_shift = {"d1": 1.0, "d2": 1.0, "d3": 1.0}[d]
            unstable_shift = {"d1": 0.1, "d2": 5.0, "d3": 20.0}[d]  # donor-driven
            for _ in range(20):
                marker = (10.0 if t == "T1" else 0.5) * donor_shift + rng.normal(0, 0.2)
                unstable = (8.0 if t == "T1" else 0.5) * unstable_shift + rng.normal(0, 0.2)
                house = 5.0 + rng.normal(0, 0.2)
                rows.append([max(marker, 0), max(unstable, 0), max(house, 0)])
                types.append(t); donors.append(d)
    X = np.array(rows)
    a = ad.AnnData(X=X)
    a.var_names = ["STABLE_MARKER", "UNSTABLE", "HOUSEKEEPING"]
    a.obs["cell_type"] = types
    a.obs["donor_id"] = donors
    return a


def test_weights_finite_nonnegative_no_removal():
    a = _toy_adata()
    genes = ["STABLE_MARKER", "UNSTABLE", "HOUSEKEEPING", "NOT_IN_REF"]
    w = compute_donor_stable_gene_weights(a, "cell_type", "donor_id", genes)
    assert list(w.index) == genes           # no removal; order preserved
    assert np.all(np.isfinite(w.to_numpy())) and (w.to_numpy() >= 0).all()


def test_unstable_donor_gene_downweighted_vs_stable_marker():
    a = _toy_adata()
    w = compute_donor_stable_gene_weights(
        a, "cell_type", "donor_id", ["STABLE_MARKER", "UNSTABLE", "HOUSEKEEPING"])
    assert w["STABLE_MARKER"] > w["UNSTABLE"]      # donor-unstable downweighted
    assert w["STABLE_MARKER"] > w["HOUSEKEEPING"]  # non-specific downweighted


def test_missing_donor_column_fallback():
    a = _toy_adata()
    w = compute_donor_stable_gene_weights(a, "cell_type", None,
                                          ["STABLE_MARKER", "UNSTABLE"])
    assert w.attrs["donor_fallback"] is True
    assert np.all(np.isfinite(w.to_numpy()))


def test_missing_gene_gets_neutral_weight_not_removed():
    a = _toy_adata()
    w = compute_donor_stable_gene_weights(a, "cell_type", "donor_id",
                                          ["STABLE_MARKER", "NOT_IN_REF"])
    assert "NOT_IN_REF" in w.index and w["NOT_IN_REF"] == 1.0
    assert w.attrs["n_missing"] == 1


def test_protocol_blacklist_downweights():
    a = _toy_adata()
    w = compute_donor_stable_gene_weights(
        a, "cell_type", "donor_id", ["STABLE_MARKER", "HOUSEKEEPING"],
        protocol_blacklist={"STABLE_MARKER"})
    w2 = compute_donor_stable_gene_weights(
        a, "cell_type", "donor_id", ["STABLE_MARKER", "HOUSEKEEPING"])
    # weights are max-normalized, so compare the relative weight: blacklisting
    # STABLE_MARKER must lower it relative to HOUSEKEEPING.
    ratio_bl = w["STABLE_MARKER"] / w["HOUSEKEEPING"]
    ratio_no = w2["STABLE_MARKER"] / w2["HOUSEKEEPING"]
    assert ratio_bl < ratio_no
