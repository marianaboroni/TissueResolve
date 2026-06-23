"""Tests for the held-out-donor pseudobulk generator (PART 2A.1 / PART 18).

Uses a tiny in-memory AnnData — fast, offline, deterministic.  Verifies the
non-negotiable property (mixture donors disjoint from reference donors),
count-level integrality, ground-truth normalisation, scenario shaping, and
reproducibility.
"""
import numpy as np
import pandas as pd
import pytest

ad = pytest.importorskip("anndata")

from benchmarks.shared import synthetic_holdout as SH


def _toy_adata(seed=0):
    rng = np.random.default_rng(seed)
    n_cells, n_genes = 720, 40
    donors = [f"D{i}" for i in range(6)]
    types = ["A", "B", "C", "Bsim"]  # Bsim ~ similar to B
    obs = pd.DataFrame({
        "cell_type": rng.choice(types, n_cells),
        "donor_id": rng.choice(donors, n_cells),
        "assay": rng.choice(["v2", "v3"], n_cells),
    })
    # type-specific mean expression so counts carry signal
    base = {t: rng.poisson(2 + 4 * i, n_genes) + 1 for i, t in enumerate(types)}
    X = np.vstack([rng.poisson(base[t]) for t in obs["cell_type"]]).astype(float)
    a = ad.AnnData(X=X, obs=obs)
    a.var_names = [f"g{j}" for j in range(n_genes)]
    return a


def test_donor_split_is_disjoint():
    a = _toy_adata()
    ref, qry = SH.split_donors(a, "donor_id", ref_frac=0.5, seed=1)
    assert set(ref).isdisjoint(set(qry))
    assert set(ref) | set(qry) == set(a.obs["donor_id"].astype(str))
    assert ref and qry


def test_generate_balanced_truth_and_disjoint_donors():
    a = _toy_adata()
    ds = SH.generate_scenario(a, celltype_col="cell_type", donor_col="donor_id",
                              scenario="balanced", n_samples=8, seed=2,
                              cells_per_sample=300)
    # mixture donors disjoint from reference donors — the core requirement
    assert set(ds.query_donors).isdisjoint(set(ds.reference_donors))
    # counts are integer, non-negative, genes × samples
    assert ds.counts.shape == (a.n_vars, 8)
    assert (ds.counts.to_numpy() >= 0).all()
    assert np.allclose(ds.counts.to_numpy(), np.round(ds.counts.to_numpy()))
    # both ground truths normalised per sample
    assert np.allclose(ds.true_cell_fractions.sum(axis=1), 1.0, atol=1e-6)
    assert np.allclose(ds.true_mrna_proportions.sum(axis=1), 1.0, atol=1e-6)
    # cell fraction and mRNA proportion differ (mRNA-content effect)
    assert not np.allclose(ds.true_cell_fractions.to_numpy(),
                           ds.true_mrna_proportions.to_numpy())


@pytest.mark.parametrize("level", [0.001, 0.005, 0.01, 0.02, 0.05])
def test_rare_scenario_sets_target_level(level):
    a = _toy_adata()
    tgt = SH.build_target_proportions(["A", "B", "C", "Bsim"], 5, "rare",
                                      seed=3, rare_type="C", rare_level=level)
    assert np.allclose(tgt["C"].to_numpy(), level)
    assert np.allclose(tgt.sum(axis=1), 1.0)


def test_similar_subtypes_concentrate_on_pair():
    tgt = SH.build_target_proportions(["A", "B", "C", "Bsim"], 6, "similar_subtypes",
                                      seed=4, similar_pair=("B", "Bsim"))
    pair_mass = tgt["B"] + tgt["Bsim"]
    assert (pair_mass >= 0.5).all()


def test_reproducible_given_seed():
    a = _toy_adata()
    kw = dict(celltype_col="cell_type", donor_col="donor_id",
              scenario="imbalanced", n_samples=6, seed=7, cells_per_sample=200)
    d1 = SH.generate_scenario(a, **kw)
    d2 = SH.generate_scenario(a, **kw)
    pd.testing.assert_frame_equal(d1.counts, d2.counts)
    pd.testing.assert_frame_equal(d1.true_mrna_proportions, d2.true_mrna_proportions)


def test_cross_platform_split_by_assay():
    a = _toy_adata()
    ref, qry = SH.split_donors(a, "donor_id", seed=0, assay_col="assay",
                               cross_platform=True)
    assert set(ref).isdisjoint(set(qry))
