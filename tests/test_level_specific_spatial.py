"""Tests for level-specific spatial smoothing + boundary metrics (Task D).

Offline/deterministic: the solver-composition logic is tested via a tiny toy
reference where deconv_spatial runs fast; boundary/rare-niche metrics tested with
synthetic arrays.
"""
import numpy as np
import pandas as pd
import pytest

from benchmarks.shared import spatial_metrics as SM


def _grid(n):
    r = np.repeat(np.arange(n), n); c = np.tile(np.arange(n), n)
    return np.c_[r, c].astype(float)


# ---- boundary metrics ----
def test_boundary_f1_perfect_when_pred_matches_domains():
    coords = _grid(10)
    left = coords[:, 1] < 5
    idx = [f"s{i}" for i in range(100)]
    # truth & pred both: family A dominant on left, B on right
    truth = pd.DataFrame({"A": np.where(left, 0.9, 0.1), "B": np.where(left, 0.1, 0.9)}, index=idx)
    pred = truth.copy()
    dom = np.where(left, "L", "R")
    bm = SM.boundary_metrics(truth, pred, coords, dom)
    assert bm["boundary_f1"] > 0.8
    assert bm["edge_blurring"] == pytest.approx(0.0, abs=1e-9)  # pred==truth → no extra edge error


def test_edge_blurring_increases_when_prediction_smoothed_across_border():
    coords = _grid(10)
    left = coords[:, 1] < 5
    idx = [f"s{i}" for i in range(100)]
    truth = pd.DataFrame({"A": np.where(left, 1.0, 0.0), "B": np.where(left, 0.0, 1.0)}, index=idx)
    # over-smoothed prediction: blurred ramp across the border (column / 9)
    ramp = 1 - coords[:, 1] / 9.0
    pred = pd.DataFrame({"A": ramp, "B": 1 - ramp}, index=idx)
    dom = np.where(left, "L", "R")
    bm_sharp = SM.boundary_metrics(truth, truth.copy(), coords, dom)
    bm_blur = SM.boundary_metrics(truth, pred, coords, dom)
    assert bm_blur["edge_blurring"] > bm_sharp["edge_blurring"]


def test_rare_niche_sensitivity():
    idx = [f"s{i}" for i in range(20)]
    dom = np.array(["niche"] * 5 + ["L"] * 15)
    pred = pd.DataFrame({"rare": [0.2] * 4 + [0.0] + [0.0] * 15, "other": [0.8] * 20}, index=idx)
    truth = pred.copy()
    s = SM.rare_niche_sensitivity(truth, pred, dom, "rare", thresh=0.01)
    assert s == pytest.approx(4 / 5)   # 4 of 5 niche spots predict rare > 0.01


# ---- level-specific module (logic, tiny toy reference) ----
def _toy_ref_and_spots():
    ad = pytest.importorskip("anndata")
    import tissueresolve as tr
    from tissueresolve.reference.build import ReferenceBuilder
    from tissueresolve.config import ReferenceConfig
    rng = np.random.default_rng(0)
    G = 50
    types = ["A1", "A2", "B1"]
    # build a tiny AnnData reference
    base = {t: rng.poisson(3 + 5 * i, G) + 1 for i, t in enumerate(types)}
    cells = np.repeat(types, 40)
    X = np.vstack([rng.poisson(base[t]) for t in cells]).astype(float)
    a = ad.AnnData(X=X, obs=pd.DataFrame({"cell_type": cells}))
    a.var_names = [f"g{j}" for j in range(G)]
    ref = ReferenceBuilder(ReferenceConfig()).build_from_adata(a, estimate_overdispersion=True)
    # synthetic spots (small grid)
    n = 16
    R = ref.as_R_cpm()
    true = rng.dirichlet(np.ones(R.shape[0]), size=n)
    lib = rng.integers(800, 1500, n).astype("float32")
    lam = true @ R; lam = lam / lam.sum(1, keepdims=True)
    Y = rng.poisson(np.clip(lam * lib[:, None], 0, None)).astype("float32")
    rows = np.repeat(np.arange(4), 4); cols = np.tile(np.arange(4), 4)
    return ref, Y, rows, cols, lib, [f"g{j}" for j in range(G)], [f"sp{i}" for i in range(n)]


def test_level_specific_mass_conservation_and_separate_lambdas():
    from tissueresolve.experimental.spatial_smoothing import fit_level_specific_spatial
    ref, Y, rows, cols, lib, genes, ids = _toy_ref_and_spots()
    mapping = {"A1": "A", "A2": "A", "B1": "B"}
    import warnings
    with warnings.catch_warnings():
        warnings.simplefilter("ignore")
        r = fit_level_specific_spatial(Y, ref, rows, cols, lib, genes, mapping,
                                       lambda_broad=0.1, lambda_fine=0.0, spot_ids=ids)
    assert r.lambda_broad == 0.1 and r.lambda_fine == 0.0
    # combined fine sums to 1 per spot (mass-consistent with broad)
    assert np.allclose(r.fine_proportions.sum(axis=1), 1.0, atol=1e-6)
    # within-family conditional sums to 1 per family per spot
    a_share = r.conditional_proportions[["A1", "A2"]].sum(axis=1)
    assert np.allclose(a_share, 1.0, atol=1e-6)
    assert (r.fine_proportions.to_numpy() >= -1e-9).all()


def test_both_zero_is_no_smoothing_runs():
    from tissueresolve.experimental.spatial_smoothing import fit_level_specific_spatial
    ref, Y, rows, cols, lib, genes, ids = _toy_ref_and_spots()
    mapping = {"A1": "A", "A2": "A", "B1": "B"}
    import warnings
    with warnings.catch_warnings():
        warnings.simplefilter("ignore")
        r = fit_level_specific_spatial(Y, ref, rows, cols, lib, genes, mapping,
                                       lambda_broad=0.0, lambda_fine=0.0, spot_ids=ids)
    assert r.lambda_broad == 0.0 and r.lambda_fine == 0.0
    assert r.fine_proportions.shape == (16, 3)
