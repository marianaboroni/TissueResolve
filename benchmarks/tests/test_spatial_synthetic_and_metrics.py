"""Tests for spatial-fidelity metrics + structured spatial generator (PART 12/18).

Offline, deterministic.  Uses a tiny in-memory AnnData and small grids.
"""
import numpy as np
import pandas as pd
import pytest

from benchmarks.shared import spatial_metrics as SM

ad = pytest.importorskip("anndata")
from benchmarks.shared import synthetic_spatial as SSP  # noqa: E402


# --- spatial metrics ---------------------------------------------------------
def _grid_coords(n):
    r = np.repeat(np.arange(n), n)
    c = np.tile(np.arange(n), n)
    return np.c_[r, c].astype(float)


def test_morans_i_high_for_smooth_gradient_low_for_noise():
    coords = _grid_coords(10)
    smooth = coords[:, 0] / 9.0                       # gradient along rows
    rng = np.random.default_rng(0)
    noise = rng.random(coords.shape[0])
    nbr, w = SM.knn_weights(coords, k=4)
    i_smooth = SM.morans_i(smooth, nbr, w)
    i_noise = SM.morans_i(noise, nbr, w)
    assert i_smooth > 0.7           # strong positive autocorrelation
    assert i_noise < 0.3            # near-zero for noise
    assert i_smooth > i_noise


def test_gearys_c_complementary_to_morans():
    coords = _grid_coords(10)
    smooth = coords[:, 1] / 9.0
    nbr, w = SM.knn_weights(coords, k=4)
    c = SM.gearys_c(smooth, nbr, w)
    assert c < 0.5                  # low Geary's C ⇒ strong positive autocorr


def test_morans_preservation_and_oversmoothing_identity():
    coords = _grid_coords(8)
    df = pd.DataFrame({"a": coords[:, 0] / 7.0, "b": 1 - coords[:, 0] / 7.0},
                      index=[f"s{i}" for i in range(64)])
    pres = SM.morans_i_preservation(df, df, coords, k=4)
    assert pres["morans_i_mae"] == pytest.approx(0.0, abs=1e-9)
    assert SM.oversmoothing_score(df, df, coords, k=4) == pytest.approx(1.0, rel=1e-6)
    assert SM.local_rmse(df, df, coords, k=4) == pytest.approx(0.0, abs=1e-9)


def test_oversmoothing_score_gt_one_when_predicted_smoother():
    coords = _grid_coords(10)
    rng = np.random.default_rng(1)
    grad = coords[:, 0] / 9.0
    # truth: weak gradient + strong noise (modest positive autocorr)
    truth = pd.DataFrame({"a": grad + rng.random(100)}, index=[f"s{i}" for i in range(100)])
    # prediction: clean gradient (higher autocorr) ⇒ over-smoothed
    smooth_pred = pd.DataFrame({"a": grad}, index=truth.index)
    assert SM.oversmoothing_score(truth, smooth_pred, coords, k=4) > 1.0


def test_domain_recovery_ari_perfect_for_separable_domains():
    coords = _grid_coords(10)
    left = coords[:, 1] < 5
    pred = pd.DataFrame({"x": np.where(left, 0.9, 0.1), "y": np.where(left, 0.1, 0.9)},
                        index=[f"s{i}" for i in range(100)])
    labels = np.where(left, "L", "R")
    assert SM.domain_recovery_ari(pred, labels, seed=0) > 0.95


# --- structured generator ----------------------------------------------------
def _toy_adata(seed=0):
    rng = np.random.default_rng(seed)
    n, g = 900, 40
    donors = [f"D{i}" for i in range(6)]
    types = ["epi1", "epi2", "Tcell", "myel", "fibro", "endo"]
    obs = pd.DataFrame({"cell_type": rng.choice(types, n),
                        "donor_id": rng.choice(donors, n)})
    base = {t: rng.poisson(2 + 3 * i, g) + 1 for i, t in enumerate(types)}
    X = np.vstack([rng.poisson(base[t]) for t in obs["cell_type"]]).astype(float)
    a = ad.AnnData(X=X, obs=obs); a.var_names = [f"g{j}" for j in range(g)]
    return a, types


def test_spatial_targets_structure_and_domains():
    _, types = _toy_adata()
    mapping = {"epi1": "Epithelial", "epi2": "Epithelial", "Tcell": "T/NK",
               "myel": "Myeloid", "fibro": "Stromal/Fibroblast", "endo": "Endothelial"}
    targets, coords, domains, meta = SSP.build_spatial_targets(
        types, mapping, n_side=12, seed=0, rare_type="Tcell")
    assert targets.shape == (144, 6)
    assert np.allclose(targets.sum(axis=1), 1.0, atol=1e-6)
    assert set(np.unique(domains)) <= {"L", "R", "niche"}
    assert "niche" in domains                         # rare niche present
    # left half is Epithelial-richer than right half (sharp border)
    lefts = coords[:, 1] < 6
    epi = targets[["epi1", "epi2"]].sum(axis=1).to_numpy()
    assert epi[lefts].mean() > epi[~lefts].mean()


def test_generate_spatial_scenario_donor_disjoint_and_shapes():
    a, types = _toy_adata()
    mapping = {"epi1": "Epithelial", "epi2": "Epithelial", "Tcell": "T/NK",
               "myel": "Myeloid", "fibro": "Stromal/Fibroblast", "endo": "Endothelial"}
    query = ["D3", "D4", "D5"]
    sc = SSP.generate_spatial_scenario(a, types, query, mapping, n_side=10,
                                       cells_per_spot=20, seed=1, rare_type="Tcell")
    assert sc["Y"].shape == (100, 40)                 # spots × genes
    assert sc["truth"].shape == (100, 6)
    assert len(sc["array_row"]) == 100 and len(sc["lib_sizes"]) == 100
    assert np.allclose(sc["truth"].sum(axis=1), 1.0, atol=1e-6)
