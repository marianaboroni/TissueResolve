"""Synthetic validation of the calibrated identifiability certificate (Stage 2B)."""
from __future__ import annotations

import inspect

import numpy as np

from tissueresolve.results import ReferenceSignature
from tissueresolve.reference.identifiability_calibration import (
    calibrated_identifiability_certificate as ccert,
)


def _ref(blocks, G=200, hi=200.0, bg=2.0, cv=0.2):
    K = len(blocks)
    R = np.full((K, G), bg)
    for k, blk in enumerate(blocks):
        R[k, list(blk)] = hi
    R = R / R.sum(1, keepdims=True) * 1e6
    return ReferenceSignature(gene_names=[f"g{i}" for i in range(G)],
                              cell_types=[f"T{k}" for k in range(K)], R_cpm=R,
                              donor_cv=np.full((G, K), cv))


def test_orthogonal_all_resolvable():
    c = ccert(_ref([range(0, 40), range(40, 80), range(80, 120)]), library_size=1e7, n_sim=120, seed=0)
    assert (c.per_type.recoverability == "RESOLVABLE").all()
    assert len(c.clusters) == 0


def test_identical_pair_group_only():
    c = ccert(_ref([range(0, 40), range(0, 40), range(80, 120)]), library_size=1e7, n_sim=120, seed=0)
    per = c.per_type.set_index("cell_type")
    assert per.loc["T0", "recoverability"] in ("GROUP_ONLY", "UNRESOLVABLE")
    assert per.loc["T1", "recoverability"] in ("GROUP_ONLY", "UNRESOLVABLE")
    assert per.loc["T2", "recoverability"] == "RESOLVABLE"
    cl = c.clusters.iloc[0]
    assert set(cl["group"].split(";")) == {"T0", "T1"}
    assert cl["swap_rmse"] > 0.15                       # within-pair split not recoverable
    assert cl["group_sum_rmse"] < 0.05                  # but the pair SUM is recoverable


def test_nonnegativity_resolves_near_collinear():
    """A near-collinear pair with a small distinct signal is RESOLVABLE under nonnegativity —
    the refinement over the unconstrained linear certificate (focus 1)."""
    G = 200
    R = np.full((2, G), 2.0)
    R[0, 0:40] = 200.0
    R[1, 0:36] = 200.0; R[1, 40:44] = 200.0             # 36/40 shared, small distinct block
    R = R / R.sum(1, keepdims=True) * 1e6
    ref = ReferenceSignature(gene_names=[f"g{i}" for i in range(G)], cell_types=["T0", "T1"],
                             R_cpm=R, donor_cv=np.full((G, 2), 0.2))
    c = ccert(ref, library_size=1e7, n_sim=150, seed=0)
    assert (c.per_type.recoverability == "RESOLVABLE").all()
    assert c.clusters.iloc[0]["swap_rmse"] < 0.10


def test_three_identical_group():
    c = ccert(_ref([range(0, 40)] * 3 + [range(80, 120)]), library_size=1e7, n_sim=150, seed=0)
    cl = c.clusters
    assert (cl.n_members == 3).any()
    grp = set(cl[cl.n_members == 3].iloc[0]["group"].split(";"))
    assert grp == {"T0", "T1", "T2"}
    assert c.per_type.set_index("cell_type").loc["T3", "recoverability"] == "RESOLVABLE"


def test_donor_uncertainty_widens_pred_sd():
    lo = ccert(_ref([range(0, 40), range(40, 80), range(80, 120)], cv=0.05), library_size=1e7, n_sim=120, seed=0)
    hi = ccert(_ref([range(0, 40), range(40, 80), range(80, 120)], cv=0.8), library_size=1e7, n_sim=120, seed=0)
    assert hi.per_type.pred_sd.mean() > lo.per_type.pred_sd.mean()


def test_depth_improves_recovery():
    """For well-separated (orthogonal) types, deeper sequencing lowers counting-noise error."""
    blocks = [range(0, 40), range(40, 80), range(80, 120)]
    shallow = ccert(_ref(blocks, cv=0.0), library_size=5e3, n_sim=120, seed=1)
    deep = ccert(_ref(blocks, cv=0.0), library_size=1e8, n_sim=120, seed=1)
    assert deep.per_type.sim_rmse.mean() < shallow.per_type.sim_rmse.mean()


def test_deterministic():
    a = ccert(_ref([range(0, 40), range(0, 40)]), library_size=1e7, n_sim=100, seed=3)
    b = ccert(_ref([range(0, 40), range(0, 40)]), library_size=1e7, n_sim=100, seed=3)
    assert a.per_type.to_dict() == b.per_type.to_dict()


def test_uses_no_truth():
    params = inspect.signature(ccert).parameters
    assert not ({"truth", "true_proportions", "held_out", "y"} & set(params))
