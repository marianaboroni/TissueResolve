"""
Controlled synthetic three-level (broad → cell type → state) benchmark.

The real breast-cancer reference has only broad + fine labels, so true
state-level inference cannot be validated there.  This module simulates a
reference with a KNOWN broad→cell_type→state hierarchy and bulk mixtures with
KNOWN proportions at all three levels, under scenarios that vary state
separability — so we can ask, honestly and offline, whether state-aware
deconvolution (a) improves true state-level accuracy where states are separable,
(b) preserves broad/cell-type accuracy, and (c) keeps non-separable states as
unresolved mass instead of inventing precision.

Everything is synthetic; no real data, no network.
"""
from __future__ import annotations

import time
from dataclasses import dataclass, field
from typing import Optional

import numpy as np
import pandas as pd

# fixed hierarchy: broad X {A(A1,A2), B(B1,B2)}, broad Y {C(C1), D(D1,D2)}
_STATE_TO_CT = {"A1": "A", "A2": "A", "B1": "B", "B2": "B",
                "C1": "C", "D1": "D", "D2": "D"}
_CT_TO_BROAD = {"A": "X", "B": "X", "C": "Y", "D": "Y"}
_STATES = list(_STATE_TO_CT)

SCENARIOS = ["separable", "weak", "non_separable", "missing_state", "imbalanced"]
# which states are genuinely sub-resolvable (truth) per scenario — used to score
# false-resolution and unresolved precision/recall.
_RESOLVABLE_STATES = {
    "separable": set(_STATES),
    "weak": set(_STATES),
    "non_separable": set(),                 # no state is truly separable
    "missing_state": set(_STATES),
    "imbalanced": set(_STATES),
}


@dataclass
class SyntheticScenario:
    name: str
    adata: object                # per-cell AnnData (broad/cell_type/state obs)
    state_cell_types: list       # states present IN the reference
    truth_state: pd.DataFrame    # samples × states (known; may include out-of-ref)
    truth_celltype: pd.DataFrame
    truth_broad: pd.DataFrame
    bulk: pd.DataFrame           # genes × samples
    missing_states: list = field(default_factory=list)
    resolvable_states: set = field(default_factory=set)


def _state_level(name: str) -> float:
    return {"separable": 160.0, "weak": 10.0, "non_separable": 0.0,
            "missing_state": 160.0, "imbalanced": 160.0}.get(name, 160.0)


def simulate_scenario(name: str, *, seed: int = 0, n_cells_per_state: int = 80,
                      n_samples: int = 12, n_genes: int = 150):
    """Simulate one scenario; returns a :class:`SyntheticScenario`."""
    if name not in SCENARIOS:
        raise ValueError(f"unknown scenario {name!r}; choose from {SCENARIOS}")
    import anndata as ad
    rng = np.random.default_rng(seed)
    genes = [f"g{i}" for i in range(n_genes)]
    gi = {g: i for i, g in enumerate(genes)}
    st_lvl = _state_level(name)

    # gene blocks
    broad_block = {"X": range(0, 15), "Y": range(15, 30)}
    ct_block = {"A": range(30, 40), "B": range(40, 50),
                "C": range(50, 60), "D": range(60, 70)}
    state_block = {"A1": range(70, 75), "A2": range(75, 80),
                   "B1": range(80, 85), "B2": range(85, 90),
                   "D1": range(90, 95), "D2": range(95, 100)}

    def cell(state: str) -> np.ndarray:
        ct = _STATE_TO_CT[state]
        broad = _CT_TO_BROAD[ct]
        v = rng.poisson(0.3, n_genes).astype(np.float64)
        v[list(broad_block[broad])] += rng.poisson(15, 15)
        v[list(ct_block[ct])] += rng.poisson(12, 10)
        if state in state_block and st_lvl > 0:
            blk = list(state_block[state])
            v[blk] += rng.poisson(st_lvl, len(blk))
        return v

    # per-cell AnnData (all 7 states); for missing_state we drop D2 from the
    # REFERENCE later but keep it in the bulk mixture.
    rows, broad_l, ct_l, st_l = [], [], [], []
    for s in _STATES:
        for _ in range(n_cells_per_state):
            rows.append(cell(s)); st_l.append(s)
            ct_l.append(_STATE_TO_CT[s]); broad_l.append(_CT_TO_BROAD[_STATE_TO_CT[s]])
    X = np.vstack(rows)
    obs = pd.DataFrame({"broad": broad_l, "cell_type": ct_l, "state": st_l},
                       index=[f"c{i}" for i in range(len(st_l))])
    adata = ad.AnnData(X=X, obs=obs, var=pd.DataFrame(index=genes))

    # reference states (drop one for missing_state scenario)
    missing = ["D2"] if name == "missing_state" else []
    ref_states = [s for s in _STATES if s not in missing]
    if missing:
        adata = adata[~adata.obs["state"].isin(missing)].copy()

    # mean state profiles (for bulk simulation) — use ALL states incl. missing
    prof = {}
    for s in _STATES:
        prof[s] = X[np.array(st_l) == s].mean(axis=0)

    # known sample proportions over ALL states (incl. missing for that scenario)
    if name == "imbalanced":
        alpha = np.array([8, 0.3, 4, 0.3, 2, 4, 0.3])      # skewed
    else:
        alpha = np.ones(len(_STATES))
    p = rng.dirichlet(alpha, size=n_samples)
    truth_state = pd.DataFrame(p, columns=_STATES,
                               index=[f"s{i}" for i in range(n_samples)])
    bulk_mat = truth_state.to_numpy() @ np.vstack([prof[s] for s in _STATES])
    bulk = pd.DataFrame(bulk_mat.T, index=genes, columns=truth_state.index)

    truth_ct = truth_state.T.groupby(_STATE_TO_CT).sum().T
    truth_broad = truth_ct.T.groupby(_CT_TO_BROAD).sum().T

    return SyntheticScenario(
        name=name, adata=adata, state_cell_types=ref_states,
        truth_state=truth_state, truth_celltype=truth_ct, truth_broad=truth_broad,
        bulk=bulk, missing_states=missing,
        resolvable_states=_RESOLVABLE_STATES[name])


# --------------------------------------------------------------------------- #
# reference + hierarchy from the simulated AnnData
# --------------------------------------------------------------------------- #
def build_state_reference(scenario: SyntheticScenario):
    """A state-level ``ReferenceSignature`` (cell_types = states) from the cells."""
    from tissueresolve.results import ReferenceSignature
    ad = scenario.adata
    X = np.asarray(ad.X, dtype=np.float64)
    states = scenario.state_cell_types
    st = ad.obs["state"].to_numpy()
    R = np.vstack([X[st == s].mean(axis=0) for s in states]).astype(np.float32)
    return ReferenceSignature(
        gene_names=list(ad.var_names), cell_types=list(states), R_cpm=R,
        R_log=np.log1p(R).astype(np.float32),
        phi_g=np.full(R.shape[1], 5.0, dtype=np.float32),
        n_cells_per_type={s: int((st == s).sum()) for s in states})


def three_level_hierarchy(scenario: SyntheticScenario):
    from tissueresolve.reference.three_level_hierarchy import build_three_level_hierarchy
    s2c = {s: _STATE_TO_CT[s] for s in scenario.state_cell_types}
    c2b = {c: _CT_TO_BROAD[c] for c in set(s2c.values())}
    return build_three_level_hierarchy(s2c, c2b)


def state_to_broad_mapping(scenario: SyntheticScenario) -> dict:
    return {s: _CT_TO_BROAD[_STATE_TO_CT[s]] for s in scenario.state_cell_types}
