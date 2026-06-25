"""Offline tests for the state-aware three-level bulk solver."""
from __future__ import annotations

import warnings

import numpy as np
import pandas as pd
import pytest

from tissueresolve.results import ReferenceSignature
from tissueresolve.reference.three_level_hierarchy import build_three_level_hierarchy
from tissueresolve.bulk.state_aware_hierarchical import (
    run_state_aware_hierarchical_bulk)

# state -> cell type -> broad
_S2C = {"A1": "A", "A2": "A", "B1": "B", "B2": "B", "C1": "C",
        "D1": "D", "E1": "E"}
_C2B = {"A": "X", "B": "X", "C": "Y", "D": "Z", "E": "Z"}


@pytest.fixture
def state_ref_and_bulk():
    rng = np.random.default_rng(0)
    G = 90
    genes = [f"g{i}" for i in range(G)]

    def prof(block, lvl=200.0):
        v = np.full(G, 2.0)
        v[list(block)] = lvl
        return v

    A1 = prof(range(0, 15))
    A2 = A1 * (1.0 + rng.normal(0, 0.003, G))     # A1≈A2 → state non-separable
    B1 = prof(range(20, 35))
    B2 = prof(range(40, 55))                       # B1,B2 separable; A vs B separable
    C1 = prof(range(60, 75))
    D1 = prof(range(76, 86))
    E1 = D1 * (1.0 + rng.normal(0, 0.003, G))      # D≈E → cell types non-separable in Z
    states = ["A1", "A2", "B1", "B2", "C1", "D1", "E1"]
    R = np.vstack([A1, A2, B1, B2, C1, D1, E1]).astype(np.float32)
    ref = ReferenceSignature(
        gene_names=genes, cell_types=states, R_cpm=R,
        R_log=np.log1p(R).astype(np.float32),
        phi_g=np.full(G, 5.0, dtype=np.float32),
        n_cells_per_type={s: 100 for s in states})

    # pseudobulk: 3 samples mixed from known state proportions (genes × samples)
    p = np.array([[0.2, 0.1, 0.15, 0.1, 0.15, 0.15, 0.15],
                  [0.1, 0.2, 0.2, 0.1, 0.1, 0.2, 0.1],
                  [0.25, 0.05, 0.1, 0.2, 0.2, 0.1, 0.1]])
    bulk = pd.DataFrame((p @ R).T, index=genes, columns=["s0", "s1", "s2"])
    return ref, bulk


@pytest.fixture
def hierarchy():
    return build_three_level_hierarchy(_S2C, _C2B)


def _run(ref, bulk, hierarchy, **kw):
    with warnings.catch_warnings():
        warnings.simplefilter("ignore")
        return run_state_aware_hierarchical_bulk(
            bulk, ref, hierarchy, min_discriminating_genes=5, **kw)


def test_three_level_proportions_returned(state_ref_and_bulk, hierarchy):
    ref, bulk = state_ref_and_bulk
    res = _run(ref, bulk, hierarchy)
    assert set(res.broad_proportions.columns) >= {"X", "Y", "Z"}
    assert res.state_proportions is not None
    # broad/cell-type/state frames all indexed by the 3 samples
    assert list(res.state_proportions.index) == ["s0", "s1", "s2"]


def test_mass_preserved_at_each_level(state_ref_and_bulk, hierarchy):
    ref, bulk = state_ref_and_bulk
    res = _run(ref, bulk, hierarchy)
    assert np.allclose(res.cell_type_proportions.sum(axis=1).to_numpy(), 1.0, atol=1e-6)
    assert np.allclose(res.state_proportions.sum(axis=1).to_numpy(), 1.0, atol=1e-6)


def test_unresolved_assigned_at_correct_level(state_ref_and_bulk, hierarchy):
    ref, bulk = state_ref_and_bulk
    res = _run(ref, bulk, hierarchy)
    levels = res.unresolved_by_level
    broad_unres = set(levels.loc[levels["level"] == "broad_family", "label"])
    ct_unres = set(levels.loc[levels["level"] == "cell_type", "label"])
    # Z's cell types (D,E) are non-separable -> unresolved at the BROAD level
    assert "unresolved_Z" in broad_unres
    # A's states (A1,A2) are non-separable -> unresolved at the CELL-TYPE level
    assert "unresolved_A" in ct_unres
    # and the unresolved_Z mass is NOT mislabelled as a cell-type-level unresolved
    assert "unresolved_Z" not in ct_unres


def test_allocation_metadata_records_panel_source(state_ref_and_bulk, hierarchy):
    ref, bulk = state_ref_and_bulk
    # pass explicit panels and check they're recorded
    panels = {"X": ref.gene_names[:30], "B": ref.gene_names[20:55]}
    res = _run(ref, bulk, hierarchy,
               celltype_panels={"X": ref.gene_names[:55], "Z": ref.gene_names[76:86]},
               state_panels={"A": ref.gene_names[:15], "B": ref.gene_names[20:55]})
    md = res.allocation_metadata
    assert {"level", "group", "resolved", "n_panel_genes", "panel_source"} <= set(md.columns)
    assert (md["panel_source"] == "cell_type_within_family").any()
    assert (md["panel_source"] == "state_within_cell_type").any()
    assert res.metadata["used_celltype_panels"] and res.metadata["used_state_panels"]


def test_two_level_fallback_no_state(state_ref_and_bulk):
    """has_states=False → only broad→cell_type, no state frame, mass preserved."""
    from tissueresolve.reference.three_level_hierarchy import build_three_level_from_two_level
    ref, bulk = state_ref_and_bulk
    # treat the 7 'states' as cell types directly (no state level)
    h = build_three_level_from_two_level(
        {s: ("X" if s in ("A1", "A2", "B1", "B2") else
             ("Y" if s == "C1" else "Z")) for s in ref.cell_types})
    res = _run(ref, bulk, h)
    assert res.state_proportions is None
    assert np.allclose(res.cell_type_proportions.sum(axis=1).to_numpy(), 1.0, atol=1e-6)


def test_outputs_written(tmp_path, state_ref_and_bulk, hierarchy):
    from tissueresolve.bulk.state_aware_hierarchical import write_state_aware_outputs
    ref, bulk = state_ref_and_bulk
    res = _run(ref, bulk, hierarchy)
    paths = write_state_aware_outputs(res, tmp_path)
    for key in ("broad_proportions", "cell_type_proportions", "state_proportions",
                "unresolved_mass_by_level", "hierarchical_allocation_metadata"):
        assert paths[key].exists()
