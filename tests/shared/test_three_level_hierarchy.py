"""Offline tests for the three-level (broad → cell type → state) hierarchy."""
from __future__ import annotations

import warnings

import pandas as pd
import pytest

from tissueresolve.reference import three_level_hierarchy as TL

_S2C = {  # state -> cell type
    "CD8_eff": "CD8 T cell", "CD8_mem": "CD8 T cell",
    "CD4_naive": "CD4 T cell", "Treg": "CD4 T cell",
    "NK_dim": "NK cell",
    "macro_M1": "macrophage", "macro_M2": "macrophage",
}
_C2B = {  # cell type -> broad family
    "CD8 T cell": "T/NK", "CD4 T cell": "T/NK", "NK cell": "T/NK",
    "macrophage": "Myeloid",
}


def test_build_and_levels():
    h = TL.build_three_level_hierarchy(_S2C, _C2B)
    assert h.has_states
    assert h.states() == sorted(_S2C)
    assert set(h.cell_types()) == set(_C2B)
    assert h.broad_families() == ["Myeloid", "T/NK"]
    # state aggregates upward correctly
    assert h.state_to_broad()["Treg"] == "T/NK"
    assert h.states_of_celltype("CD8 T cell") == ["CD8_eff", "CD8_mem"]
    assert h.celltypes_of_broad("T/NK") == ["CD4 T cell", "CD8 T cell", "NK cell"]


def test_missing_celltype_mapping_raises():
    bad_c2b = dict(_C2B)
    del bad_c2b["NK cell"]            # a state's cell type now unmapped
    with pytest.raises(ValueError):
        TL.build_three_level_hierarchy(_S2C, bad_c2b)


def test_two_level_fallback_has_no_states():
    h = TL.build_three_level_from_two_level({"CD8 T cell": "T/NK",
                                             "macrophage": "Myeloid"})
    assert h.has_states is False
    # states identical to cell types (1:1)
    assert set(h.states()) == set(h.cell_types())
    assert h.state_to_celltype["CD8 T cell"] == "CD8 T cell"


def test_from_obs_three_level():
    obs = pd.DataFrame({
        "broad": ["T/NK", "T/NK", "Myeloid", "Myeloid"],
        "cell_type": ["CD8 T cell", "CD4 T cell", "macrophage", "macrophage"],
        "state": ["CD8_eff", "Treg", "macro_M1", "macro_M2"],
    })
    h = TL.build_three_level_from_obs(obs, "broad", "cell_type", "state")
    assert h.has_states
    assert h.state_to_celltype["macro_M2"] == "macrophage"
    assert h.celltype_to_broad["CD4 T cell"] == "T/NK"


def test_from_obs_inconsistent_broad_raises():
    obs = pd.DataFrame({"broad": ["T/NK", "Myeloid"],
                        "cell_type": ["X", "X"]})   # one cell type, two broads
    with pytest.raises(ValueError):
        TL.build_three_level_from_obs(obs, "broad", "cell_type")


def test_validation_flags_low_cells_and_single_donor():
    h = TL.build_three_level_hierarchy(_S2C, _C2B)
    counts = {s: 100 for s in _S2C}
    counts["NK_dim"] = 5            # low support
    donors = {ct: {"d1", "d2"} for ct in _C2B}
    donors["macrophage"] = {"d1"}  # single donor
    with warnings.catch_warnings():
        warnings.simplefilter("ignore")
        v = TL.validate_three_level_hierarchy(
            h, state_counts=counts, celltype_donors=donors,
            min_cells_per_state=20)
    assert bool(v.set_index("state").loc["NK_dim", "low_support"]) is True
    macro_rows = v[v["cell_type"] == "macrophage"]
    assert macro_rows["single_donor"].all()
    assert {"state", "cell_type", "broad_family", "low_support",
            "single_donor"} <= set(v.columns)


def test_validation_emits_warnings():
    h = TL.build_three_level_hierarchy(_S2C, _C2B)
    counts = {s: 5 for s in _S2C}   # all low
    with pytest.warns(UserWarning):
        TL.validate_three_level_hierarchy(h, state_counts=counts,
                                          min_cells_per_state=20)


def test_summary_and_outputs(tmp_path):
    h = TL.build_three_level_hierarchy(_S2C, _C2B)
    counts = {s: 50 for s in _S2C}
    summary = TL.summarize_three_level_hierarchy(h, state_counts=counts)
    tnk = summary.set_index("broad_family").loc["T/NK"]
    assert int(tnk["n_cell_types"]) == 3 and int(tnk["n_states"]) == 5
    with warnings.catch_warnings():
        warnings.simplefilter("ignore")
        v = TL.validate_three_level_hierarchy(h, state_counts=counts)
    paths = TL.write_three_level_outputs(h, v, summary, tmp_path)
    for key in ("hierarchy_three_level", "hierarchy_validation", "hierarchy_summary"):
        assert paths[key].exists()
    frame = pd.read_csv(paths["hierarchy_three_level"], sep="\t")
    assert {"state", "cell_type", "broad_family"} == set(frame.columns)
