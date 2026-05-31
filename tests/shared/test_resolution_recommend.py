"""
Tests for the resolution-aware recommendation system (Stage 8).

Synthetic, offline.  Covers the non-separable graph, connected-component merge
families, family-name heuristics, the "no merge when separable" case, post-hoc
aggregation (mass preservation, fine estimates untouched), and the improved
separability warning text.
"""
from __future__ import annotations

import warnings

import numpy as np
import pandas as pd
import pytest

from tissueresolve.reference.hierarchy import (
    aggregate_predictions_by_family,
    compare_fine_vs_merged_predictions,
)
from tissueresolve.reference.resolution import (
    assign_resolution_families,
    build_nonseparable_graph,
    family_label_for_members,
    recommend_cell_type_merges,
)
from tissueresolve.reference.separability import SeparabilityWarning, compute_separability
from tissueresolve.results import ReferenceSignature


def _ref(profiles: dict[str, np.ndarray]) -> ReferenceSignature:
    cts = list(profiles)
    G = len(next(iter(profiles.values())))
    R = np.vstack([profiles[c] for c in cts]).astype(np.float32)
    return ReferenceSignature(
        gene_names=[f"G{i:03d}" for i in range(G)], cell_types=cts,
        R_cpm=R, R_log=np.log1p(R).astype(np.float32),
        phi_g=np.full(G, 5.0, dtype=np.float32),
        n_cells_per_type={c: 100 for c in cts})


def _grouped_confusable_ref():
    """3 confusable families with disjoint blocks across families."""
    G = 30
    base = np.full(G, 5.0)

    def block(lo, hi, val=1000.0):
        v = base.copy()
        v[lo:hi] = val
        return v

    return _ref({
        # T/NK family — block 0:10, nearly identical
        "CD4-positive helper T cell": block(0, 10, 1000),
        "CD8-positive, alpha-beta T cell": block(0, 10, 1005),
        "natural killer cell": block(0, 10, 995),
        # myeloid family — block 10:20
        "macrophage": block(10, 20, 1000),
        "monocyte": block(10, 20, 1004),
        "dendritic cell": block(10, 20, 996),
        # endothelial family — block 20:30
        "capillary endothelial cell": block(20, 30, 1000),
        "vein endothelial cell": block(20, 30, 1003),
    })


def _separable_ref():
    G = 30
    base = np.full(G, 5.0)
    profiles = {}
    for k, name in enumerate(["A", "B", "C"]):
        v = base.copy()
        v[k * 10:(k + 1) * 10] = 1000.0
        profiles[name] = v
    return _ref(profiles)


# ---------------------------------------------------------------------------
# Family-name heuristics
# ---------------------------------------------------------------------------


@pytest.mark.parametrize("members,expected", [
    (["CD4-positive helper T cell", "CD8-positive, alpha-beta T cell",
      "regulatory T cell", "natural killer cell", "lymphocyte"], "T/NK lymphocytes"),
    (["macrophage", "monocyte", "dendritic cell"], "myeloid cells"),
    (["capillary endothelial cell", "vein endothelial cell",
      "endothelial cell of artery"], "endothelial cells"),
    (["luminal epithelial cell of mammary gland", "basal cell",
      "mammary gland epithelial cell"], "epithelial cells"),
    (["pericyte", "vascular associated smooth muscle cell"], "mural cells"),
])
def test_family_label_heuristics(members, expected):
    assert family_label_for_members(members) == expected


def test_family_label_singleton_is_self():
    assert family_label_for_members(["macrophage"]) == "macrophage"


# ---------------------------------------------------------------------------
# Graph + components + recommendations
# ---------------------------------------------------------------------------


def test_nonseparable_graph_links_within_families():
    ref = _grouped_confusable_ref()
    with warnings.catch_warnings():
        warnings.simplefilter("ignore")
        sep = compute_separability(ref)
    g = build_nonseparable_graph(sep, list(ref.cell_types))
    # CD4 and CD8 are linked; CD4 and macrophage are not.
    assert "CD8-positive, alpha-beta T cell" in g["CD4-positive helper T cell"]
    assert "macrophage" not in g["CD4-positive helper T cell"]


def test_recommend_merges_components_become_families():
    ref = _grouped_confusable_ref()
    with warnings.catch_warnings():
        warnings.simplefilter("ignore")
        sep = compute_separability(ref)
    merges = recommend_cell_type_merges(sep, list(ref.cell_types))
    names = set(merges["family_name"])
    assert {"T/NK lymphocytes", "myeloid cells", "endothelial cells"} <= names
    assert (merges["n_members"] >= 2).all()


def test_assign_families_mapping_covers_all_types():
    ref = _grouped_confusable_ref()
    with warnings.catch_warnings():
        warnings.simplefilter("ignore")
        sep = compute_separability(ref)
    mapping = assign_resolution_families(sep, list(ref.cell_types))
    assert set(mapping) == set(ref.cell_types)
    assert mapping["CD4-positive helper T cell"] == "T/NK lymphocytes"
    assert mapping["macrophage"] == "myeloid cells"


def test_no_merge_when_separable():
    ref = _separable_ref()
    with warnings.catch_warnings():
        warnings.simplefilter("ignore")
        sep = compute_separability(ref)
    merges = recommend_cell_type_merges(sep, list(ref.cell_types))
    assert len(merges) == 0
    mapping = assign_resolution_families(sep, list(ref.cell_types))
    assert mapping == {c: c for c in ref.cell_types}  # everyone maps to self


# ---------------------------------------------------------------------------
# Post-hoc aggregation
# ---------------------------------------------------------------------------


def test_aggregate_preserves_total_mass():
    rng = np.random.default_rng(0)
    props = pd.DataFrame(rng.dirichlet(np.ones(4), size=5),
                         index=[f"s{i}" for i in range(5)],
                         columns=["a", "b", "c", "d"])
    mapping = {"a": "X", "b": "X", "c": "Y", "d": "Y"}
    fam = aggregate_predictions_by_family(props, mapping)
    assert list(fam.columns) == ["X", "Y"]
    np.testing.assert_allclose(fam.sum(axis=1).to_numpy(),
                               props.sum(axis=1).to_numpy(), atol=1e-9)


def test_compare_does_not_overwrite_fine():
    props = pd.DataFrame([[0.6, 0.4]], index=["s0"], columns=["a", "b"])
    out = compare_fine_vs_merged_predictions(props, {"a": "X", "b": "X"})
    pd.testing.assert_frame_equal(out["fine"], props)  # unchanged
    assert out["merge_stage"] == "post_hoc_aggregation"
    assert out["family"].loc["s0", "X"] == pytest.approx(1.0)


# ---------------------------------------------------------------------------
# Improved warning text
# ---------------------------------------------------------------------------


def test_separability_warning_is_actionable():
    ref = _grouped_confusable_ref()
    with pytest.warns(SeparabilityWarning) as rec:
        compute_separability(ref)
    msg = str(rec[0].message)
    assert "Recommended action" in msg
    assert "recommended_merges.tsv" in msg
    assert "--resolution-mode" in msg
    assert "merge manually" not in msg  # no bare "merge manually" instruction
