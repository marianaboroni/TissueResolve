"""Synthetic validation of the identifiability certificate (Stage 2A)."""
from __future__ import annotations

import inspect

import numpy as np

from tissueresolve.reference.identifiability import identifiability_certificate as cert


def _R(blocks, G=80, hi=200.0, bg=2.0):
    """Build (K, G) CPM-normalised reference from per-type high-expression gene blocks."""
    rows = []
    for blk in blocks:
        v = np.full(G, bg); v[list(blk)] = hi
        rows.append(v / v.sum() * 1e6)
    return np.vstack(rows)


GENES = [f"g{i}" for i in range(80)]


def test_orthogonal_types_resolvable():
    R = _R([range(0, 20), range(20, 40), range(40, 60)])
    c = cert(R, ["A", "B", "C"], GENES, library_size=1e7)
    cls = dict(zip(c.per_type.cell_type, c.per_type.recoverability))
    assert all(cls[t] == "RESOLVABLE" for t in ("A", "B", "C")), cls
    assert c.condition_number < 5


def test_identical_profiles_unresolvable_and_merged():
    R = _R([range(0, 20), range(20, 40), range(0, 20)])   # A and X share the exact block
    c = cert(R, ["A", "B", "X"], GENES, library_size=1e7)
    per = c.per_type.set_index("cell_type")
    assert per.loc["A", "recoverability"] in ("GROUP_ONLY", "UNRESOLVABLE")
    assert per.loc["X", "recoverability"] in ("GROUP_ONLY", "UNRESOLVABLE")
    merge = set(per.loc["A", "recommended_merge"].split(";"))
    assert {"A", "X"} <= merge, per.loc["A"].to_dict()
    assert any({"A", "X"} <= set(g["group"].split(";")) for g in c.confounded_groups)


def test_near_collinear_depth_dependent():
    """A near-collinear pair is high-condition-number, and its resolvability depends on depth
    (the certificate is noise-aware, not a pure geometry test)."""
    R = _R([range(0, 20), range(20, 40)])
    R[1] = 0.99 * R[0] + 0.01 * R[1]                       # cos ~ 0.9999 collinear pair
    orth = cert(_R([range(0, 20), range(20, 40)]), ["A", "B"], GENES, library_size=1e5)
    shallow = cert(R, ["A", "B"], GENES, library_size=1e5)
    deep = cert(R, ["A", "B"], GENES, library_size=1e7)
    assert orth.condition_number < 5 < 50 < shallow.condition_number   # collinearity ⇒ ill-conditioned
    scls = dict(zip(shallow.per_type.cell_type, shallow.per_type.recoverability))
    assert scls["B"] in ("WEAKLY_RESOLVABLE", "GROUP_ONLY", "UNRESOLVABLE"), scls
    assert deep.noise_adjusted_rank >= shallow.noise_adjusted_rank      # depth recovers the small axis


def test_three_type_confounded_group():
    R = _R([range(0, 20)] * 3 + [range(40, 60)])          # A,B,C identical; D unique
    for i in range(3):                                     # tiny jitter so not exactly equal
        R[i] = R[i] * (1 + 0.001 * i)
    c = cert(R, ["A", "B", "C", "D"], GENES, library_size=1e7)
    grp = [set(g["group"].split(";")) for g in c.confounded_groups]
    assert any({"A", "B", "C"} <= g for g in grp), grp
    assert c.per_type.set_index("cell_type").loc["D", "recoverability"] == "RESOLVABLE"


def test_query_gene_loss_reduces_recoverability():
    """A and B share a block and differ ONLY in g70/g71; losing those two query genes makes them
    identical in query space, collapsing A from RESOLVABLE to a confounded group."""
    rows = []
    for extra in (70, 71):
        v = np.full(80, 2.0); v[list(range(0, 20))] = 200.0; v[extra] = 200.0
        rows.append(v / v.sum() * 1e6)
    v = np.full(80, 2.0); v[list(range(40, 60))] = 200.0; rows.append(v / v.sum() * 1e6)  # C unique
    R = np.vstack(rows)
    full = cert(R, ["A", "B", "C"], GENES, library_size=1e7)
    detect = [g for g in GENES if g not in {"g70", "g71"}]  # lose the only discriminating genes
    lost = cert(R, ["A", "B", "C"], GENES, query_detectable_genes=detect, library_size=1e7)
    assert lost.query_gene_retention["retention_frac"] < 1.0
    pf = full.per_type.set_index("cell_type"); pl = lost.per_type.set_index("cell_type")
    assert pf.loc["A", "recoverability"] == "RESOLVABLE"
    assert pl.loc["A", "recoverability"] in ("GROUP_ONLY", "UNRESOLVABLE")
    assert pl.loc["A", "rec_fraction"] < pf.loc["A", "rec_fraction"]


def test_depth_improves_recoverability():
    R = _R([range(0, 20), range(20, 40)])
    R[1] = 0.97 * R[0] + 0.03 * R[1]                       # collinear pair
    low = cert(R, ["A", "B"], GENES, library_size=1e5)
    high = cert(R, ["A", "B"], GENES, library_size=1e9)
    assert high.noise_adjusted_rank >= low.noise_adjusted_rank
    assert (high.per_type.theoretical_detection_floor.mean()
            <= low.per_type.theoretical_detection_floor.mean())


def test_donor_uncertainty_not_testable():
    R = _R([range(0, 20), range(20, 40), range(40, 60)])
    c = cert(R, ["A", "B", "C"], GENES, library_size=1e7,
             n_donors_by_type={"A": 1, "B": 5, "C": 5})
    per = c.per_type.set_index("cell_type")
    assert per.loc["A", "recoverability"] == "NOT_TESTABLE"
    assert per.loc["B", "recoverability"] == "RESOLVABLE"


def test_deterministic():
    R = _R([range(0, 20), range(20, 40), range(0, 20)])
    a = cert(R, ["A", "B", "X"], GENES); b = cert(R, ["A", "B", "X"], GENES)
    np.testing.assert_array_equal(a.singular_values, b.singular_values)
    assert a.per_type.to_dict() == b.per_type.to_dict()


def test_no_truth_used():
    """Certificate must depend only on R + query genes + depth — never on any truth/labels."""
    params = inspect.signature(cert).parameters
    assert "truth" not in params and "true_proportions" not in params and "y" not in params
