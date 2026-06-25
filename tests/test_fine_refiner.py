"""Unit + invariant tests for the experimental FineGranularityRefiner.

Offline, deterministic, no generated data.  Covers the mandatory constraints:
mass conservation (within-family + total), simplex outputs, residual centering,
finite non-negative contrast weights, missing-donor handling, support-based
fallback, calibration train/test separation + simplex projection, default-off
behaviour, and refinement-status metadata.  Also verifies the canonical order
(refine BEFORE soft gating, soft gating applied exactly once).
"""
import numpy as np
import pandas as pd
import pytest

from tissueresolve.results import ReferenceSignature
from tissueresolve.experimental.soft_hierarchy.fine_refiner import (
    FineGranularityRefiner, FineRefinerConfig, SpilloverCalibrator,
    project_to_simplex, family_program_and_residual, contrast_weights,
)
from tissueresolve.experimental.soft_hierarchy import apply_partial_confidence_gating


# --------------------------------------------------------------------- fixtures
def _ref():
    """Toy fine ref: FamX (X1/X2 separable), FamP (P1 distinct, P2≈P3 collinear)."""
    rng = np.random.default_rng(0)
    G = 120
    genes = [f"g{i}" for i in range(G)]

    def prof(a, lvl=200.0):
        v = np.full(G, 2.0)
        v[list(a)] = lvl
        return v

    X1 = prof(range(0, 20))
    X2 = prof(range(20, 40))
    P1 = prof(range(40, 60))
    P2 = prof(range(60, 80))
    P3 = P2 * (1 + rng.normal(0, 0.004, G))
    R = np.vstack([X1, X2, P1, P2, P3]).astype(np.float32)
    cts = ["X1", "X2", "P1", "P2", "P3"]
    ref = ReferenceSignature(gene_names=genes, cell_types=cts, R_cpm=R,
                             R_log=np.log1p(R).astype(np.float32),
                             n_cells_per_type={c: 100 for c in cts})
    mapping = {"X1": "FamX", "X2": "FamX", "P1": "FamP", "P2": "FamP", "P3": "FamP"}
    return ref, mapping


def _query(ref, n=6, seed=1):
    """A small genes × samples query built from random subtype mixtures."""
    rng = np.random.default_rng(seed)
    R = ref.as_R_cpm()                       # (K, G)
    W = rng.dirichlet(np.ones(R.shape[0]), size=n)   # (n, K)
    mu = W @ R                                # (n, G)
    counts = rng.poisson(np.clip(mu, 0, None)).astype(float).T   # (G, n)
    return pd.DataFrame(counts, index=ref.gene_names,
                        columns=[f"s{i}" for i in range(n)])


def _raw_conditional(ref, mapping, query):
    from tissueresolve.api import deconv_bulk
    import warnings
    with warnings.catch_warnings():
        warnings.simplefilter("ignore")
        h = deconv_bulk(query, ref, solver="nnls", resolution_mode="hierarchical",
                        hierarchy_mapping=mapping)
    est = h.estimates
    return est.conditional_proportions, est


# --------------------------------------------------------------------- pure helpers
def test_project_to_simplex_rows_sum_to_one_nonneg():
    V = np.array([[0.2, -0.5, 1.3], [-1.0, -1.0, -1.0], [0.0, 0.0, 0.0]])
    W = project_to_simplex(V)
    assert np.all(W >= -1e-12)
    np.testing.assert_allclose(W.sum(axis=1), 1.0, atol=1e-9)
    # all-nonpositive / all-zero rows → uniform
    np.testing.assert_allclose(W[1], 1 / 3, atol=1e-9)
    np.testing.assert_allclose(W[2], 1 / 3, atol=1e-9)


def test_residual_contrast_is_centered():
    ref, _ = _ref()
    S = ref.as_R_cpm()[2:5]                  # FamP members
    omega = np.array([100.0, 100.0, 100.0])
    B, D = family_program_and_residual(S, omega)
    w = omega / omega.sum()
    # Σ_k ω_k D_{k,g} == 0 for every gene
    centered = (w[:, None] * D).sum(axis=0)
    np.testing.assert_allclose(centered, 0.0, atol=1e-6)


def test_contrast_weights_finite_nonnegative():
    ref, _ = _ref()
    S = ref.as_R_cpm()[2:5]
    w = contrast_weights(S, np.array([100.0, 100, 100]),
                         gene_names=list(ref.gene_names))
    assert np.all(np.isfinite(w))
    assert np.all(w >= 0.0)
    assert w.shape[0] == ref.as_R_cpm().shape[1]


def test_contrast_weights_handle_missing_donor_labels():
    ref, _ = _ref()
    S = ref.as_R_cpm()[2:5]
    # donor_residuals=None must be handled safely (donor stability neutral)
    w = contrast_weights(S, None, gene_names=list(ref.gene_names), donor_residuals=None)
    assert np.all(np.isfinite(w)) and np.all(w >= 0)


# --------------------------------------------------------------------- refiner invariants
def _refine(mode, **cfg_kw):
    ref, mapping = _ref()
    query = _query(ref)
    raw_cond, est = _raw_conditional(ref, mapping, query)
    cfg = FineRefinerConfig(mode=mode, **cfg_kw)
    r = FineGranularityRefiner(cfg).refine(raw_cond, ref, mapping, query=query)
    return r, raw_cond, est, mapping, ref


@pytest.mark.parametrize("mode", ["contrast_weighted", "residual_contrast"])
def test_refined_conditional_nonneg_and_sums_to_one(mode):
    r, raw_cond, est, mapping, ref = _refine(mode)
    cond = r.refined_conditional
    assert (cond.to_numpy() >= -1e-9).all()
    # within each family the conditional sums to 1
    fams: dict = {}
    for st in cond.columns:
        fams.setdefault(mapping[st], []).append(st)
    for fam, mem in fams.items():
        np.testing.assert_allclose(cond[mem].sum(axis=1), 1.0, atol=1e-6)
    assert r.metadata["within_family_mass_error"] < 1e-6


@pytest.mark.parametrize("mode", ["contrast_weighted", "residual_contrast"])
def test_refiner_does_not_change_broad_family_mass(mode):
    from tissueresolve.reference.hierarchy import combine_family_and_conditional_estimates
    r, raw_cond, est, mapping, ref = _refine(mode)
    fam_props = est.family_proportions
    abs_refined = combine_family_and_conditional_estimates(
        fam_props, r.refined_conditional, mapping)
    # per-family mass (sum of member subtypes) equals broad-family mass exactly
    fams: dict = {}
    for st in abs_refined.columns:
        fams.setdefault(mapping[st], []).append(st)
    for fam, mem in fams.items():
        np.testing.assert_allclose(abs_refined[mem].sum(axis=1), fam_props[fam], atol=1e-9)
    # total mass equals the family-level total
    np.testing.assert_allclose(abs_refined.sum(axis=1), fam_props.sum(axis=1), atol=1e-9)


def test_refiner_output_is_conditional_not_gated():
    # The refiner returns conditional proportions only — no unresolved columns,
    # i.e. soft gating has NOT been applied inside the refiner.
    r, *_ = _refine("contrast_weighted")
    assert not any(str(c).startswith("unresolved_") for c in r.refined_conditional.columns)


def test_order_refine_then_softgate_conserves_mass_once():
    from tissueresolve.reference.hierarchy import combine_family_and_conditional_estimates
    r, raw_cond, est, mapping, ref = _refine("contrast_weighted")
    fam_props = est.family_proportions
    # canonical order: refine → θ = π q_refined → soft gate (applied exactly once)
    abs_refined = combine_family_and_conditional_estimates(fam_props, r.refined_conditional, mapping)
    gated = apply_partial_confidence_gating(fam_props, abs_refined, mapping, confidence=0.6)
    assert gated.mass_conservation_error < 1e-6
    # total after gating equals broad total (resolved + unresolved)
    total = gated.combined.sum(axis=1)
    np.testing.assert_allclose(total, fam_props.sum(axis=1), atol=1e-6)


def test_default_mode_none_is_passthrough():
    ref, mapping = _ref()
    query = _query(ref)
    raw_cond, _ = _raw_conditional(ref, mapping, query)
    r = FineGranularityRefiner(FineRefinerConfig()).refine(raw_cond, ref, mapping, query=query)
    pd.testing.assert_frame_equal(r.refined_conditional, raw_cond.astype(float))
    assert r.metadata["fine_refinement"] == "none"
    assert r.refined_families == []


def test_insufficient_family_support_falls_back():
    ref, mapping = _ref()
    query = _query(ref)
    raw_cond, _ = _raw_conditional(ref, mapping, query)
    cfg = FineRefinerConfig(mode="contrast_weighted", min_family_donors=5)
    # claim FamP has too few donors → must be skipped and left identical to raw
    r = FineGranularityRefiner(cfg).refine(
        raw_cond, ref, mapping, query=query,
        family_donor_counts={"FamP": 2, "FamX": 10})
    assert "FamP" in r.skipped_families
    pd.testing.assert_series_equal(r.refined_conditional["P2"], raw_cond["P2"].astype(float),
                                   check_names=False)


def test_metadata_records_refinement_status():
    r, *_ = _refine("contrast_weighted")
    m = r.metadata
    assert m["fine_refinement"] == "contrast_weighted"
    assert m["version"].startswith("fine_refiner")
    assert m["feature_status"] == "experimental"
    assert isinstance(m["refined_families"], list)
    assert "within_family_mass_error" in m
    assert isinstance(m["contrast_components"], dict)


# --------------------------------------------------------------------- calibration
def _cond_frames(ref, mapping, members, n, seed):
    """Synthetic (pred, true) conditional frames for one family."""
    rng = np.random.default_rng(seed)
    true = pd.DataFrame(rng.dirichlet(np.ones(len(members)), size=n),
                        columns=members, index=[f"c{i}" for i in range(n)])
    # predicted = true blurred toward uniform (systematic spillover)
    pred = 0.6 * true + 0.4 / len(members)
    pred = pred.div(pred.sum(axis=1), axis=0)
    return pred, true


def test_calibration_projects_to_simplex_and_uses_only_fitted_data():
    members = ["P1", "P2", "P3"]
    pred_cal, true_cal = _cond_frames(None, None, members, 30, seed=2)
    cal = SpilloverCalibrator(kind="ridge", ridge_alpha=1.0).fit(
        {"FamP": pred_cal}, {"FamP": true_cal})
    assert "FamP" in cal.maps_
    # apply to *unseen* test rows
    pred_test, _ = _cond_frames(None, None, members, 8, seed=999)
    out = cal.transform_family("FamP", pred_test)
    assert (out.to_numpy() >= -1e-9).all()
    np.testing.assert_allclose(out.sum(axis=1), 1.0, atol=1e-9)
    # the fitted map must depend only on calibration data: refitting on the same
    # calibration data reproduces it exactly (no hidden state from transform).
    cal2 = SpilloverCalibrator(kind="ridge", ridge_alpha=1.0).fit(
        {"FamP": pred_cal}, {"FamP": true_cal})
    np.testing.assert_allclose(cal.maps_["FamP"][1], cal2.maps_["FamP"][1], atol=1e-12)


def test_calibration_none_is_identity():
    members = ["P1", "P2", "P3"]
    pred, _ = _cond_frames(None, None, members, 8, seed=3)
    cal = SpilloverCalibrator(kind="none").fit({"FamP": pred}, {"FamP": pred})
    out = cal.transform_family("FamP", pred)
    pd.testing.assert_frame_equal(out, pred)


def test_calibration_skips_small_families():
    members = ["P1", "P2", "P3"]
    pred, true = _cond_frames(None, None, members, 4, seed=4)   # only 4 rows < min_rows
    cal = SpilloverCalibrator(kind="full", min_rows=6).fit({"FamP": pred}, {"FamP": true})
    assert "FamP" in cal.skipped_ and "FamP" not in cal.maps_
    # transform falls back to identity for a skipped family
    pd.testing.assert_frame_equal(cal.transform_family("FamP", pred), pred)
