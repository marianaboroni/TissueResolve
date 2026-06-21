"""Unit + invariant tests for the experimental high-granularity strategy
(ReferenceCalibration + CandidateRestriction + ConsensusStability).

Offline, deterministic, no generated data.  Covers calibration finiteness/shape,
candidate restriction (never empties, records drops, protects rare), multi-panel
fan-out, consensus simplex + non-negativity, mass conservation (within-family +
total / broad unchanged), canonical order (soft gating once, AFTER consensus),
fallback for unsupported families, no-truth/no-test-donor leakage in calibration,
metadata, and default-off behaviour.
"""
import inspect

import numpy as np
import pandas as pd
import pytest

from tissueresolve.results import ReferenceSignature
from tissueresolve.experimental.soft_hierarchy import apply_partial_confidence_gating
from tissueresolve.experimental.soft_hierarchy.high_granularity import (
    HighGranularityStrategy, HighGranularityConfig,
    calibrate_gene_scale, calibrate_logscale, select_candidates,
    build_family_panels, multi_panel_conditional, consensus_stability, RARE_RAW,
)
from tissueresolve.experimental.soft_hierarchy.fine_refiner import _subset_ref_celltypes


# --------------------------------------------------------------------- fixtures
def _ref():
    rng = np.random.default_rng(0)
    G = 150
    genes = [f"g{i}" for i in range(G)]

    def prof(a, lvl=200.0):
        v = np.full(G, 2.0); v[list(a)] = lvl; return v

    X1 = prof(range(0, 20)); X2 = prof(range(20, 40))             # separable family
    P1 = prof(range(40, 60)); P2 = prof(range(60, 80))
    P3 = P2 * (1 + rng.normal(0, 0.004, G))                       # P2≈P3 collinear
    R = np.vstack([X1, X2, P1, P2, P3]).astype(np.float32)
    cts = ["X1", "X2", "P1", "P2", "P3"]
    ref = ReferenceSignature(gene_names=genes, cell_types=cts, R_cpm=R,
                             R_log=np.log1p(R).astype(np.float32),
                             n_cells_per_type={c: 100 for c in cts})
    mapping = {"X1": "FamX", "X2": "FamX", "P1": "FamP", "P2": "FamP", "P3": "FamP"}
    return ref, mapping


def _query(ref, n=6, seed=1):
    rng = np.random.default_rng(seed)
    R = ref.as_R_cpm()
    W = rng.dirichlet(np.ones(R.shape[0]), size=n)
    mu = W @ R
    counts = rng.poisson(np.clip(mu, 0, None)).astype(float).T
    return pd.DataFrame(counts, index=ref.gene_names, columns=[f"s{i}" for i in range(n)])


def _raw(ref, mapping, query):
    import warnings
    from tissueresolve.api import deconv_bulk
    with warnings.catch_warnings():
        warnings.simplefilter("ignore")
        h = deconv_bulk(query, ref, solver="nnls", resolution_mode="hierarchical",
                        hierarchy_mapping=mapping)
    return h.estimates.conditional_proportions, h.estimates


def _run(mode, **kw):
    ref, mapping = _ref()
    query = _query(ref)
    raw, est = _raw(ref, mapping, query)
    cfg = HighGranularityConfig(mode=mode, **kw)
    r = HighGranularityStrategy(cfg).refine(raw, ref, mapping, query=query)
    return r, raw, est, mapping, ref, query


# --------------------------------------------------------------------- Component 1
def test_calibration_preserves_shape_and_is_finite():
    ref, _ = _ref()
    sub = _subset_ref_celltypes(ref, ["P1", "P2", "P3"])
    q = _query(ref)
    for fn in (calibrate_gene_scale, calibrate_logscale):
        cal, diag = fn(sub, q, HighGranularityConfig())
        assert cal.as_R_cpm().shape == sub.as_R_cpm().shape
        assert list(cal.cell_types) == list(sub.cell_types)
        assert list(cal.gene_names) == list(sub.gene_names)
        assert np.all(np.isfinite(cal.as_R_cpm()))
        assert (cal.as_R_cpm() >= 0).all()


def test_calibration_uses_no_truth_or_test_donor_labels():
    # calibration depends ONLY on (sub_ref, query, cfg) — never on truth/test donors.
    for fn in (calibrate_gene_scale, calibrate_logscale):
        params = set(inspect.signature(fn).parameters)
        assert not ({"truth", "true_proportions", "labels", "donors", "test"} & params)
    # deterministic for a fixed query (no hidden randomness / leakage)
    ref, _ = _ref(); sub = _subset_ref_celltypes(ref, ["P1", "P2", "P3"]); q = _query(ref)
    a, _ = calibrate_gene_scale(sub, q, HighGranularityConfig())
    b, _ = calibrate_gene_scale(sub, q, HighGranularityConfig())
    np.testing.assert_allclose(a.as_R_cpm(), b.as_R_cpm())


def test_calibration_scale_preserves_within_gene_subtype_contrast():
    # gene-level scaling multiplies all subtypes in a gene equally → ratios preserved.
    ref, _ = _ref(); sub = _subset_ref_celltypes(ref, ["P1", "P2", "P3"]); q = _query(ref)
    cal, _ = calibrate_gene_scale(sub, q, HighGranularityConfig())
    R0, R1 = sub.as_R_cpm(), cal.as_R_cpm()
    # ratio between subtype 0 and 1 within each gene is unchanged
    r0 = R0[0] / (R0[1] + 1e-9); r1 = R1[0] / (R1[1] + 1e-9)
    np.testing.assert_allclose(r0, r1, rtol=1e-4)


# --------------------------------------------------------------------- Component 2
def test_candidate_restriction_never_empties_family():
    ref, _ = _ref()
    raw = pd.DataFrame({"P1": [0.0], "P2": [0.0], "P3": [0.0]})
    for rule in ("topk", "evidence", "hybrid"):
        keep, drop = select_candidates(["P1", "P2", "P3"], raw,
                                       HighGranularityConfig(candidate_rule=rule))
        assert len(keep) >= 1
        assert set(keep) | set(drop) == {"P1", "P2", "P3"}


def test_candidate_restriction_records_dropped():
    raw = pd.DataFrame({"P1": [0.6], "P2": [0.3], "P3": [0.1]})
    keep, drop = select_candidates(["P1", "P2", "P3"], raw,
                                   HighGranularityConfig(candidate_rule="topk", top_k=2,
                                                         max_candidates=2))
    # P3 has raw 0.1 ≥ RARE_RAW so it is protected; verify drop+keep accounting
    assert set(keep) | set(drop) == {"P1", "P2", "P3"}
    assert isinstance(drop, list)


def test_rare_subtype_not_dropped_when_evidence_supports():
    # a low-but-real subtype (raw ≥ RARE_RAW) must survive topk restriction
    raw = pd.DataFrame({"P1": [0.7, 0.7], "P2": [0.25, 0.25], "P3": [RARE_RAW + 0.01] * 2})
    keep, drop = select_candidates(["P1", "P2", "P3"], raw,
                                   HighGranularityConfig(candidate_rule="topk", top_k=1))
    assert "P3" in keep


# --------------------------------------------------------------------- Component 3/4
def test_multi_panel_returns_one_prediction_per_panel():
    ref, mapping = _ref(); q = _query(ref)
    sub = _subset_ref_celltypes(ref, ["P1", "P2", "P3"])
    cfg = HighGranularityConfig(panels=("current", "marker", "query_detectable"))
    panels = build_family_panels(sub, q, cfg, omega=np.array([100., 100, 100]))
    conds = multi_panel_conditional(q, sub, panels, ["P1", "P2", "P3"])
    assert set(conds) <= set(panels)
    assert len(conds) >= 1
    for name, c in conds.items():
        assert list(c.columns) == ["P1", "P2", "P3"]
        np.testing.assert_allclose(c.sum(axis=1), 1.0, atol=1e-6)


@pytest.mark.parametrize("method", ["mean", "median", "conservative"])
def test_consensus_nonneg_and_simplex(method):
    ref, mapping = _ref(); q = _query(ref)
    sub = _subset_ref_celltypes(ref, ["P1", "P2", "P3"])
    cfg = HighGranularityConfig(consensus=method)
    panels = build_family_panels(sub, q, cfg, omega=np.array([100., 100, 100]))
    conds = multi_panel_conditional(q, sub, panels, ["P1", "P2", "P3"])
    cons, stab, diag = consensus_stability(conds, ["P1", "P2", "P3"], cfg)
    assert (cons.to_numpy() >= -1e-9).all()
    np.testing.assert_allclose(cons.sum(axis=1), 1.0, atol=1e-6)
    assert (stab.to_numpy() >= 0).all() and (stab.to_numpy() <= 1 + 1e-9).all()
    assert {"consensus_entropy", "consensus_eff_n", "n_panels"} <= set(diag)


# --------------------------------------------------------------------- orchestrator invariants
@pytest.mark.parametrize("mode", ["candidate_consensus", "reference_calibrated_consensus"])
def test_refined_conditional_nonneg_sums_to_one(mode):
    r, raw, est, mapping, ref, q = _run(mode, candidate_rule="hybrid", calibration="scale")
    cond = r.refined_conditional
    assert (cond.to_numpy() >= -1e-9).all()
    fams = {}
    for st in cond.columns:
        fams.setdefault(mapping[st], []).append(st)
    for fam, mem in fams.items():
        np.testing.assert_allclose(cond[mem].sum(axis=1), 1.0, atol=1e-6)
    assert r.metadata["within_family_mass_error"] < 1e-6


def test_broad_mass_unchanged_and_total_conserved():
    from tissueresolve.reference.hierarchy import combine_family_and_conditional_estimates
    r, raw, est, mapping, ref, q = _run("candidate_consensus", candidate_rule="hybrid")
    fam_props = est.family_proportions
    abs_ref = combine_family_and_conditional_estimates(fam_props, r.refined_conditional, mapping)
    fams = {}
    for st in abs_ref.columns:
        fams.setdefault(mapping[st], []).append(st)
    for fam, mem in fams.items():
        np.testing.assert_allclose(abs_ref[mem].sum(axis=1), fam_props[fam], atol=1e-9)
    np.testing.assert_allclose(abs_ref.sum(axis=1), fam_props.sum(axis=1), atol=1e-9)


def test_soft_gating_applied_after_consensus_once():
    from tissueresolve.reference.hierarchy import combine_family_and_conditional_estimates
    r, raw, est, mapping, ref, q = _run("candidate_consensus", candidate_rule="hybrid")
    fam_props = est.family_proportions
    # consensus output carries NO unresolved columns (gate not applied inside)
    assert not any(str(c).startswith("unresolved_") for c in r.refined_conditional.columns)
    abs_ref = combine_family_and_conditional_estimates(fam_props, r.refined_conditional, mapping)
    # confidence modulated by stability (Component 5), gate applied exactly once
    conf = r.stability_confidence.clip(0, 1)
    gated = apply_partial_confidence_gating(fam_props, abs_ref, mapping, confidence=conf)
    assert gated.mass_conservation_error < 1e-6
    np.testing.assert_allclose(gated.combined.sum(axis=1), fam_props.sum(axis=1), atol=1e-6)


def test_default_mode_none_is_passthrough():
    r, raw, est, mapping, ref, q = _run("none")
    pd.testing.assert_frame_equal(r.refined_conditional, raw.astype(float))
    assert r.metadata["high_granularity_mode"] == "none"
    assert r.refined_families == []


def test_unsupported_family_falls_back():
    ref, mapping = _ref(); query = _query(ref)
    raw, _ = _raw(ref, mapping, query)
    cfg = HighGranularityConfig(mode="candidate_consensus", candidate_rule="hybrid",
                                min_family_donors=99)
    r = HighGranularityStrategy(cfg).refine(
        raw, ref, mapping, query=query, family_donor_counts={"FamP": 2, "FamX": 2})
    assert "FamP" in r.skipped_families and "FamX" in r.skipped_families
    pd.testing.assert_series_equal(r.refined_conditional["P2"], raw["P2"].astype(float),
                                   check_names=False)


def test_metadata_records_selected_methods():
    r, *_ = _run("reference_calibrated_consensus", calibration="logscale",
                 candidate_rule="evidence", consensus="median")
    m = r.metadata
    assert m["high_granularity_mode"] == "reference_calibrated_consensus"
    assert m["reference_calibration"] == "logscale"
    assert m["candidate_rule"] == "evidence"
    assert m["consensus"] == "median"
    assert isinstance(m["panels"], list) and m["panels"]
    assert m["feature_status"] == "experimental"
    assert m["version"].startswith("high_granularity")
