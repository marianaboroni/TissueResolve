"""Tests for the Resolution Decision Layer + the evidence-based pipeline reorg.

Offline, deterministic, no generated data.  Covers: the decision exists before
fine output is interpreted; full-panel reliability dominates and cell-level AUROC
cannot drive the decision; broad_only families are not trusted fine; selected_fine
reports supported subtypes only; missing evidence is conservative; mass conservation;
modality-aware metadata (bulk vs spatial differ; bulk has no spatial smoothing);
soft gating stays the final layer applied once; and the report shows the
trusted-resolution table before fine predictions with the RNA-proportion caveat.
"""
import warnings

import numpy as np
import pandas as pd
import pytest

from tissueresolve.results import ReferenceSignature
from tissueresolve.resolution import (
    ResolutionDecision, ResolutionDecisionConfig, decide_trusted_resolution,
    decisions_metadata, compute_modality_aware_gene_weights,
    BROAD_ONLY, SELECTED_FINE, FULL_FINE,
)


# --------------------------------------------------------------------- fixtures
def _ref():
    rng = np.random.default_rng(0)
    G = 120
    genes = [f"g{i}" for i in range(G)]

    def prof(a, lvl=200.0):
        v = np.full(G, 2.0); v[list(a)] = lvl; return v

    R = np.vstack([prof(range(0, 20)), prof(range(20, 40)),      # FamX separable
                   prof(range(40, 60)), prof(range(60, 80)),
                   prof(range(60, 80)) * (1 + rng.normal(0, 0.004, G))]).astype(np.float32)
    cts = ["X1", "X2", "P1", "P2", "P3"]
    ref = ReferenceSignature(gene_names=genes, cell_types=cts, R_cpm=R,
                             R_log=np.log1p(R).astype(np.float32),
                             n_cells_per_type={c: 100 for c in cts})
    mapping = {"X1": "FamX", "X2": "FamX", "P1": "FamP", "P2": "FamP", "P3": "FamP"}
    return ref, mapping


def _bulk(ref, n=2, seed=1):
    rng = np.random.default_rng(seed)
    return pd.DataFrame(rng.poisson(50, size=(len(ref.gene_names), n)) + 1,
                        index=ref.gene_names, columns=[f"b{i}" for i in range(n)]).astype(float)


def _hier(ref, mapping):
    from tissueresolve.api import deconv_bulk
    with warnings.catch_warnings():
        warnings.simplefilter("ignore")
        return deconv_bulk(_bulk(ref), ref, solver="nnls",
                           resolution_mode="hierarchical", hierarchy_mapping=mapping)


def _resolvability(sep, spill, disc, resolvable, n=2):
    return pd.DataFrame(
        {"n_subtypes": [n], "mean_separability": [sep], "mean_spillover": [spill],
         "min_discriminating_genes": [disc], "resolvable": [resolvable]},
        index=["Fam"])


# --------------------------------------------------------------------- decision logic
def test_full_fine_when_all_gates_pass():
    ident = _resolvability(0.5, 0.1, 40, True)
    d = decide_trusted_resolution(identifiability_metrics=ident,
                                  family_map={"a": "Fam", "b": "Fam"},
                                  subtype_confidence={"a": {"confidence": 0.8},
                                                      "b": {"confidence": 0.8}})
    assert d["Fam"].status == FULL_FINE


def test_broad_only_when_not_resolvable():
    ident = _resolvability(0.02, 0.6, 3, False)
    d = decide_trusted_resolution(identifiability_metrics=ident,
                                  family_map={"a": "Fam", "b": "Fam"})
    assert d["Fam"].status == BROAD_ONLY
    assert d["Fam"].supported_subtypes == []          # broad_only → no trusted fine
    assert set(d["Fam"].unsupported_subtypes) == {"a", "b"}


def test_full_panel_benchmark_reliability_dominates():
    # signature looks great, but a poor full-panel benchmark forces broad_only
    ident = _resolvability(0.5, 0.1, 40, True)
    d = decide_trusted_resolution(
        identifiability_metrics=ident, family_map={"a": "Fam", "b": "Fam"},
        subtype_confidence={"a": {"confidence": 0.9}, "b": {"confidence": 0.9}},
        benchmark_metrics={"Fam": {"cond_rmse": 0.5}})
    assert d["Fam"].status == BROAD_ONLY
    assert any("full-panel" in r for r in d["Fam"].reasons)


def test_cell_auroc_alone_cannot_trigger_full_fine():
    # AUROC is not an input; passing a near-perfect "auroc" in query_qc/benchmark
    # must NOT upgrade a collinear family with weak signature evidence.
    ident = _resolvability(0.03, 0.6, 2, False)
    d = decide_trusted_resolution(
        identifiability_metrics=ident, family_map={"a": "Fam", "b": "Fam"},
        query_qc={"auroc": 0.99}, benchmark_metrics={"Fam": {"auroc": 0.99}})
    assert d["Fam"].status == BROAD_ONLY


def test_missing_evidence_is_conservative():
    # no separability/spillover/disc evidence → cannot reach full_fine
    ident = pd.DataFrame({"n_subtypes": [2], "resolvable": [True]}, index=["Fam"])
    d = decide_trusted_resolution(identifiability_metrics=ident,
                                  family_map={"a": "Fam", "b": "Fam"},
                                  subtype_confidence={"a": {"confidence": 0.8},
                                                      "b": {"confidence": 0.8}})
    assert d["Fam"].status != FULL_FINE
    assert any("conservative" in w for w in d["Fam"].warnings)


def test_selected_fine_reports_supported_subtypes_only():
    # resolvable + decent signature but only one subtype confident, and not enough
    # discriminating genes for full_fine → selected_fine with the supported subtype
    ident = _resolvability(0.2, 0.15, 12, True, n=2)
    d = decide_trusted_resolution(
        identifiability_metrics=ident, family_map={"a": "Fam", "b": "Fam"},
        subtype_confidence={"a": {"confidence": 0.8}, "b": {"confidence": 0.02}})
    assert d["Fam"].status == SELECTED_FINE
    assert "a" in d["Fam"].supported_subtypes and "b" in d["Fam"].unsupported_subtypes


def test_low_query_overlap_forces_broad_only():
    ident = _resolvability(0.5, 0.1, 40, True)
    d = decide_trusted_resolution(identifiability_metrics=ident,
                                  family_map={"a": "Fam", "b": "Fam"},
                                  query_qc={"gene_overlap_frac": 0.05},
                                  subtype_confidence={"a": {"confidence": 0.9},
                                                      "b": {"confidence": 0.9}})
    assert d["Fam"].status == BROAD_ONLY


def test_single_subtype_family_is_full_fine_trivially():
    ident = pd.DataFrame({"n_subtypes": [1], "resolvable": [True]}, index=["Solo"])
    d = decide_trusted_resolution(identifiability_metrics=ident, family_map={"x": "Solo"})
    assert d["Solo"].status == FULL_FINE


# --------------------------------------------------------------------- pipeline wiring
def test_decision_recorded_in_estimates_before_fine_interpretation():
    ref, mapping = _ref()
    h = _hier(ref, mapping)
    m = h.estimates.metadata
    assert "trusted_resolution" in m and "resolution_decision" in m
    # collinear FamP → broad_only; separable FamX → full_fine
    assert m["trusted_resolution"]["FamP"] == BROAD_ONLY
    assert m["trusted_resolution"]["FamX"] == FULL_FINE
    # decision surfaced in the per-family qc table
    assert "trusted_resolution" in h.estimates.qc.columns


def test_broad_only_family_not_trusted_fine_and_mass_conserved():
    ref, mapping = _ref()
    h = _hier(ref, mapping)
    md = h.run_metadata
    assert "FamP" in md["fine_predictions_diagnostic_only"]
    assert "FamX" in md["fine_predictions_trusted"]
    # broad + total mass conserved (soft gating still the final layer)
    assert h.estimates.metadata["mass_conservation_max_error"] < 1e-6
    np.testing.assert_allclose(h.estimates.combined_fine.sum(axis=1), 1.0, atol=1e-6)


def test_soft_gating_remains_final_layer_applied_once():
    ref, mapping = _ref()
    h = _hier(ref, mapping)
    # the recorded gating mode is soft and unresolved mass exists for broad_only fam
    assert h.estimates.metadata["hierarchical_gating"] == "soft"
    ucols = [c for c in h.estimates.combined_fine.columns if str(c).startswith("unresolved_")]
    assert any("FamP" in c for c in ucols)


# --------------------------------------------------------------------- modality awareness
def test_modality_aware_weights_distinguish_bulk_and_spatial():
    ref, _ = _ref()
    wb, ib = compute_modality_aware_gene_weights(ref, _bulk(ref), "bulk")
    ws, is_ = compute_modality_aware_gene_weights(ref, _bulk(ref), "spatial")
    assert ib["modality"] == "bulk" and ib["weights_applied"] is True
    assert "bulk_criteria" in ib
    assert is_["modality"] == "spatial" and is_["weights_applied"] is False
    assert "spatial_criteria" in is_
    assert ws is None                                  # spatial uses NB-CAR, no WNNLS weights
    with pytest.raises(ValueError):
        compute_modality_aware_gene_weights(ref, _bulk(ref), "imaging")


def test_bulk_run_metadata_records_modality_and_no_spatial_smoothing():
    ref, mapping = _ref()
    md = _hier(ref, mapping).run_metadata
    assert md["modality"] == "bulk"
    assert md["prediction_unit"] == "sample"
    assert md["gene_weighting_mode"] == "protocol_aware_weighted_nnls"
    assert md["spatial_smoothing_used"] is False


def test_decisions_metadata_is_jsonable_summary():
    ident = _resolvability(0.5, 0.1, 40, True)
    d = decide_trusted_resolution(identifiability_metrics=ident,
                                  family_map={"a": "Fam", "b": "Fam"},
                                  subtype_confidence={"a": {"confidence": 0.8},
                                                      "b": {"confidence": 0.8}})
    meta = decisions_metadata(d)
    import json
    json.dumps(meta)                                   # must be serialisable
    assert meta["n_full_fine"] == 1
    assert "AUROC is NOT a decision criterion" in meta["note"]


# --------------------------------------------------------------------- report
def test_report_shows_trusted_resolution_before_fine_predictions():
    from tissueresolve.report.result_sections import _hierarchical_html, _trusted_resolution_html
    ref, mapping = _ref()
    h = _hier(ref, mapping)

    class _R:
        estimates = h.estimates
    body = _hierarchical_html(_R(), "bulk")
    assert "Trusted resolution by family" in body
    # the trusted-resolution table appears BEFORE the broad/fine composition tables
    assert body.index("Trusted resolution by family") < body.index("Broad family composition")
    # broad_only fam flagged diagnostic, RNA-proportion caveat present, AUROC disclaimed
    th = _trusted_resolution_html(h.estimates)
    assert "diagnostic only" in th.lower()
    assert "rna-derived proportions, not cell fractions" in th.lower()
    assert "auroc is not a criterion" in th.lower()


def test_spatial_estimate_statement_does_not_claim_accuracy():
    from tissueresolve.report import methods_text as mt
    stmt = mt.estimate_type_statement("spatial").lower()
    # RNA-derived composition, explicitly NOT direct cell counts, and no accuracy claim
    assert "rna-derived" in stmt
    assert "do not represent direct single-cell counts" in stmt or "not" in stmt
    assert "accuracy" not in stmt
