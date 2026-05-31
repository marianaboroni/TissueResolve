"""Offline tests for reference suitability and spatial auto-parameter selection."""
from __future__ import annotations

import warnings

import numpy as np
import pandas as pd
import pytest

from tissueresolve.results import ReferenceSignature


def _ref(separable=True, n_donor_cv=True):
    rng = np.random.default_rng(0)
    G = 80
    genes = [f"g{i}" for i in range(G)]

    def prof(a, lvl=200.0):
        v = np.full(G, 2.0); v[list(a)] = lvl; return v
    if separable:
        R = np.vstack([prof(range(0, 20)), prof(range(20, 40)),
                       prof(range(40, 60)), prof(range(60, 80))]).astype(np.float32)
    else:
        base = prof(range(0, 20))
        R = np.vstack([base, base * 1.001, base * 0.999, base * 1.002]).astype(np.float32)
    cts = ["A", "B", "C", "D"]
    return ReferenceSignature(
        gene_names=genes, cell_types=cts, R_cpm=R,
        R_log=np.log1p(R).astype(np.float32),
        donor_cv=(rng.uniform(0.05, 0.2, (G, 4)).astype(np.float32) if n_donor_cv else None),
        n_cells_per_type={c: 200 for c in cts})


# --- reference suitability --------------------------------------------------

def test_suitability_high_quality_passes():
    from tissueresolve.reference import suitability as su
    ref = _ref(separable=True)
    with warnings.catch_warnings():
        warnings.simplefilter("ignore")
        res = su.compute_reference_suitability_score(
            ref, query_genes=list(ref.gene_names),
            mapping={"A": "F1", "B": "F1", "C": "F2", "D": "F2"},
            protocol_info={"risk": "low"},
            library_info={"overall": "scrna_10x_3p"})
    assert res.classification in ("PASS", "CAUTION")
    assert res.overall_score >= 0.6


def test_suitability_poor_gene_overlap_fails():
    from tissueresolve.reference import suitability as su
    ref = _ref()
    c = su.score_gene_overlap(ref, query_genes=["zzz1", "zzz2"])  # ~no overlap
    assert c.status in ("WARNING", "FAIL")
    with warnings.catch_warnings():
        warnings.simplefilter("ignore")
        res = su.compute_reference_suitability_score(ref, query_genes=["zzz1"])
    assert res.classification in ("WARNING", "FAIL")


def test_suitability_missing_hierarchy_is_unknown():
    from tissueresolve.reference import suitability as su
    ref = _ref()
    c = su.score_hierarchy_quality(ref, mapping=None)
    assert c.status == "UNKNOWN"


def test_suitability_failing_component_caps_classification():
    """A FAILing component must not be hidden behind a high mean."""
    from tissueresolve.reference import suitability as su
    # imbalanced reference (one tiny type) → balance FAIL
    ref = _ref()
    ref.n_cells_per_type = {"A": 1000, "B": 1000, "C": 1000, "D": 2}
    with warnings.catch_warnings():
        warnings.simplefilter("ignore")
        res = su.compute_reference_suitability_score(
            ref, query_genes=list(ref.gene_names),
            mapping={"A": "F1", "B": "F1", "C": "F2", "D": "F2"})
    bal = next(c for c in res.components if c.name == "celltype_balance")
    assert bal.status == "FAIL"
    assert res.classification != "PASS"


def test_suitability_save(tmp_path):
    from tissueresolve.reference import suitability as su
    ref = _ref()
    with warnings.catch_warnings():
        warnings.simplefilter("ignore")
        res = su.compute_reference_suitability_score(ref, query_genes=list(ref.gene_names))
        written = su.save_reference_suitability(res, tmp_path)
    for k in ("score", "components", "report", "warnings"):
        assert written[k].exists()


# --- spatial auto-parameter selection ---------------------------------------

def test_over_smoothing_penalty_monotone():
    from tissueresolve.spatial.auto_params import over_smoothing_penalty
    base = pd.DataFrame({"A": [0.1, 0.9, 0.5], "B": [0.9, 0.1, 0.5]})
    flat = pd.DataFrame({"A": [0.5, 0.5, 0.5], "B": [0.5, 0.5, 0.5]})
    assert over_smoothing_penalty(base, base) == pytest.approx(0.0, abs=1e-9)
    assert over_smoothing_penalty(flat, base) > 0.9  # fully smoothed → high penalty


def test_select_lambda_spatial_runs_and_records(tmp_path):
    from benchmarks.shared.synthetic import toy_reference, toy_spatial
    from tissueresolve.spatial.auto_params import (
        select_lambda_spatial, save_spatial_auto_outputs)
    ref, _ = toy_reference()
    sc, _ = toy_spatial(ref, n_rows=6, n_cols=6)
    with warnings.catch_warnings():
        warnings.simplefilter("ignore")
        best, comp, sel = select_lambda_spatial(
            sc["Y"], ref, sc["array_row"], sc["array_col"], sc["lib_sizes"],
            sc["gene_names"], sc["spot_ids"], candidates=(0.0, 0.5), max_iter=6)
    assert best in (0.0, 0.5)
    assert "objective" in comp.columns and "over_smoothing_penalty" in comp.columns
    assert "selected_lambda_spatial" in sel
    written = save_spatial_auto_outputs(comp, sel, tmp_path)
    assert written["comparison"].exists() and written["selected"].exists()


def test_graph_diagnostics():
    from benchmarks.shared.synthetic import toy_reference, toy_spatial
    from tissueresolve.spatial.auto_params import graph_parameter_diagnostics
    ref, _ = toy_reference()
    sc, _ = toy_spatial(ref, n_rows=5, n_cols=5)
    d = graph_parameter_diagnostics(sc["array_row"], sc["array_col"])
    assert d["n_spots"] == 25
