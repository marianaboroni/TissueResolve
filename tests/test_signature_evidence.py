"""Offline tests for independent sibling / rare-confirmation evidence scores (Stage 1.5)."""
from __future__ import annotations

import numpy as np
import pandas as pd

from tissueresolve.reference.signature_evidence import (
    to_log_cpm, rare_confirmation_evidence, sibling_evidence)


def _query():
    """Two bulk samples (genes×samples): s_A dominated by subtype-A genes (0-9),
    s_B dominated by subtype-B genes (10-19). Background genes 20-59 low."""
    G = 60
    genes = [f"g{i}" for i in range(G)]
    a = np.full(G, 2.0); a[0:10] = 400.0
    b = np.full(G, 2.0); b[10:20] = 400.0
    return pd.DataFrame({"s_A": a, "s_B": b}, index=genes)


def test_rare_confirmation_evidence_specific():
    q = _query()
    rc = {"A": [f"g{i}" for i in range(10)], "B": [f"g{i}" for i in range(10, 20)]}
    ev = rare_confirmation_evidence(q, rc)
    assert ev.loc["s_A", "A"] > 0.8 and ev.loc["s_A", "B"] < 0.3   # A sample -> A genes on
    assert ev.loc["s_B", "B"] > 0.8 and ev.loc["s_B", "A"] < 0.3
    assert (ev.to_numpy() >= 0).all() and (ev.to_numpy() <= 1).all()


def test_sibling_evidence_direction():
    q = _query()
    sib = {"A": [f"g{i}" for i in range(10)], "B": [f"g{i}" for i in range(10, 20)]}
    mapping = {"A": "fam1", "B": "fam1"}
    ev = sibling_evidence(q, sib, mapping)
    assert ev.loc["s_A", "A"] > 0 > ev.loc["s_A", "B"]            # A sample favours A over B
    assert ev.loc["s_B", "B"] > 0 > ev.loc["s_B", "A"]


def test_deterministic_and_missing_genes():
    q = _query()
    rc = {"A": ["g0", "g1", "gZZZ"]}          # one gene absent from query
    a = rare_confirmation_evidence(q, rc); b = rare_confirmation_evidence(q, rc)
    pd.testing.assert_frame_equal(a, b)        # deterministic
    assert np.isfinite(a.loc["s_A", "A"])      # absent gene ignored, not crash
    empty = rare_confirmation_evidence(q, {"X": ["gNONE"]})
    assert empty["X"].isna().all()             # no observable genes -> NaN, not fake 0


def test_log_cpm_shape():
    q = _query()
    lc = to_log_cpm(q)
    assert lc.shape == q.shape and (lc.to_numpy() >= 0).all()
