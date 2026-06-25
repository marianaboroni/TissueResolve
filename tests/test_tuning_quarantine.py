"""
Regression test for the quarantined tuning stub.

The previous ``tissueresolve.tuning`` stub wrote *fabricated* tuning grids and
scores (e.g. ``marker_genes=50 -> 0.8``) to disk without running any search.
It must never silently return invented metrics again: the functions now raise
``NotImplementedError`` and must not write any artifacts.
"""
from __future__ import annotations

import pytest

from tissueresolve import tuning


def test_tune_bulk_parameters_raises_and_writes_nothing(tmp_path):
    out = tmp_path / "bulk_tuning"
    with pytest.raises(NotImplementedError):
        tuning.tune_bulk_parameters(out)
    # no fabricated artifacts written
    assert not out.exists() or not any(out.iterdir())


def test_tune_spatial_parameters_raises_and_writes_nothing(tmp_path):
    out = tmp_path / "spatial_tuning"
    with pytest.raises(NotImplementedError):
        tuning.tune_spatial_parameters(out)
    assert not out.exists() or not any(out.iterdir())
