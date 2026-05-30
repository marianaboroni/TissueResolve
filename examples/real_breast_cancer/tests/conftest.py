"""
Offline test fixtures for the real-data validation harness.

These tests use tiny synthetic AnnData and never access the network or download
any dataset.  The scripts directory is put on sys.path so ``import _harness``
and loading the numbered scripts works.
"""
from __future__ import annotations

import importlib.util
import sys
from pathlib import Path

import numpy as np
import pytest

SCRIPTS_DIR = Path(__file__).resolve().parent.parent / "scripts"
if str(SCRIPTS_DIR) not in sys.path:
    sys.path.insert(0, str(SCRIPTS_DIR))


@pytest.fixture
def harness():
    import _harness

    return _harness


@pytest.fixture
def load_script():
    """Return a loader for the digit-prefixed numbered scripts."""
    def _load(filename: str):
        safe = "rbc_" + filename.replace(".py", "").replace("-", "_")
        spec = importlib.util.spec_from_file_location(safe, SCRIPTS_DIR / filename)
        mod = importlib.util.module_from_spec(spec)
        spec.loader.exec_module(mod)
        return mod

    return _load


@pytest.fixture
def tiny_sc_adata():
    """A small single-cell AnnData with cell-type labels and integer counts.

    3 well-populated cell types (60 cells each) plus a rare type (5 cells) so
    the < min_cells drop can be exercised.  40 genes, block-structured.
    """
    import anndata as ad
    import pandas as pd

    rng = np.random.default_rng(0)
    types = ["Tcell"] * 60 + ["Bcell"] * 60 + ["Myeloid"] * 60 + ["Rare"] * 5
    n_cells = len(types)
    n_genes = 40
    block = n_genes // 3

    X = rng.poisson(0.5, (n_cells, n_genes)).astype(np.float32)
    type_to_block = {"Tcell": 0, "Bcell": 1, "Myeloid": 2, "Rare": 0}
    for i, t in enumerate(types):
        b = type_to_block[t]
        X[i, b * block:(b + 1) * block] += rng.poisson(8, block)

    obs = pd.DataFrame({"cell_type": types},
                       index=[f"cell_{i}" for i in range(n_cells)])
    var = pd.DataFrame(index=[f"GENE_{i:03d}" for i in range(n_genes)])
    adata = ad.AnnData(X=X, obs=obs, var=var)
    adata.layers["counts"] = X.copy()
    return adata
