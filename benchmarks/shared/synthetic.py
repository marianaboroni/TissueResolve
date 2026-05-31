"""
Synthetic, offline benchmark data with known ground truth.

Builds a toy reference with broad families (each containing separable and
near-identical subtypes), pseudobulk mixtures, and a small spatial grid — used
by the toy benchmark and the offline tests.  No network, fully deterministic.
"""
from __future__ import annotations

import numpy as np
import pandas as pd

from tissueresolve.results import ReferenceSignature


def toy_reference(*, n_genes: int = 120, seed: int = 0):
    """Reference: FamX (2 separable), FamY (2 near-identical), FamZ (singleton)."""
    rng = np.random.default_rng(seed)
    genes = [f"g{i}" for i in range(n_genes)]

    def prof(active, lvl=200.0):
        v = np.full(n_genes, 2.0)
        v[list(active)] = lvl
        return v

    X1 = prof(range(0, 20)); X2 = prof(range(20, 40))
    Y1 = prof(range(40, 60)); Y2 = Y1 * (1.0 + rng.normal(0, 0.004, n_genes))
    Z = prof(range(60, 80))
    R = np.vstack([X1, X2, Y1, Y2, Z]).astype(np.float32)
    cts = ["X_sub1", "X_sub2", "Y_sub1", "Y_sub2", "Z_solo"]
    ref = ReferenceSignature(
        gene_names=genes, cell_types=cts, R_cpm=R,
        R_log=np.log1p(R).astype(np.float32),
        phi_g=np.full(n_genes, 5.0, dtype=np.float32),
        donor_cv=rng.uniform(0.05, 0.4, size=(n_genes, 5)).astype(np.float32),
        n_cells_per_type={c: 100 for c in cts})
    mapping = {"X_sub1": "FamX", "X_sub2": "FamX",
               "Y_sub1": "FamY", "Y_sub2": "FamY", "Z_solo": "FamZ"}
    return ref, mapping


def toy_bulk(ref, *, n_samples: int = 6, seed: int = 1):
    """Pseudobulk mixtures + ground-truth proportions (samples × cell types)."""
    rng = np.random.default_rng(seed)
    R = ref.as_R_cpm()
    true = rng.dirichlet(np.ones(R.shape[0]), size=n_samples)
    mat = true @ R
    bulk = pd.DataFrame(mat.T, index=list(ref.gene_names),
                        columns=[f"s{i}" for i in range(n_samples)])
    truth = pd.DataFrame(true, index=[f"s{i}" for i in range(n_samples)],
                         columns=list(ref.cell_types))
    return bulk, truth


def toy_spatial(ref, *, n_rows: int = 6, n_cols: int = 6, seed: int = 2):
    """Synthetic spatial grid with Poisson counts + ground-truth proportions."""
    rng = np.random.default_rng(seed)
    n = n_rows * n_cols
    R = ref.as_R_cpm()
    true = rng.dirichlet(np.ones(R.shape[0]), size=n)
    lib = rng.integers(800, 1500, size=n).astype("float32")
    lam = true @ R
    lam = lam / lam.sum(axis=1, keepdims=True)
    Y = rng.poisson(np.clip(lam * lib[:, None], 0, None)).astype("float32")
    rows = np.repeat(np.arange(n_rows), n_cols)
    cols = np.tile(np.arange(n_cols), n_rows)
    spot_ids = [f"spot{i}" for i in range(n)]
    truth = pd.DataFrame(true, index=spot_ids, columns=list(ref.cell_types))
    return {"Y": Y, "array_row": rows, "array_col": cols, "lib_sizes": lib,
            "gene_names": list(ref.gene_names), "spot_ids": spot_ids}, truth
