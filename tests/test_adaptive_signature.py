"""Offline synthetic validation for the separability-aware adaptive budget (Stage 1.6)."""
from __future__ import annotations

import numpy as np
import pandas as pd
import pytest

anndata = pytest.importorskip("anndata")


@pytest.fixture
def synthetic():
    """6 types × 5 donors with designed difficulty:
      Easy   — unique block g0-19 (highly separable)
      CollA  — g20-24 unique + g25-29 shared (moderate)
      CollB  — g30-34 unique + g25-29 shared (moderate, confounds CollA)
      Rare   — unique block g40-49 but only 6 cells/donor
      TwinA  — block g60-69, NO unique markers (shared with TwinB)
      TwinB  — block g60-69, identical to TwinA -> structurally unresolvable
    """
    rng = np.random.default_rng(0)
    G = 80
    genes = [f"g{i}" for i in range(G)]
    blocks = {"Easy": range(0, 20), "CollA": list(range(20, 25)) + list(range(25, 30)),
              "CollB": list(range(30, 35)) + list(range(25, 30)),
              "Rare": range(40, 50), "TwinA": range(60, 70), "TwinB": range(60, 70)}
    n_per = {"Easy": 40, "CollA": 40, "CollB": 40, "Rare": 6, "TwinA": 40, "TwinB": 40}
    X, ct, dn = [], [], []
    for d in [f"d{i}" for i in range(5)]:
        for t, blk in blocks.items():
            for _ in range(n_per[t]):
                v = rng.poisson(2.0, size=G).astype(float)
                v[list(blk)] += rng.poisson(40.0, size=len(list(blk)))
                X.append(v); ct.append(t); dn.append(d)
    ad = anndata.AnnData(np.vstack(X).astype("float32"))
    ad.var_names = genes
    ad.obs["cell_type"] = pd.Categorical(ct)
    ad.obs["donor_id"] = pd.Categorical(dn)
    return ad


MAP = {t: "fam" for t in ("Easy", "CollA", "CollB", "Rare", "TwinA", "TwinB")}


def _alloc(**kw):
    from tissueresolve.reference.adaptive_signature import SeparabilityAwareBudgetAllocator
    return SeparabilityAwareBudgetAllocator(
        "cell_type", "donor_id", MAP, min_cells=5, sizes=(5, 10, 15, 20, 30),
        max_genes=30, unresolvable_rmse=0.25, **kw)


def test_all_types_represented(synthetic):
    a = _alloc().allocate(synthetic)
    assert set(a.per_type) == set(MAP)
    assert a.genes and len(a.genes) == len(set(a.genes))    # deduped global signature


def test_easy_is_compact_and_resolvable(synthetic):
    a = _alloc().allocate(synthetic)
    easy = a.per_type["Easy"]
    assert easy["status"] in ("EASY_COMPACT", "MODERATE")
    assert easy["best_rmse"] is not None and easy["best_rmse"] < 0.2


def test_twins_are_unresolvable(synthetic):
    """Types with no unique markers (identical twin) must be flagged, not expanded forever."""
    a = _alloc().allocate(synthetic)
    assert a.per_type["TwinA"]["status"] == "UNRESOLVABLE"
    assert a.per_type["TwinB"]["status"] == "UNRESOLVABLE"


def test_easy_easier_than_collinear(synthetic):
    a = _alloc().allocate(synthetic)
    assert a.per_type["Easy"]["best_rmse"] <= a.per_type["CollA"]["best_rmse"] + 1e-9


def test_adaptive_differentiates_by_difficulty(synthetic):
    """The allocator must differentiate populations by difficulty. In this clean toy the
    resolvable types all resolve compactly, so the coherent adaptive signal is STATUS:
    easily-separable types → EASY_COMPACT; the identical twins → UNRESOLVABLE.
    (Gene-count adaptivity is exercised on the real-tissue benchmark.)"""
    a = _alloc().allocate(synthetic)
    statuses = {t: v["status"] for t, v in a.per_type.items()}
    assert len(set(statuses.values())) > 1, f"no differentiation: {statuses}"
    assert "UNRESOLVABLE" in statuses.values() and "EASY_COMPACT" in statuses.values()


def test_deterministic(synthetic):
    a = _alloc().allocate(synthetic); b = _alloc().allocate(synthetic)
    assert a.genes == b.genes
    assert {k: v["gene_count"] for k, v in a.per_type.items()} == \
           {k: v["gene_count"] for k, v in b.per_type.items()}


def test_no_donor_fallback(synthetic):
    from tissueresolve.reference.adaptive_signature import SeparabilityAwareBudgetAllocator
    ad = synthetic.copy(); del ad.obs["donor_id"]
    import warnings
    with warnings.catch_warnings(record=True) as w:
        warnings.simplefilter("always")
        a = SeparabilityAwareBudgetAllocator("cell_type", "donor_id", MAP,
                                             min_cells=5).allocate(ad)
    assert all(v["status"] == "INSUFFICIENT_REFERENCE" for v in a.per_type.values())
    assert any("donor" in str(x.message).lower() for x in w)


def test_save_outputs(tmp_path, synthetic):
    a = _alloc().allocate(synthetic)
    a.save(tmp_path / "adapt")
    for f in ("adaptive_signature.tsv", "adaptive_gene_budget_by_celltype.tsv",
              "celltype_budget_curve.tsv", "adaptive_signature_manifest.json"):
        assert (tmp_path / "adapt" / f).exists(), f"missing {f}"
