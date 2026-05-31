"""
Gene-masking cross-validation.

No external ground truth is needed: hide a fraction of shared genes, deconvolve
on the rest, reconstruct the held-out genes from estimated proportions ×
reference signatures, and score reconstruction against the observed held-out
expression.  Used to select the solver backbone, marker count, and spatial
smoothing — choosing parameters that *generalise*, not that overfit one panel.
"""
from __future__ import annotations

from dataclasses import dataclass, field
from typing import Any, Callable, Optional

import numpy as np
import pandas as pd

from tissueresolve.results import ReferenceSignature

__all__ = [
    "split_genes_for_masking", "reconstruct_masked_genes",
    "score_masked_gene_reconstruction", "run_gene_masking_cv_bulk",
    "compare_solvers_by_gene_masking", "GeneMaskingResult",
]


@dataclass
class GeneMaskingResult:
    score: float                       # mean masked-gene Pearson across folds
    rmse: float
    per_fold: list[dict] = field(default_factory=list)
    diagnostics: dict[str, Any] = field(default_factory=dict)


def split_genes_for_masking(genes: list[str], *, mask_fraction: float = 0.15,
                            n_splits: int = 3, seed: int = 0):
    """Deterministic gene folds: list of (train_genes, masked_genes)."""
    genes = list(map(str, genes))
    rng = np.random.default_rng(seed)
    n_mask = max(1, int(round(len(genes) * mask_fraction)))
    folds = []
    for s in range(n_splits):
        # seed varies by fold but is reproducible (no Math.random equivalent)
        r = np.random.default_rng(seed + s + 1)
        masked = set(r.choice(len(genes), size=n_mask, replace=False).tolist())
        folds.append(([g for i, g in enumerate(genes) if i not in masked],
                      [g for i, g in enumerate(genes) if i in masked]))
    return folds


def reconstruct_masked_genes(proportions: pd.DataFrame, ref: ReferenceSignature,
                             masked_genes: list[str]) -> pd.DataFrame:
    """Reconstruct samples×masked_genes = proportions · reference[masked]."""
    sub = ref.subset_genes(masked_genes)
    # align cell-type columns
    cols = [c for c in proportions.columns if c in set(sub.cell_types)]
    P = proportions[cols].to_numpy(float)
    idx = [list(sub.cell_types).index(c) for c in cols]
    R = sub.as_R_cpm()[idx]            # K' × masked_genes
    recon = P @ R                      # samples × masked_genes
    return pd.DataFrame(recon, index=proportions.index, columns=list(sub.gene_names))


def score_masked_gene_reconstruction(observed: pd.DataFrame,
                                     reconstructed: pd.DataFrame) -> dict:
    """Pearson/Spearman/RMSE between observed and reconstructed masked genes.

    Both are scaled to per-sample relative expression before comparison so
    library-size differences do not dominate."""
    from scipy.stats import pearsonr, spearmanr
    genes = [g for g in observed.columns if g in set(reconstructed.columns)]
    obs = observed[genes].to_numpy(float)
    rec = reconstructed[genes].to_numpy(float)
    obs = obs / (obs.sum(axis=1, keepdims=True) + 1e-9)
    rec = rec / (rec.sum(axis=1, keepdims=True) + 1e-9)
    ov, rv = obs.ravel(), rec.ravel()
    if ov.std() < 1e-12 or rv.std() < 1e-12:
        pear = spear = float("nan")
    else:
        pear = float(pearsonr(ov, rv)[0])
        spear = float(spearmanr(ov, rv)[0])
    rmse = float(np.sqrt(np.mean((ov - rv) ** 2)))
    return {"masked_gene_pearson": pear, "masked_gene_spearman": spear,
            "masked_gene_rmse": rmse}


def run_gene_masking_cv_bulk(bulk: pd.DataFrame, ref: ReferenceSignature,
                             solve_fn: Callable, *, mask_fraction: float = 0.15,
                             n_splits: int = 3, seed: int = 0) -> GeneMaskingResult:
    """Run gene-masking CV for a bulk *solve_fn*.

    *solve_fn(bulk_train, ref, train_genes) -> proportions DataFrame*.
    Reconstructs masked genes from those proportions and scores them.
    """
    bulk = bulk.copy()
    bulk.index = bulk.index.map(str)
    ref_genes = set(map(str, ref.gene_names))
    shared = [g for g in bulk.index if g in ref_genes]
    folds = split_genes_for_masking(shared, mask_fraction=mask_fraction,
                                    n_splits=n_splits, seed=seed)
    per_fold, pears, rmses = [], [], []
    for train_genes, masked_genes in folds:
        props = solve_fn(bulk.loc[train_genes], ref, train_genes)
        recon = reconstruct_masked_genes(props, ref, masked_genes)
        observed = bulk.loc[masked_genes].T  # samples × masked_genes
        sc = score_masked_gene_reconstruction(observed, recon)
        per_fold.append(sc)
        if not np.isnan(sc["masked_gene_pearson"]):
            pears.append(sc["masked_gene_pearson"])
        rmses.append(sc["masked_gene_rmse"])
    return GeneMaskingResult(
        score=float(np.mean(pears)) if pears else float("nan"),
        rmse=float(np.mean(rmses)) if rmses else float("nan"),
        per_fold=per_fold,
        diagnostics={"n_splits": n_splits, "mask_fraction": mask_fraction,
                     "n_shared_genes": len(shared)})


def compare_solvers_by_gene_masking(bulk: pd.DataFrame, ref: ReferenceSignature,
                                    solvers: dict, *, mask_fraction: float = 0.15,
                                    n_splits: int = 3, seed: int = 0) -> pd.DataFrame:
    """Score each named solver by gene-masking CV.  Returns a sorted DataFrame.

    *solvers* maps name → ``solve_fn(bulk_train, ref, train_genes)``.
    """
    rows = []
    for name, fn in solvers.items():
        try:
            res = run_gene_masking_cv_bulk(bulk, ref, fn, mask_fraction=mask_fraction,
                                           n_splits=n_splits, seed=seed)
            rows.append({"solver": name, "masked_gene_pearson": res.score,
                         "masked_gene_rmse": res.rmse})
        except Exception as exc:  # noqa: BLE001
            rows.append({"solver": name, "masked_gene_pearson": float("nan"),
                         "masked_gene_rmse": float("nan"), "error": str(exc)})
    df = pd.DataFrame(rows).set_index("solver")
    return df.sort_values("masked_gene_pearson", ascending=False)
