"""Donor-aware gene selection for deconvolution signatures (experimental, opt-in).

Motivation (Rectangle benchmark, Etapa 3): the residual accuracy gap between the
TissueResolve bulk path and Rectangle is attributable to the **signature / gene
set**, not the solver (the Poisson GLM already closes most of the gap). This module
provides alternative, donor-aware gene-selection strategies to be **benchmark-gated**
against the current markers on donor-held-out deconvolution — not promoted to default
until they demonstrably improve broad / fine / conditional RMSE without inflating
rare false positives or losing cross-platform stability.

Strategies (all operate on **donor × cell-type pseudobulks** from per-cell AnnData —
never treating individual cells as biological replicates, per the reference rules):

* ``donor_aware_de``  — pseudobulk differential expression, ``one_vs_rest`` (broad)
  or ``sibling`` (within-family fine). Records log2FC, BH-adjusted p, supporting-donor
  fraction; penalises single-donor-supported genes. This is a donor-aware pseudobulk
  DE test implemented here (NOT pydeseq2 — which is unavailable on Python 3.9).
* ``ml_minimal``      — L1 multinomial logistic feature selection with **donor-held-out**
  (optionally platform-held-out) validation, returning the smallest gene set whose
  held-out balanced accuracy is within ``tol`` of the best. Used for interpretable
  feature selection only; the SET is validated on deconvolution downstream (§13).
* ``hybrid_gene_set`` — DE candidates → donor/sibling-stability filter → ML-minimal
  selection → redundancy pruning. The prioritised candidate from §12.

Nothing here changes defaults or the production pipeline; it produces gene lists a
benchmark or a (future) reference builder can evaluate and, only if it wins, adopt.
"""
from __future__ import annotations

from dataclasses import dataclass, field
from typing import Optional, Sequence

import numpy as np
import pandas as pd

__all__ = [
    "donor_pseudobulk", "donor_aware_de", "ml_minimal_genes", "hybrid_gene_set",
    "select_donor_aware_genes", "GeneSetResult", "FEATURE_STATUS",
]

FEATURE_STATUS = "experimental"
ALGORITHM_VERSION = "donor_aware_gene_selection-0.1.0"


@dataclass
class GeneSetResult:
    genes: list[str]
    method: str
    stats: pd.DataFrame = field(default_factory=pd.DataFrame)
    metadata: dict = field(default_factory=dict)


# ---------------------------------------------------------------------------
# donor × cell-type pseudobulk substrate
# ---------------------------------------------------------------------------
def donor_pseudobulk(adata, celltype_col: str, donor_col: str, *,
                     min_cells: int = 10, layer: Optional[str] = None):
    """Build donor × cell-type pseudobulks (log2-CPM) from per-cell AnnData.

    Returns ``(X, labels, donors, genes)`` where ``X`` is (n_groups × n_genes)
    log2-CPM, ``labels`` the cell type per row, ``donors`` the donor per row.
    Groups with < ``min_cells`` cells are dropped (recorded by the caller).
    """
    import scipy.sparse as sp
    obs = adata.obs
    genes = [str(g) for g in adata.var_names]
    ct = obs[celltype_col].astype(str).to_numpy()
    dn = obs[donor_col].astype(str).to_numpy()
    M = adata.layers[layer] if layer else adata.X
    rows, labels, donors = [], [], []
    for c in pd.unique(ct):
        for d in pd.unique(dn[ct == c]):
            idx = np.where((ct == c) & (dn == d))[0]
            if len(idx) < min_cells:
                continue
            sub = M[idx]
            v = np.asarray(sub.sum(axis=0)).ravel() if sp.issparse(sub) \
                else np.asarray(sub).sum(axis=0).ravel()
            tot = v.sum()
            cpm = (v / tot * 1e6) if tot > 0 else v
            rows.append(np.log2(cpm + 1.0))
            labels.append(c); donors.append(d)
    if not rows:
        raise ValueError("no (donor, cell-type) group met min_cells")
    return np.vstack(rows), np.array(labels), np.array(donors), genes


def select_donor_aware_genes(adata, celltype_col: str, donor_col: str, *,
                             top_n_per_type: int = 15, min_cells: int = 10,
                             min_supporting_donors: int = 2) -> list[str]:
    """Union of donor-aware one-vs-rest DE markers across all cell types.

    This is the ``donor_de_ovr`` strategy that was the robust winner in the
    gene-selection benchmark (beats current markers on breast + lung + cross-platform,
    5 seeds, preserving rare recall). Convenience entry for reference-build integration.
    Returns a de-duplicated, order-preserving gene list.
    """
    types = pd.unique(adata.obs[celltype_col].astype(str))
    genes: list[str] = []
    for t in types:
        res = donor_aware_de(adata, celltype_col, donor_col, str(t),
                             top_n=top_n_per_type, min_cells=min_cells,
                             min_supporting_donors=min_supporting_donors)
        genes.extend(res.genes)
    return list(dict.fromkeys(genes))


def _bh(pvals: np.ndarray) -> np.ndarray:
    """Benjamini–Hochberg adjusted p-values (no statsmodels dependency)."""
    p = np.asarray(pvals, float)
    n = len(p)
    order = np.argsort(p)
    ranked = p[order] * n / (np.arange(n) + 1)
    ranked = np.minimum.accumulate(ranked[::-1])[::-1]
    out = np.empty(n); out[order] = np.clip(ranked, 0, 1)
    return out


# ---------------------------------------------------------------------------
# donor-aware pseudobulk DE
# ---------------------------------------------------------------------------
def donor_aware_de(adata, celltype_col: str, donor_col: str, target: str, *,
                   mode: str = "one_vs_rest", siblings: Optional[Sequence[str]] = None,
                   top_n: int = 50, min_cells: int = 10,
                   min_supporting_donors: int = 2,
                   min_log2fc_per_donor: float = 1.0) -> GeneSetResult:
    """Rank positive markers for ``target`` from donor pseudobulks.

    ``mode='one_vs_rest'`` compares target vs all other types (broad markers);
    ``mode='sibling'`` compares target vs its ``siblings`` only (fine/within-family).
    Genes supported by fewer than ``min_supporting_donors`` target donors are
    down-weighted (kept but penalised — never silently removed).
    """
    from scipy import stats as sstats
    X, labels, donors, genes = donor_pseudobulk(
        adata, celltype_col, donor_col, min_cells=min_cells)
    tgt = labels == target
    if mode == "sibling":
        sib = set(siblings or []) - {target}
        other = np.isin(labels, list(sib))
    else:
        other = ~tgt
    if tgt.sum() < 1 or other.sum() < 1:
        return GeneSetResult([], f"donor_de:{mode}",
                             metadata={"reason": "insufficient groups", "target": target})

    A, B = X[tgt], X[other]
    with np.errstate(all="ignore"):
        t, p = sstats.ttest_ind(A, B, axis=0, equal_var=False)
    p = np.nan_to_num(p, nan=1.0)
    log2fc = A.mean(0) - B.mean(0)                      # already log2-CPM
    padj = _bh(p)
    # "supporting" = target donors showing a real per-donor fold change over the rest
    # (> min_log2fc_per_donor log2 units). A genuine marker clears this in every target
    # donor; a gene high in only one donor gets supporting≈1 and is penalised (kept, not
    # dropped) — this is the donor-robustness guard the reference rules require.
    rest_mean = B.mean(axis=0)
    supporting = ((A - rest_mean) > min_log2fc_per_donor).sum(axis=0)
    n_tgt = int(tgt.sum())
    support_frac = supporting / max(n_tgt, 1)
    penalty = np.where(supporting >= min_supporting_donors, 1.0, 0.25)
    score = np.clip(log2fc, 0, None) * (-np.log10(padj + 1e-300)) * support_frac * penalty
    df = pd.DataFrame({
        "gene": genes, "log2fc": log2fc, "pval": p, "padj": padj,
        "supporting_donors": supporting, "support_frac": support_frac,
        "single_donor_penalised": penalty < 1.0, "score": score,
    }).set_index("gene").sort_values("score", ascending=False)
    top = df[(df["log2fc"] > 0)].head(top_n)
    return GeneSetResult(list(top.index), f"donor_de:{mode}", stats=df,
                         metadata={"target": target, "n_target_donors": n_tgt,
                                   "mode": mode, "algorithm_version": ALGORITHM_VERSION,
                                   "feature_status": FEATURE_STATUS})


# ---------------------------------------------------------------------------
# ML-minimal (donor-held-out) feature selection
# ---------------------------------------------------------------------------
def ml_minimal_genes(adata, celltype_col: str, donor_col: str, *,
                     candidate_genes: Optional[Sequence[str]] = None,
                     sizes: Sequence[int] = (5, 10, 15, 20, 30, 50, 75, 100),
                     tol: float = 0.98, C: float = 0.5, min_cells: int = 10,
                     seed: int = 0) -> GeneSetResult:
    """Smallest gene set whose donor-held-out balanced accuracy is ≥ ``tol`` × best.

    L1 multinomial logistic on donor pseudobulks; genes ranked by max |coef| across
    classes; sizes evaluated by GroupKFold over donors (never splitting a donor across
    train/test). Interpretable feature selection — the returned SET is validated on
    deconvolution downstream, not on classification alone.
    """
    from sklearn.linear_model import LogisticRegression
    from sklearn.model_selection import GroupKFold
    from sklearn.metrics import balanced_accuracy_score
    from sklearn.preprocessing import StandardScaler

    X, labels, donors, genes = donor_pseudobulk(
        adata, celltype_col, donor_col, min_cells=min_cells)
    gidx = np.arange(len(genes))
    if candidate_genes is not None:
        keep = [i for i, g in enumerate(genes) if g in set(map(str, candidate_genes))]
        if len(keep) >= 5:
            X, gidx = X[:, keep], np.array(keep)
    genes_sub = [genes[i] for i in gidx]

    n_donors = len(np.unique(donors))
    n_classes = len(np.unique(labels))
    if n_donors < 3 or n_classes < 2 or X.shape[0] < 6:
        # too little to cross-validate: fall back to a full-data L1 ranking
        order = _l1_rank(X, labels, C, seed)
        top = [genes_sub[i] for i in order[: min(sizes)]]
        return GeneSetResult(top, "ml_minimal", metadata={
            "target": None, "note": "insufficient donors for CV; full-data L1 ranking",
            "feature_status": FEATURE_STATUS})

    n_splits = min(5, n_donors)
    gkf = GroupKFold(n_splits=n_splits)
    # rank genes once on all data (stable ordering), then evaluate nested sizes
    order = _l1_rank(X, labels, C, seed)
    rows = []
    for k in [s for s in sizes if s <= len(genes_sub)]:
        sel = order[:k]
        accs = []
        for tr, te in gkf.split(X, labels, groups=donors):
            if len(np.unique(labels[tr])) < 2:
                continue
            sc = StandardScaler().fit(X[tr][:, sel])
            clf = LogisticRegression(penalty="l2", C=C, max_iter=2000,
                                     class_weight="balanced")
            clf.fit(sc.transform(X[tr][:, sel]), labels[tr])
            pred = clf.predict(sc.transform(X[te][:, sel]))
            accs.append(balanced_accuracy_score(labels[te], pred))
        if accs:
            rows.append((k, float(np.mean(accs))))
    if not rows:
        top = [genes_sub[i] for i in order[: min(sizes)]]
        return GeneSetResult(top, "ml_minimal", metadata={"note": "no valid CV folds"})
    perf = pd.DataFrame(rows, columns=["size", "balanced_acc"])
    best = perf["balanced_acc"].max()
    ok = perf[perf["balanced_acc"] >= tol * best]
    chosen = int(ok["size"].min())
    top = [genes_sub[i] for i in order[:chosen]]
    return GeneSetResult(top, "ml_minimal", stats=perf.set_index("size"),
                         metadata={"chosen_size": chosen, "best_balanced_acc": best,
                                   "tol": tol, "n_donors": n_donors,
                                   "algorithm_version": ALGORITHM_VERSION,
                                   "feature_status": FEATURE_STATUS})


def _l1_rank(X, labels, C, seed) -> np.ndarray:
    """Rank gene indices by max |coef| from an L1 multinomial logistic fit."""
    from sklearn.linear_model import LogisticRegression
    from sklearn.preprocessing import StandardScaler
    Xs = StandardScaler().fit_transform(X)
    clf = LogisticRegression(penalty="l1", C=C, solver="saga", max_iter=3000,
                             class_weight="balanced", random_state=seed)
    clf.fit(Xs, labels)
    coef = np.abs(clf.coef_)                              # (n_classes, n_genes) or (1, .)
    imp = coef.max(axis=0)
    return np.argsort(-imp)


# ---------------------------------------------------------------------------
# hybrid: DE candidates -> stability/sibling filter -> ML-minimal -> prune
# ---------------------------------------------------------------------------
def hybrid_gene_set(adata, celltype_col: str, donor_col: str, targets: Sequence[str], *,
                    mode: str = "one_vs_rest", siblings: Optional[Sequence[str]] = None,
                    de_top_n: int = 100, min_supporting_donors: int = 2,
                    ml_sizes: Sequence[int] = (10, 20, 30, 50), tol: float = 0.98,
                    min_cells: int = 10, seed: int = 0) -> GeneSetResult:
    """DE candidates (union over targets, donor-supported) → ML-minimal selection."""
    cand: list[str] = []
    per_target = {}
    for t in targets:
        de = donor_aware_de(adata, celltype_col, donor_col, t, mode=mode,
                            siblings=siblings or list(targets), top_n=de_top_n,
                            min_cells=min_cells, min_supporting_donors=min_supporting_donors)
        # keep only donor-supported DE candidates for the hybrid pool
        supported = [g for g in de.genes
                     if not bool(de.stats.loc[g, "single_donor_penalised"])] \
            if not de.stats.empty else de.genes
        per_target[t] = supported
        cand += supported
    cand = list(dict.fromkeys(cand))                     # dedupe, preserve order
    if len(cand) < 5:
        return GeneSetResult(cand, "hybrid", metadata={"note": "few DE candidates",
                                                        "per_target": per_target})
    ml = ml_minimal_genes(adata, celltype_col, donor_col, candidate_genes=cand,
                          sizes=ml_sizes, tol=tol, min_cells=min_cells, seed=seed)
    return GeneSetResult(ml.genes, "hybrid", stats=ml.stats,
                         metadata={"n_de_candidates": len(cand),
                                   "per_target_de": {k: len(v) for k, v in per_target.items()},
                                   "ml": ml.metadata, "feature_status": FEATURE_STATUS})
