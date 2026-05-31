"""Additional internal bulk baselines: weighted NNLS and marker-only NNLS.

These give more than one non-TissueResolve executable comparator without any
external dependency, so the benchmark always compares ≥3 non-TissueResolve
methods locally (plus any imported external results)."""
from __future__ import annotations

import numpy as np
import pandas as pd

from benchmarks.shared.base import BenchmarkMethod


def _shared(ref, bulk):
    ref_genes = set(map(str, ref.gene_names))
    genes = [g for g in bulk.index if g in ref_genes]
    if not genes:
        raise ValueError("no shared genes between bulk and reference")
    sub = ref.subset_genes(genes)
    return sub, bulk.loc[[g for g in sub.gene_names]]


def _solve(R, B, weights=None):
    from scipy.optimize import nnls
    if weights is not None:
        w = np.sqrt(np.clip(weights, 1e-9, None))[:, None]
        R = R * w
        B = B * w
    out = []
    for j in range(B.shape[1]):
        x, _ = nnls(R, B[:, j])
        s = x.sum()
        out.append(x / s if s > 0 else np.full(len(x), 1.0 / len(x)))
    return np.array(out)


class WeightedNNLSBaseline(BenchmarkMethod):
    name = "WNNLS_baseline"
    modality = "bulk"
    requires_raw_counts = False
    external = False

    def _run(self, scenario: dict) -> pd.DataFrame:
        ref, bulk = _shared(scenario["reference"], scenario["bulk"])
        R = ref.as_R_cpm().T              # genes × K
        # weight genes by specificity: max/sum across cell types (∈(0,1])
        row = R / (R.sum(axis=1, keepdims=True) + 1e-9)
        weights = row.max(axis=1)
        props = _solve(R, bulk.to_numpy(float), weights=weights)
        return pd.DataFrame(props, index=list(bulk.columns), columns=list(ref.cell_types))


class MarkerOnlyNNLSBaseline(BenchmarkMethod):
    name = "MarkerOnly_NNLS"
    modality = "bulk"
    requires_raw_counts = False
    external = False

    def __init__(self, top_n: int = 25):
        self.top_n = top_n

    def _run(self, scenario: dict) -> pd.DataFrame:
        ref, bulk = _shared(scenario["reference"], scenario["bulk"])
        R = ref.as_R_cpm()                # K × genes
        # top-N marker genes per cell type by log2FC vs the mean of the rest
        K, G = R.shape
        logR = np.log1p(R)
        marker_idx = set()
        for k in range(K):
            others = np.delete(logR, k, axis=0).mean(axis=0)
            lfc = logR[k] - others
            marker_idx.update(np.argsort(lfc)[::-1][:self.top_n].tolist())
        idx = sorted(marker_idx)
        genes = [ref.gene_names[i] for i in idx]
        Rm = R[:, idx].T                  # markers × K
        Bm = bulk.loc[genes].to_numpy(float)
        props = _solve(Rm, Bm)
        return pd.DataFrame(props, index=list(bulk.columns), columns=list(ref.cell_types))
