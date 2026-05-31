"""
Reference-based deconvolution solver backbones.

A *solver* maps a query (bulk genes×samples, or spatial spot pseudo-profiles)
plus a :class:`~tissueresolve.results.ReferenceSignature` to per-observation
proportions.  These are pure, dependency-light backbones used by ``solver=auto``
selection, the ensemble, and the benchmark.  They do **not** replace the
protocol-aware ``BulkPipeline`` (which is the ``weighted_nnls``-style backbone);
they make alternative backbones available so cross-validation can pick the best.
"""
from __future__ import annotations

from dataclasses import dataclass, field
from typing import Any, Optional

import numpy as np
import pandas as pd

from tissueresolve.results import ReferenceSignature

__all__ = ["SolverResult", "BaseSolver", "align_query_to_reference"]


@dataclass
class SolverResult:
    proportions: pd.DataFrame                 # obs × cell_type, rows sum to 1
    genes_used: list[str]
    diagnostics: dict[str, Any] = field(default_factory=dict)


def align_query_to_reference(query: pd.DataFrame, ref: ReferenceSignature,
                             genes: Optional[list[str]] = None):
    """Return (sub_ref, B) with B as genes×samples aligned to the reference.

    *query* is genes×samples.  *genes* optionally restricts the panel; genes not
    in the reference are dropped (never silently zero-filled)."""
    ref_genes = set(map(str, ref.gene_names))
    q = query.copy()
    q.index = q.index.map(str)
    candidate = [g for g in (genes or list(q.index)) if g in ref_genes and g in set(q.index)]
    if not candidate:
        raise ValueError("no shared genes between query and reference")
    sub = ref.subset_genes(candidate)
    B = q.loc[list(sub.gene_names)].to_numpy(float)
    return sub, B


def _normalize(x: np.ndarray) -> np.ndarray:
    s = x.sum()
    return x / s if s > 0 else np.full(len(x), 1.0 / len(x))


class BaseSolver:
    name = "base"

    def solve(self, query: pd.DataFrame, ref: ReferenceSignature) -> SolverResult:
        raise NotImplementedError

    # shared helpers -----------------------------------------------------
    @staticmethod
    def _frame(props: np.ndarray, samples, cell_types) -> pd.DataFrame:
        return pd.DataFrame(props, index=list(samples), columns=list(cell_types))

    @staticmethod
    def condition_number(R: np.ndarray) -> float:
        try:
            return float(np.linalg.cond(R))
        except Exception:
            return float("inf")
