"""Independent sibling / rare-confirmation evidence scores (Stage 1.5, experimental).

Transparent, deterministic, INFERENCE-computable scores that use the `sibling` and
`rare_confirmation` panels as **independent evidence** about a fine subtype's presence —
NOT as an abundance estimator and NOT as a calibrated probability. They exist to TEST
(Strategy E) whether these panels carry information beyond the `fine_global` abundance
estimate. Nothing here changes proportions or defaults.

Scores are per (sample, subtype), in interpretable units:
  * rare_confirmation_evidence — fraction of the subtype's confirmation genes observed
    above the sample background (0..1). "Do the specific confirmation genes show up?"
  * sibling_evidence — mean log2-CPM of the subtype's sibling-discriminating markers
    minus that of its siblings' markers. ">0 = query looks like this subtype within its
    family; <0 = looks like a sibling (confounded)."

We deliberately do NOT call these probabilities.
"""
from __future__ import annotations

from typing import Mapping, Optional, Sequence

import numpy as np
import pandas as pd

__all__ = ["to_log_cpm", "rare_confirmation_evidence", "sibling_evidence", "FEATURE_STATUS"]

FEATURE_STATUS = "experimental"


def to_log_cpm(query: pd.DataFrame) -> pd.DataFrame:
    """query genes×samples counts -> genes×samples log2-CPM (str gene index)."""
    q = query.copy()
    q.index = q.index.map(str)
    col = q.sum(axis=0).replace(0, np.nan)
    return np.log2(q.div(col, axis=1).fillna(0.0) * 1e6 + 1.0)


def rare_confirmation_evidence(query: pd.DataFrame,
                               rare_confirmation: Mapping[str, Sequence[str]],
                               *, background_quantile: float = 0.5,
                               already_logcpm: bool = False) -> pd.DataFrame:
    """Per (sample, subtype): fraction of the subtype's confirmation genes expressed
    above the sample's background (``background_quantile`` of nonzero log-CPM)."""
    lc = query if already_logcpm else to_log_cpm(query)
    gidx = set(lc.index)
    samples = list(lc.columns)
    out = {}
    # per-sample background threshold from nonzero genes
    bg = {}
    for s in samples:
        v = lc[s].to_numpy(float); nz = v[v > 0]
        bg[s] = float(np.quantile(nz, background_quantile)) if nz.size else 0.0
    for sub, panel in rare_confirmation.items():
        present = [g for g in map(str, panel) if g in gidx]
        if not present:
            out[sub] = pd.Series(np.nan, index=samples); continue
        sub_lc = lc.loc[present]
        out[sub] = pd.Series(
            {s: float((sub_lc[s].to_numpy(float) > bg[s]).mean()) for s in samples})
    return pd.DataFrame(out).reindex(index=samples)


def sibling_evidence(query: pd.DataFrame,
                     sibling_genes_by_subtype: Mapping[str, Sequence[str]],
                     mapping: Mapping[str, str],
                     *, already_logcpm: bool = False) -> pd.DataFrame:
    """Per (sample, subtype): mean log2-CPM of the subtype's sibling markers minus the
    mean of its siblings' markers (within the same broad family). >0 favours the subtype."""
    lc = query if already_logcpm else to_log_cpm(query)
    gidx = set(lc.index)
    samples = list(lc.columns)
    # group subtypes by broad family
    fam = {}
    for s in sibling_genes_by_subtype:
        fam.setdefault(str(mapping.get(s, s)), []).append(s)
    out = {}
    for f, members in fam.items():
        for s in members:
            own = [g for g in map(str, sibling_genes_by_subtype.get(s, [])) if g in gidx]
            others = [g for m in members if m != s
                      for g in map(str, sibling_genes_by_subtype.get(m, [])) if g in gidx]
            if not own:
                out[s] = pd.Series(np.nan, index=samples); continue
            own_mean = lc.loc[own].mean(axis=0)
            oth_mean = lc.loc[others].mean(axis=0) if others else pd.Series(0.0, index=samples)
            out[s] = (own_mean - oth_mean).reindex(samples)
    return pd.DataFrame(out).reindex(index=samples)
