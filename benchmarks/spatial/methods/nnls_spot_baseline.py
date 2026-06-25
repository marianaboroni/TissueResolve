"""Internal spot-level NNLS spatial baseline — always available."""
from __future__ import annotations

import numpy as np
import pandas as pd

from benchmarks.shared.base import BenchmarkMethod


class NNLSSpotBaseline(BenchmarkMethod):
    name = "NNLS_spot_baseline"
    modality = "spatial"
    requires_raw_counts = True
    supports_hierarchical_reference = False
    external = False

    def _run(self, scenario: dict) -> pd.DataFrame:
        from scipy.optimize import nnls
        import scipy.sparse as sp

        ref = scenario["reference"]
        Y = scenario["Y"]
        gene_names = [str(g) for g in scenario["gene_names"]]
        spot_ids = scenario.get("spot_ids") or [f"spot{i}" for i in range(Y.shape[0])]

        ref_genes = set(map(str, ref.gene_names))
        shared = [g for g in gene_names if g in ref_genes]
        if not shared:
            raise ValueError("no shared genes between spatial and reference")
        col_idx = [gene_names.index(g) for g in shared]
        sub = ref.subset_genes(shared)
        R = sub.as_R_cpm().T  # genes × K
        Yd = (Y.tocsr()[:, col_idx].toarray() if sp.issparse(Y)
              else np.asarray(Y)[:, col_idx]).astype(float)
        props = []
        for i in range(Yd.shape[0]):
            x, _ = nnls(R, Yd[i])
            s = x.sum()
            props.append(x / s if s > 0 else np.full(len(x), 1.0 / len(x)))
        return pd.DataFrame(props, index=list(spot_ids), columns=list(sub.cell_types))
