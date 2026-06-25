"""Internal NNLS bulk baseline — always available, no external dependencies."""
from __future__ import annotations

import numpy as np
import pandas as pd

from benchmarks.shared.base import BenchmarkMethod


class NNLSBaseline(BenchmarkMethod):
    name = "NNLS_baseline"
    modality = "bulk"
    requires_raw_counts = False
    supports_hierarchical_reference = False
    external = False

    def _run(self, scenario: dict) -> pd.DataFrame:
        from scipy.optimize import nnls

        ref = scenario["reference"]            # ReferenceSignature
        bulk = scenario["bulk"]                # genes × samples DataFrame
        genes = [g for g in bulk.index if g in set(map(str, ref.gene_names))]
        if not genes:
            raise ValueError("no shared genes between bulk and reference")
        sub = ref.subset_genes(genes)
        R = sub.as_R_cpm().T                   # genes × K
        gene_idx = {g: i for i, g in enumerate(sub.gene_names)}
        B = bulk.loc[[g for g in sub.gene_names]].to_numpy(float)  # genes × samples
        props = []
        for j in range(B.shape[1]):
            x, _ = nnls(R, B[:, j])
            s = x.sum()
            props.append(x / s if s > 0 else np.full(len(x), 1.0 / len(x)))
        return pd.DataFrame(props, index=list(bulk.columns),
                            columns=list(sub.cell_types))
