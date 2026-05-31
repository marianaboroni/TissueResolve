"""Shared input export for spatial external methods (manual/web execution)."""
from __future__ import annotations

import numpy as np
import pandas as pd
import scipy.sparse as sp

from benchmarks.shared.io import write_tsv, OUTPUTS_DIR


def export_spatial_inputs(scenario: dict, method: str) -> list[str]:
    out = OUTPUTS_DIR / "spatial" / "exports" / method
    files = []
    Y = scenario.get("Y")
    gene_names = [str(g) for g in scenario.get("gene_names", [])]
    spot_ids = scenario.get("spot_ids") or (
        [f"spot{i}" for i in range(Y.shape[0])] if Y is not None else [])
    if Y is not None:
        dense = Y.toarray() if sp.issparse(Y) else np.asarray(Y)
        files.append(str(write_tsv(
            pd.DataFrame(dense, index=spot_ids, columns=gene_names),
            out / "spatial_counts.tsv")))
    ref = scenario.get("reference")
    if ref is not None:
        sig = pd.DataFrame(ref.as_R_cpm().T, index=list(ref.gene_names),
                           columns=list(ref.cell_types))
        files.append(str(write_tsv(sig, out / "reference_signature_cpm.tsv")))
    if "array_row" in scenario and "array_col" in scenario:
        coords = pd.DataFrame({"array_row": scenario["array_row"],
                               "array_col": scenario["array_col"]}, index=spot_ids)
        files.append(str(write_tsv(coords, out / "spot_coordinates.tsv")))
    return files
