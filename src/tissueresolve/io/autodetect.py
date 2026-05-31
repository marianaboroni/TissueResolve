from __future__ import annotations

from pathlib import Path
from typing import Literal

import anndata as ad
import pandas as pd

InputType = Literal[
    "single_cell_h5ad_reference",
    "bulk_counts_table",
    "spatial_visium_h5ad",
    "spatial_visium_folder",
    "unknown",
]


def inspect_anndata(path: Path) -> dict:
    adata = ad.read_h5ad(path)
    info = {
        "n_obs": adata.n_obs,
        "n_vars": adata.n_vars,
        "has_spatial": "spatial" in adata.obsm or any(k.lower().startswith("spatial") for k in adata.obsm.keys()),
        "has_celltype": any(c.lower().startswith("cell") or "cell_type" in c.lower() for c in adata.obs.columns),
        "uns_keys": list(adata.uns.keys()),
    }
    return info


def inspect_visium_folder(path: Path) -> bool:
    # Look for typical Visium outputs
    if not path.is_dir():
        return False
    for fname in ("filtered_feature_bc_matrix.h5", "filtered_feature_bc_matrix.h5.h5", "filtered_feature_bc_matrix", "matrix.mtx"):
        if (path / fname).exists():
            return True
    if (path / "spatial").exists():
        return True
    return False


def inspect_table(path: Path) -> dict:
    # Read a small slice and infer orientation
    try:
        df = pd.read_csv(path, sep="\t", comment="#", nrows=10, index_col=0)
    except Exception:
        try:
            df = pd.read_csv(path, sep=",", comment="#", nrows=10, index_col=0)
        except Exception:
            return {"ok": False}
    return {"ok": True, "nrows": df.shape[0], "ncols": df.shape[1]}


def detect_input_type(path: str | Path) -> InputType:
    p = Path(path)
    if not p.exists():
        return "unknown"
    if p.suffix == ".h5ad":
        try:
            info = inspect_anndata(p)
            if info.get("has_spatial"):
                return "spatial_visium_h5ad"
            if info.get("has_celltype"):
                return "single_cell_h5ad_reference"
            return "unknown"
        except Exception:
            return "unknown"
    if p.is_dir():
        if inspect_visium_folder(p):
            return "spatial_visium_folder"
        return "unknown"
    # assume table
    tbl = inspect_table(p)
    if tbl.get("ok"):
        return "bulk_counts_table"
    return "unknown"


def summarize_input_detection(path: str | Path) -> dict:
    t = detect_input_type(path)
    return {"path": str(path), "detected_type": t}
