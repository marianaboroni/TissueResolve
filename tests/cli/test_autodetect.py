from pathlib import Path

import anndata as ad
import numpy as np

from tissueresolve.io import autodetect


def test_detect_h5ad_spatial(tmp_path):
    adata = ad.AnnData(X=np.ones((3, 4)))
    adata.obsm["spatial"] = np.zeros((3, 2))
    p = tmp_path / "spatial.h5ad"
    adata.write(p)
    t = autodetect.detect_input_type(p)
    assert t == "spatial_visium_h5ad"


def test_detect_h5ad_reference(tmp_path):
    adata = ad.AnnData(X=np.ones((10, 5)))
    adata.obs["cell_type"] = ["A"] * 10
    p = tmp_path / "ref.h5ad"
    adata.write(p)
    t = autodetect.detect_input_type(p)
    assert t == "single_cell_h5ad_reference"


def test_detect_table(tmp_path):
    p = tmp_path / "counts.tsv"
    p.write_text("# comment\tline\ngene1\t1\n")
    t = autodetect.detect_input_type(p)
    assert t == "bulk_counts_table"
