import os
from pathlib import Path

import pytest

from tissueresolve import cli


def test_run_help():
    with pytest.raises(SystemExit) as e:
        cli.run(["--help"])
    assert e.type is SystemExit


def test_run_dry(tmp_path):
    # create tiny reference and query files
    ref = tmp_path / "ref.h5ad"
    import anndata as ad
    import numpy as np

    adata = ad.AnnData(X=np.ones((3, 2)))
    adata.obs["cell_type"] = ["A"] * 3
    adata.write(ref)
    query = tmp_path / "counts.tsv"
    query.write_text("gene\ts1\nG1\t1\n")
    out = tmp_path / "out"
    rc = cli.run(["--reference", str(ref), "--query", str(query), "--out", str(out), "--dry-run"]) 
    assert rc == 0
    assert (out / "analysis_plan.json").exists()
