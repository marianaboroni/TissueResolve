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


def test_run_dry_with_resolution_mode(tmp_path):
    ref = tmp_path / "ref.h5ad"
    import anndata as ad
    import numpy as np

    adata = ad.AnnData(X=np.ones((3, 2)))
    adata.obs["cell_type"] = ["A"] * 3
    adata.write(ref)
    query = tmp_path / "counts.tsv"
    query.write_text("gene\ts1\nG1\t1\n")
    out = tmp_path / "out"
    rc = cli.run([
        "--reference", str(ref),
        "--query", str(query),
        "--out", str(out),
        "--resolution-mode", "hierarchical",
        "--dry-run",
    ])
    assert rc == 0
    assert (out / "analysis_plan.json").exists()


def test_cli_exposes_hierarchical_flags():
    """The hierarchical broad/fine flags are exposed on `tissueresolve run`."""
    from click.testing import CliRunner

    result = CliRunner().invoke(cli.cli, ["run", "--help"])
    assert result.exit_code == 0
    for flag in ("--broad-cell-type-col", "--fine-cell-type-col",
                 "--cell-type-hierarchy", "--allow-unresolved",
                 "--resolution-mode"):
        assert flag in result.output


def test_run_dry_hierarchical_records_plan(tmp_path):
    """Dry-run in hierarchical mode records broad/fine columns in the plan."""
    import json

    import anndata as ad
    import numpy as np

    ref = tmp_path / "ref.h5ad"
    adata = ad.AnnData(X=np.ones((4, 2)))
    adata.obs["broad_cell_type"] = ["T/NK", "T/NK", "Myeloid", "Myeloid"]
    adata.obs["sub_cell_type"] = ["CD4 T", "CD8 T", "macro", "mono"]
    adata.write(ref)
    query = tmp_path / "counts.tsv"
    query.write_text("gene\ts1\nG1\t1\n")
    out = tmp_path / "out"
    rc = cli.run([
        "--reference", str(ref), "--query", str(query), "--out", str(out),
        "--mode", "bulk", "--resolution-mode", "hierarchical",
        "--broad-cell-type-col", "broad_cell_type",
        "--fine-cell-type-col", "sub_cell_type", "--dry-run",
    ])
    assert rc == 0
    plan = json.loads((out / "analysis_plan.json").read_text())
    assert plan["resolution_mode"] == "hierarchical"
    assert plan["hierarchical"]["broad_cell_type_col"] == "broad_cell_type"
    assert plan["hierarchical"]["fine_cell_type_col"] == "sub_cell_type"
