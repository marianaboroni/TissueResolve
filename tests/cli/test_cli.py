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


def _write_ref_with_broad_fine(tmp_path):
    import anndata as ad
    import numpy as np
    ref = tmp_path / "ref.h5ad"
    adata = ad.AnnData(X=np.ones((4, 3)))
    adata.obs["broad_cell_type"] = ["T/NK", "T/NK", "Myeloid", "Myeloid"]
    adata.obs["sub_cell_type"] = ["CD4 T", "CD8 T", "macro", "mono"]
    adata.write(ref)
    return ref


def _write_ref_fine_only(tmp_path):
    import anndata as ad
    import numpy as np
    ref = tmp_path / "ref_fine.h5ad"
    adata = ad.AnnData(X=np.ones((3, 3)))
    adata.obs["cell_type"] = ["A", "B", "C"]
    adata.write(ref)
    return ref


def _bulk_query(tmp_path):
    q = tmp_path / "counts.tsv"
    q.write_text("gene\ts1\nG1\t1\n")
    return q


def test_auto_selects_hierarchical_when_broad_fine_present(tmp_path):
    import json
    ref = _write_ref_with_broad_fine(tmp_path)
    out = tmp_path / "out"
    rc = cli.run(["--reference", str(ref), "--query", str(_bulk_query(tmp_path)),
                  "--out", str(out), "--mode", "bulk",
                  "--resolution-mode", "auto", "--dry-run"])
    assert rc == 0
    plan = json.loads((out / "analysis_plan.json").read_text())
    assert plan["requested_resolution_mode"] == "auto"
    assert plan["resolution_mode"] == "hierarchical"
    assert "reason" in plan["resolution_mode_reason"].lower() or plan["resolution_mode_reason"]


def test_auto_falls_back_to_flat_without_labels(tmp_path):
    import json
    ref = _write_ref_fine_only(tmp_path)
    out = tmp_path / "out2"
    rc = cli.run(["--reference", str(ref), "--query", str(_bulk_query(tmp_path)),
                  "--out", str(out), "--mode", "bulk",
                  "--resolution-mode", "auto", "--preset", "standard", "--dry-run"])
    assert rc == 0
    plan = json.loads((out / "analysis_plan.json").read_text())
    assert plan["resolution_mode"] == "none"  # flat fallback


def test_flat_is_explicit(tmp_path):
    import json
    ref = _write_ref_with_broad_fine(tmp_path)
    out = tmp_path / "out3"
    rc = cli.run(["--reference", str(ref), "--query", str(_bulk_query(tmp_path)),
                  "--out", str(out), "--mode", "bulk",
                  "--resolution-mode", "flat", "--dry-run"])
    assert rc == 0
    plan = json.loads((out / "analysis_plan.json").read_text())
    # flat is honoured even though broad/fine labels exist
    assert plan["resolution_mode"] == "none"
    assert plan["requested_resolution_mode"] == "flat"


def test_auto_publication_without_labels_stops(tmp_path):
    ref = _write_ref_fine_only(tmp_path)
    out = tmp_path / "out4"
    with pytest.raises(Exception):
        cli.run(["--reference", str(ref), "--query", str(_bulk_query(tmp_path)),
                 "--out", str(out), "--mode", "bulk",
                 "--resolution-mode", "auto", "--preset", "publication", "--dry-run"])


def test_mapping_file_makes_auto_hierarchical(tmp_path):
    import json
    ref = _write_ref_fine_only(tmp_path)
    mapping = tmp_path / "map.tsv"
    mapping.write_text("fine_cell_type\tbroad_cell_type\nA\tFam1\nB\tFam1\nC\tFam2\n")
    out = tmp_path / "out5"
    rc = cli.run(["--reference", str(ref), "--query", str(_bulk_query(tmp_path)),
                  "--out", str(out), "--mode", "bulk", "--resolution-mode", "auto",
                  "--cell-type-hierarchy", str(mapping), "--dry-run"])
    assert rc == 0
    plan = json.loads((out / "analysis_plan.json").read_text())
    assert plan["resolution_mode"] == "hierarchical"
