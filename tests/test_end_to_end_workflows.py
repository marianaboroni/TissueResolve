"""
Offline end-to-end workflow regression tests.

These run the *real* bulk and spatial pipelines through the CLI on tiny
deterministic synthetic fixtures (no downloads, no network) and assert the
**actual** output contract of ``tissueresolve run`` — predictions + QC +
metadata JSON.  They also confirm report generation works independently from a
report-shaped results directory.

``run`` writes ``deconv/``, ``qc/``, ``analysis_plan.json``,
``run_metadata.json``, plus the report bundle ``methods.txt``, ``warnings.json``
and ``report.html`` (from the in-memory result; no standalone figure files are
rendered by ``run`` — those come from the report layer / harness).  See
``docs/END_TO_END_WORKFLOW_AUDIT.md``.
"""
from __future__ import annotations

import json

import numpy as np
import pandas as pd
import pytest

from tissueresolve import cli


# --- deterministic synthetic fixtures ---------------------------------------

N_GENES = 60
CELL_TYPES = ["A", "B", "C"]


def _profiles(rng):
    prof = {}
    for k, ct in enumerate(CELL_TYPES):
        base = rng.gamma(1.0, 1.0, N_GENES)
        base[k * 15:(k + 1) * 15] *= 8.0  # block markers → separable
        prof[ct] = base
    return prof


def _make_reference(tmp_path):
    import anndata as ad
    rng = np.random.default_rng(0)
    prof = _profiles(rng)
    rows, labels = [], []
    for ct in CELL_TYPES:
        for _ in range(30):
            rows.append(rng.poisson(prof[ct] * 5))
            labels.append(ct)
    adata = ad.AnnData(X=np.asarray(rows, dtype=float))
    adata.var_names = [f"G{i:03d}" for i in range(N_GENES)]
    adata.obs["cell_type"] = labels
    ref = tmp_path / "ref.h5ad"
    adata.write(ref)
    return ref


def _make_bulk(tmp_path):
    rng = np.random.default_rng(1)
    prof = _profiles(np.random.default_rng(0))
    genes = [f"G{i:03d}" for i in range(N_GENES)]
    true = rng.dirichlet(np.ones(3), size=4)
    mat = np.zeros((N_GENES, 4))
    for j in range(4):
        mix = sum(true[j, k] * prof[ct] for k, ct in enumerate(CELL_TYPES))
        mat[:, j] = rng.poisson(mix * 1000)
    df = pd.DataFrame(mat, index=genes, columns=[f"s{j}" for j in range(4)])
    q = tmp_path / "bulk.tsv"
    df.to_csv(q, sep="\t")
    return q


def _make_visium(tmp_path):
    import anndata as ad
    import scipy.sparse as sp
    rng = np.random.default_rng(2)
    prof = {ct: (p / p.sum()) for ct, p in _profiles(np.random.default_rng(0)).items()}
    genes = [f"G{i:03d}" for i in range(N_GENES)]
    rpc = 6
    rows, cols = [], []
    for r in range(rpc):
        for c in range(rpc):
            rows.append(r)
            cols.append(c * 2 + (r % 2))
    N = len(rows)
    true = rng.dirichlet(np.ones(3), size=N)
    X = np.zeros((N, N_GENES))
    for i in range(N):
        mix = sum(true[i, k] * prof[ct] for k, ct in enumerate(CELL_TYPES))
        X[i] = rng.poisson(mix * 800)
    adata = ad.AnnData(X=sp.csr_matrix(X))
    adata.var_names = genes
    adata.obs_names = [f"spot{i}" for i in range(N)]
    adata.obs["array_row"] = np.asarray(rows, dtype=np.int32)
    adata.obs["array_col"] = np.asarray(cols, dtype=np.int32)
    adata.obs["total_counts"] = np.asarray(X.sum(1)).ravel()
    adata.obsm["spatial"] = np.c_[cols, rows].astype(float)
    q = tmp_path / "visium.h5ad"
    adata.write(q)
    return q


def _assert_metadata(out):
    assert (out / "analysis_plan.json").exists()
    assert (out / "run_metadata.json").exists()
    plan = json.loads((out / "analysis_plan.json").read_text())
    assert plan["mode"] in ("bulk", "spatial")


# --- bulk happy path --------------------------------------------------------

def test_bulk_end_to_end_writes_expected_outputs(tmp_path):
    ref = _make_reference(tmp_path)
    bulk = _make_bulk(tmp_path)
    out = tmp_path / "bulk_out"
    rc = cli.run(["--reference", str(ref), "--query", str(bulk),
                  "--out", str(out), "--mode", "bulk", "--preset", "quick"])
    assert rc == 0
    # metadata JSON
    _assert_metadata(out)
    # prediction table
    props = out / "deconv" / "proportions.tsv"
    assert props.exists()
    df = pd.read_csv(props, sep="\t", index_col=0, comment="#")
    assert df.shape[0] == 4 and df.shape[1] == 3
    # rows are (near) simplex
    assert np.allclose(df.sum(axis=1).to_numpy(), 1.0, atol=1e-3)
    # QC table present
    assert (out / "qc" / "metadata.json").exists()
    # report bundle produced by `run` (methods.txt, warnings.json, report.html)
    assert (out / "methods.txt").exists()
    warns = json.loads((out / "warnings.json").read_text())
    assert any(w["category"] == "estimate_type" for w in warns)
    report = (out / "report.html").read_text()
    assert "mRNA" in report  # estimate-type-aware, predictions populated


# --- spatial happy path -----------------------------------------------------

def test_spatial_end_to_end_writes_expected_outputs(tmp_path):
    ref = _make_reference(tmp_path)
    visium = _make_visium(tmp_path)
    out = tmp_path / "spatial_out"
    rc = cli.run(["--reference", str(ref), "--query", str(visium),
                  "--out", str(out), "--mode", "spatial", "--preset", "quick"])
    assert rc == 0
    _assert_metadata(out)
    props = out / "deconv" / "proportions.tsv"
    assert props.exists()
    df = pd.read_csv(props, sep="\t", index_col=0, comment="#")
    assert df.shape[0] == 36
    # QC + Moran's I present
    assert (out / "qc" / "metadata.json").exists()
    assert (out / "qc" / "morans_i.tsv").exists()
    # report bundle produced by `run`
    assert (out / "methods.txt").exists()
    warns = json.loads((out / "warnings.json").read_text())
    assert any(w["category"] == "estimate_type" for w in warns)
    assert (out / "report.html").exists()


# --- report generation independent of a run --------------------------------

def test_report_generation_from_results_dir(tmp_path):
    """Report generation works from a report-shaped results directory."""
    from tissueresolve.report import generate_report

    rdir = tmp_path / "results"
    (rdir / "tables").mkdir(parents=True)
    (rdir / "figures").mkdir(parents=True)
    pd.DataFrame(np.eye(3), index=["s0", "s1", "s2"],
                 columns=["A", "B", "C"]).to_csv(
        rdir / "bulk_estimated_proportions.tsv", sep="\t")
    pd.DataFrame({"value": [40, 8, 3]},
                 index=["n_cells", "n_genes", "n_cell_types"]).to_csv(
        rdir / "tables" / "reference_summary.tsv", sep="\t")
    (rdir / "tables" / "warnings.json").write_text(
        json.dumps({"qc_recommendations": ["sample s2 has low R2"]}))
    out = generate_report("bulk", rdir, rdir / "report.html")
    assert out.exists()
    doc = out.read_text()
    assert "Methods" in doc
    assert "low R" in doc  # warning surfaced, not hidden
