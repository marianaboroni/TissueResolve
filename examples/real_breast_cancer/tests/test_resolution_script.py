"""
Offline test for scripts/06_resolution_spillover_analysis.py.

Builds a tiny saved reference, runs the analysis end-to-end (no network), and
checks the required outputs are written.
"""
from __future__ import annotations

import numpy as np
import pandas as pd


def _make_saved_reference(tmp_path):
    from tissueresolve.results import ReferenceSignature

    # Two near-identical types (T0/T1) + a distinct one (T2).
    G = 30
    R = np.full((3, G), 5.0, dtype=np.float32)
    R[0, :10] = 1000.0
    R[1, :10] = 980.0
    R[2, 20:] = 1000.0
    ref = ReferenceSignature(
        gene_names=[f"G{i:03d}" for i in range(G)],
        cell_types=["T0", "T1", "T2"],
        R_cpm=R, R_log=np.log1p(R).astype(np.float32),
        phi_g=np.full(G, 5.0, dtype=np.float32),
        n_cells_per_type={"T0": 100, "T1": 100, "T2": 100},
    )
    ref_dir = tmp_path / "reference"
    ref.save(ref_dir)
    return ref_dir


def test_06_imports(load_script):
    m06 = load_script("06_resolution_spillover_analysis.py")
    assert hasattr(m06, "main") and hasattr(m06, "run_analysis")


def test_06_end_to_end_writes_outputs(load_script, harness, monkeypatch, tmp_path):
    ref_dir = _make_saved_reference(tmp_path)
    out_dir = tmp_path / "resolution"
    monkeypatch.setattr(harness, "SAVED_REFERENCE_DIR", ref_dir)
    monkeypatch.setattr(harness, "OUT_RESOLUTION_DIR", out_dir)
    monkeypatch.setattr(harness, "ALL_DIRS", [out_dir])
    # No bulk/spatial outputs present → annotation step is skipped cleanly.
    monkeypatch.setattr(harness, "OUT_BULK_DIR", tmp_path / "nobulk")
    monkeypatch.setattr(harness, "OUT_SPATIAL_DIR", tmp_path / "nospatial")

    m06 = load_script("06_resolution_spillover_analysis.py")
    rc = m06.main(["--spillover-method", "expression"])
    assert rc == 0

    for fname in ("cell_type_families.tsv", "pairwise_resolvability.tsv",
                  "spillover_matrix.tsv", "spillover_risk_by_celltype.tsv",
                  "pairwise_spillover_report.tsv", "recommended_merges.tsv",
                  "unresolved_families.tsv", "resolution_summary.md"):
        assert (out_dir / fname).exists(), fname

    # Spillover matrix is square with rows summing to ~1.
    M = pd.read_csv(out_dir / "spillover_matrix.tsv", sep="\t", index_col=0)
    assert M.shape == (3, 3)
    np.testing.assert_allclose(M.to_numpy().sum(axis=1), 1.0, atol=1e-6)

    # The near-identical T0/T1 pair should be flagged unresolved/poorly-resolved.
    unresolved = pd.read_csv(out_dir / "unresolved_families.tsv", sep="\t")
    assert len(unresolved) >= 1


def test_06_deconvolution_method(load_script, harness, monkeypatch, tmp_path):
    ref_dir = _make_saved_reference(tmp_path)
    out_dir = tmp_path / "resolution"
    monkeypatch.setattr(harness, "SAVED_REFERENCE_DIR", ref_dir)
    monkeypatch.setattr(harness, "OUT_RESOLUTION_DIR", out_dir)
    monkeypatch.setattr(harness, "ALL_DIRS", [out_dir])
    monkeypatch.setattr(harness, "OUT_BULK_DIR", tmp_path / "nobulk")
    monkeypatch.setattr(harness, "OUT_SPATIAL_DIR", tmp_path / "nospatial")

    m06 = load_script("06_resolution_spillover_analysis.py")
    assert m06.main([]) == 0  # default deconvolution method
    M = pd.read_csv(out_dir / "spillover_matrix.tsv", sep="\t", index_col=0)
    np.testing.assert_allclose(M.to_numpy().sum(axis=1), 1.0, atol=1e-4)
