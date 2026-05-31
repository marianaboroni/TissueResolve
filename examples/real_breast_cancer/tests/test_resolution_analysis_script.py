"""
Offline test for scripts/08_resolution_analysis.py.

Builds a tiny saved reference with confusable families + tiny fine bulk/spatial
estimates, runs the analysis, and checks the family-level outputs (mass
preserved, fine estimates untouched, recommended merges written).
"""
from __future__ import annotations

import numpy as np
import pandas as pd


def _make_saved_reference(tmp_path):
    from tissueresolve.results import ReferenceSignature

    G = 30
    base = np.full(G, 5.0, dtype=np.float32)

    def block(lo, hi, val):
        v = base.copy()
        v[lo:hi] = val
        return v

    profiles = {
        "CD4-positive helper T cell": block(0, 10, 1000),
        "CD8-positive, alpha-beta T cell": block(0, 10, 1005),
        "macrophage": block(10, 20, 1000),
        "monocyte": block(10, 20, 1004),
        "capillary endothelial cell": block(20, 30, 1000),
    }
    cts = list(profiles)
    R = np.vstack([profiles[c] for c in cts]).astype(np.float32)
    ref = ReferenceSignature(
        gene_names=[f"G{i:03d}" for i in range(G)], cell_types=cts,
        R_cpm=R, R_log=np.log1p(R).astype(np.float32),
        phi_g=np.full(G, 5.0, dtype=np.float32),
        n_cells_per_type={c: 100 for c in cts})
    d = tmp_path / "reference"
    ref.save(d)
    return d, cts


def test_08_help_exposes_resolution_mode(load_script):
    m08 = load_script("08_resolution_analysis.py")
    # argparse exposes --resolution-mode
    import argparse
    import contextlib
    import io

    buf = io.StringIO()
    with contextlib.redirect_stdout(buf):
        try:
            m08.main(["--help"])
        except SystemExit:
            pass
    assert "--resolution-mode" in buf.getvalue()


def test_08_end_to_end(load_script, harness, monkeypatch, tmp_path):
    ref_dir, cts = _make_saved_reference(tmp_path)
    bulk_dir = tmp_path / "bulk"
    spatial_dir = tmp_path / "spatial"
    res_dir = tmp_path / "resolution"
    bulk_dir.mkdir(); spatial_dir.mkdir()
    monkeypatch.setattr(harness, "SAVED_REFERENCE_DIR", ref_dir)
    monkeypatch.setattr(harness, "OUT_BULK_DIR", bulk_dir)
    monkeypatch.setattr(harness, "OUT_SPATIAL_DIR", spatial_dir)
    monkeypatch.setattr(harness, "OUT_RESOLUTION_DIR", res_dir)
    monkeypatch.setattr(harness, "ALL_DIRS", [res_dir])

    rng = np.random.default_rng(0)
    bulk = pd.DataFrame(rng.dirichlet(np.ones(len(cts)), size=3),
                        index=["s0", "s1", "s2"], columns=cts)
    bulk.to_csv(bulk_dir / "bulk_estimated_proportions.tsv", sep="\t")
    spatial = pd.DataFrame(rng.dirichlet(np.ones(len(cts)), size=4),
                           index=[f"sp{i}" for i in range(4)], columns=cts)
    spatial.to_csv(spatial_dir / "spatial_spot_proportions.tsv", sep="\t")

    m08 = load_script("08_resolution_analysis.py")
    rc = m08.main([])
    assert rc == 0

    for f in ("pairwise_separability.tsv", "cell_type_families.tsv",
              "recommended_merges.tsv", "unresolved_families.tsv",
              "bulk_family_proportions.tsv", "spatial_family_proportions.tsv",
              "resolution_summary.md", "resolution_mapping.json"):
        assert (res_dir / f).exists(), f

    # Family aggregation preserves total mass; fine file untouched.
    fam = pd.read_csv(res_dir / "bulk_family_proportions.tsv", sep="\t", index_col=0)
    np.testing.assert_allclose(fam.sum(axis=1).to_numpy(), 1.0, atol=1e-6)
    fine_after = pd.read_csv(bulk_dir / "bulk_estimated_proportions.tsv",
                             sep="\t", index_col=0)
    pd.testing.assert_frame_equal(fine_after, bulk)

    merges = pd.read_csv(res_dir / "recommended_merges.tsv", sep="\t")
    assert "T/NK lymphocytes" in set(merges["family_name"])

    import json
    mapping = json.loads((res_dir / "resolution_mapping.json").read_text())
    assert mapping["fine_estimates_overwritten"] is False
    assert mapping["merge_stage"] == "post_hoc_aggregation"
