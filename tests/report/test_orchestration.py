"""
Tests for the canonical report orchestration layer (report/orchestration.py).

Verifies the single entry point renders the unified single-page shell from
both a results directory and an in-memory result, and that the public report
package + deprecated html shims all route through it.
"""
from __future__ import annotations

import numpy as np
import pandas as pd

from tissueresolve.report import orchestration as O
from tissueresolve.report.orchestration import render_sections


def test_render_sections_writes_unified_shell(tmp_path):
    out = render_sections(tmp_path / "r.html", "Title",
                          [("Reference quality", "<p>ref</p>"),
                           ("Predictions", "<p>preds</p>")])
    doc = out.read_text()
    assert out.name == "r.html"
    # unified single-page shell: sidebar nav + section anchors
    assert "sidebar-nav" in doc
    assert "Reference quality" in doc and "Predictions" in doc
    assert "id='reference-quality-0'" in doc  # deterministic anchor


def _bulk_results_dir(tmp_path):
    rdir = tmp_path / "bulk"
    (rdir / "tables").mkdir(parents=True)
    (rdir / "figures").mkdir(parents=True)
    pd.DataFrame(np.eye(3), index=["s0", "s1", "s2"],
                 columns=["A", "B", "C"]).to_csv(
        rdir / "bulk_estimated_proportions.tsv", sep="\t")
    pd.DataFrame({"value": [40, 8, 3]},
                 index=["n_cells", "n_genes", "n_cell_types"]).to_csv(
        rdir / "tables" / "reference_summary.tsv", sep="\t")
    return rdir


def test_generate_report_from_results_dir(tmp_path):
    rdir = _bulk_results_dir(tmp_path)
    out = O.generate_report("bulk", rdir)
    assert out.exists() and out.name == "report.html"
    doc = out.read_text()
    assert "sidebar-nav" in doc            # unified shell
    assert "Single-cell reference quality" in doc
    assert "cell fractions" in doc          # estimate-type caveat


def test_in_memory_result_requires_out_path(tmp_path):
    import pytest

    class _Fake:
        pass

    with pytest.raises(ValueError):
        O.generate_report("bulk", _Fake())  # no out path → clear error


def test_bad_modality_rejected(tmp_path):
    import pytest
    with pytest.raises(ValueError):
        O.generate_report("nonsense", tmp_path)


def test_public_package_entry_routes_through_orchestration():
    import tissueresolve.report as R
    from tissueresolve.report import orchestration
    assert R.generate_report is orchestration.generate_report


def test_html_shims_delegate(tmp_path):
    """The deprecated html.generate_* shims still work (via orchestration)."""
    from tissueresolve.report import html
    rdir = _bulk_results_dir(tmp_path)
    out = html.generate_bulk_report(str(rdir), out=rdir / "r.html")
    assert out.exists()
    assert "sidebar-nav" in out.read_text()  # rendered via the unified shell
