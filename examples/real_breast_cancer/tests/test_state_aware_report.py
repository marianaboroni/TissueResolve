"""State-aware report subsection checks (skip when the report isn't generated)."""
from __future__ import annotations

from pathlib import Path

import pytest

REPORT = Path("examples/real_breast_cancer/outputs/report.html")


@pytest.fixture
def rep(load_script):
    return load_script("07_generate_reports.py")


def test_subsection_helper_states_status_when_not_run(rep, monkeypatch, tmp_path):
    """When no signature outputs exist, the helper states the experimental path
    was not run and that no third-level state labels are available."""
    import _harness as H
    monkeypatch.setattr(H, "OUTPUTS_DIR", tmp_path)   # empty outputs dir
    from tissueresolve.report import components as C
    html = rep._state_aware_subsection(C)
    assert "State-aware / multi-granularity deconvolution" in html
    assert "experimental" in html.lower()
    assert "not run for this dataset" in html
    assert "no third-level cell-state labels" in html


@pytest.mark.skipif(not REPORT.exists(), reason="report not generated")
def test_generated_report_has_state_aware_section():
    html = REPORT.read_text()
    assert "State-aware / multi-granularity deconvolution" in html
    # honest about state-label availability + experimental status
    assert "no third-level cell-state labels are" in html
    assert "experimental" in html.lower()
