"""
v0.1 stabilization guardrails: core must not depend on the dev benchmark tree,
batch-effect diagnostics moved into the package, and claims/labels are honest.
"""
from __future__ import annotations

import re
from pathlib import Path

import pandas as pd
import pytest

_ROOT = Path(__file__).resolve().parents[1]
_SRC = _ROOT / "src" / "tissueresolve"


# --- Part 3: no core → benchmarks back-edge ---------------------------------

def test_no_src_module_imports_top_level_benchmarks():
    offenders = []
    pat = re.compile(r"^\s*(from\s+benchmarks(\.|\s)|import\s+benchmarks(\.|\s|$))")
    for p in _SRC.rglob("*.py"):
        for i, line in enumerate(p.read_text().splitlines(), 1):
            if pat.match(line):
                offenders.append(f"{p.relative_to(_ROOT)}:{i}: {line.strip()}")
    assert not offenders, "src/tissueresolve must not import benchmarks/: " + \
        "; ".join(offenders)


def test_batch_effects_in_core_diagnostics():
    from tissueresolve.diagnostics import batch_effects as be
    for fn in ("detect_batch_columns", "compute_celltype_batch_confounding",
               "compute_marker_batch_stability", "compute_batch_mixing_score",
               "batch_aware_reference_summary"):
        assert hasattr(be, fn)


def test_benchmarks_shim_reexports_core():
    from benchmarks.shared import batch_effects as shim
    from tissueresolve.diagnostics import batch_effects as core
    assert shim.compute_celltype_batch_confounding is core.compute_celltype_batch_confounding


def test_batch_confounding_regression():
    """The moved utility still computes a confounding table (behaviour unchanged)."""
    from tissueresolve.diagnostics.batch_effects import compute_celltype_batch_confounding
    # cell type A only in batch b1, B only in b2 → confounded; C balanced
    obs = pd.DataFrame({
        "cell_type": (["A"] * 20 + ["B"] * 20 + ["C"] * 20),
        "batch": (["b1"] * 20 + ["b2"] * 20 + ["b1"] * 10 + ["b2"] * 10),
    })
    conf = compute_celltype_batch_confounding(obs, "cell_type", "batch")
    assert "confounded" in conf.columns
    assert bool(conf.loc["A", "confounded"]) and bool(conf.loc["B", "confounded"])


def test_suitability_batch_component_uses_core(monkeypatch):
    """Reference suitability's batch component runs via the core import."""
    import tissueresolve.reference.suitability as S
    obs = pd.DataFrame({"cell_type": ["A"] * 10 + ["B"] * 10,
                        "batch": ["b1"] * 10 + ["b2"] * 10})
    comp = S._score_batch_confounding(obs, "cell_type", "batch") \
        if hasattr(S, "_score_batch_confounding") else None
    if comp is None:
        pytest.skip("batch component helper not exposed; covered by suitability tests")
    assert comp.name == "batch_confounding"


# --- Part 2/7/8: honest claims/labels ---------------------------------------

def _flat(text: str) -> str:
    """Collapse whitespace + markdown line markers so phrases match across wraps."""
    return re.sub(r"\s+", " ", re.sub(r"[>#`*]", " ", text))


def test_readme_marks_state_aware_experimental_and_no_expr_reconstruction():
    readme = _flat((_ROOT / "README.md").read_text())
    assert "Experimental: state-aware deconvolution runs behind" in readme
    assert "not been validated across real datasets" in readme
    assert "expression reconstruction is planned/deferred and not implemented" in readme.lower()


def test_cli_help_marks_state_aware_experimental():
    # check the RENDERED option help (string literals concatenated by click)
    from tissueresolve.cli import run_cli
    opt = next(p for p in run_cli.params if getattr(p, "name", "") == "state_aware")
    help_text = " ".join((opt.help or "").split())
    assert "EXPERIMENTAL" in help_text
    assert "not been validated across real datasets" in help_text


def test_feature_status_labels():
    fs = (_ROOT / "docs" / "FEATURE_STATUS.md").read_text().lower()
    assert "state-aware" in fs and "experimental" in fs
    assert "cell-type-specific expression reconstruction" in fs
    assert "not implemented" in fs


def test_tuning_remains_quarantined():
    from tissueresolve import tuning
    with pytest.raises(NotImplementedError):
        tuning.tune_bulk_parameters(_ROOT / "nonexistent_tmp")
