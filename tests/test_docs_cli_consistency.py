"""
CLI/documentation consistency and documentation-honesty regression tests.

Offline.  Guards the audit fixes in ``docs/CLI_DOCUMENTATION_AUDIT.md`` and
``docs/DOCUMENTATION_HONESTY_AUDIT.md`` so they cannot silently regress.
"""
from __future__ import annotations

import re
from pathlib import Path

from click.testing import CliRunner

from tissueresolve import cli

_ROOT = Path(__file__).resolve().parents[1]


def _flat(text: str) -> str:
    return re.sub(r"\s+", " ", re.sub(r"[>#`*]", " ", text))


# --- CLI commands exist ------------------------------------------------------

def test_documented_cli_commands_exist():
    """Every top-level command the README advertises exists in the CLI."""
    r = CliRunner().invoke(cli.cli, ["--help"])
    assert r.exit_code == 0
    for cmd in ("run", "bulk", "spatial", "report", "info"):
        assert cmd in r.output, cmd


def test_run_has_mode_bulk_and_spatial():
    r = CliRunner().invoke(cli.cli, ["run", "--help"])
    assert r.exit_code == 0
    assert "--mode" in r.output
    # the documented bulk/spatial entrypoint
    assert "auto|bulk|spatial" in r.output or "bulk" in r.output


def test_bulk_run_stub_is_marked_not_implemented():
    """`tissueresolve bulk run` is a pseudocode/stub; it must fail clearly."""
    r = CliRunner().invoke(cli.cli, ["bulk", "run"])
    assert r.exit_code != 0
    assert "not yet implemented" in (r.output + str(r.stderr_bytes or b"")).lower() \
        or "not yet implemented" in r.output.lower()


def test_solver_and_resolution_choices_match_docs():
    r = CliRunner().invoke(cli.cli, ["run", "--help"])
    for tok in ("auto", "nnls", "weighted_nnls", "marker_nnls", "ridge_nnls",
                "ensemble_nnls", "pipeline"):
        assert tok in r.output, tok
    for tok in ("hierarchical", "flat", "none", "suggest"):
        assert tok in r.output, tok


# --- experimental flags labelled --------------------------------------------

def test_state_aware_flag_labelled_experimental():
    r = CliRunner().invoke(cli.cli, ["run", "--help"])
    assert "--state-aware" in r.output
    opt = next(p for p in cli.run_cli.params if getattr(p, "name", "") == "state_aware")
    help_text = " ".join((opt.help or "").split())
    assert "EXPERIMENTAL" in help_text


# --- README does not reference non-existent commands ------------------------

def test_readme_only_references_real_top_level_commands():
    readme = (_ROOT / "README.md").read_text()
    # collect `tissueresolve <word>` occurrences in fenced/inline code
    cmds = set(re.findall(r"tissueresolve\s+([a-z][a-z\-]*)", readme))
    valid = {"run", "bulk", "spatial", "report", "info", "combine-report"}
    unknown = cmds - valid
    assert not unknown, f"README references unknown commands: {unknown}"


def test_combine_report_command_exists():
    r = CliRunner().invoke(cli.cli, ["--help"])
    assert "combine-report" in r.output
    h = CliRunner().invoke(cli.cli, ["combine-report", "--help"])
    assert h.exit_code == 0
    for flag in ("--bulk-dir", "--spatial-dir", "--out"):
        assert flag in h.output, flag


def test_docs_mention_separate_output_dirs():
    for f in ("README.md", "docs/tutorial.md", "docs/output_interpretation.md"):
        txt = (_ROOT / f).read_text()
        assert "results/bulk" in txt and "results/spatial" in txt, f
        assert "combine-report" in txt, f


def test_benchmark_include_imported_flag_registered():
    """run_all.py must accept --include-imported (documented workflow)."""
    src = (_ROOT / "benchmarks" / "run_all.py").read_text()
    assert "--include-imported" in src


# --- documentation honesty ---------------------------------------------------

def test_no_claim_expression_reconstruction_implemented():
    fs = _flat((_ROOT / "docs" / "FEATURE_STATUS.md").read_text()).lower()
    assert "cell-type-specific expression reconstruction" in fs
    assert "not implemented" in fs
    readme = _flat((_ROOT / "README.md").read_text()).lower()
    assert "expression reconstruction is planned/deferred and not implemented" in readme


def test_spatial_benchmark_no_accuracy_without_ground_truth():
    bm = _flat((_ROOT / "docs" / "benchmarking.md").read_text()).lower()
    assert "no absolute accuracy is claimed without ground truth" in bm
    oi = _flat((_ROOT / "docs" / "output_interpretation.md").read_text()).lower()
    assert "spatial benchmark is not accuracy" in oi


def test_external_benchmark_only_executed_or_imported_ranked():
    fs = _flat((_ROOT / "docs" / "FEATURE_STATUS.md").read_text()).lower()
    assert "only executed or imported tools are ranked" in fs


def test_publication_status_is_alpha():
    for f in ("README.md", "docs/FEATURE_STATUS.md", "docs/dev/V0_1_SCOPE.md"):
        txt = _flat((_ROOT / f).read_text()).lower()
        assert "alpha" in txt and "early-access" in txt, f
        assert "not yet" in txt and "publication-ready" in txt, f


def test_estimate_type_wording_present():
    fs = _flat((_ROOT / "docs" / "FEATURE_STATUS.md").read_text()).lower()
    assert "rna-derived proportions, not absolute cell fractions" in fs
    assert "spot-level rna-derived composition, not single-cell labels" in fs
