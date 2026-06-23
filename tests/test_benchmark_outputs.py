"""Tests for the external-benchmark metadata + a loader/summarizer for tool outputs.

Offline and independent of git-ignored generated metric files: the registry/manifest
metadata under benchmarks/external_tools/ are small versioned files; the loader is
tested on a synthetic in-memory metrics frame.
"""
from pathlib import Path

import numpy as np
import pandas as pd
import pytest

REPO = Path(__file__).resolve().parents[1]
EXT = REPO / "benchmarks" / "external_tools"


def summarize_external_metrics(df: pd.DataFrame, lower_is_better=None) -> pd.DataFrame:
    """Best-method-per-metric summary over a (method × metric) frame.

    Mirrors the benchmark's summary logic; pure + reusable so it can be tested
    without depending on generated TSVs.
    """
    lower = set(lower_is_better or {"rmse", "spillover", "false_positive_rate", "jsd"})
    g = df.groupby("method").mean(numeric_only=True)
    rows = []
    for metric in g.columns:
        asc = any(k in metric for k in lower)
        order = g[metric].sort_values(ascending=asc)
        rows.append({"metric": metric, "best_method": order.index[0],
                     "best_value": float(order.iloc[0])})
    return pd.DataFrame(rows)


# --------------------------------------------------------------------- loader
def test_loader_summarizes_external_metrics():
    df = pd.DataFrame({
        "method": ["A", "A", "B", "B"],
        "fine_rmse": [0.10, 0.12, 0.20, 0.22],
        "false_positive_rate": [0.30, 0.28, 0.05, 0.06],
        "fine_pearson": [0.80, 0.82, 0.50, 0.52],
    })
    s = summarize_external_metrics(df).set_index("metric")
    assert s.loc["fine_rmse", "best_method"] == "A"             # lower is better
    assert s.loc["false_positive_rate", "best_method"] == "B"   # lower is better
    assert s.loc["fine_pearson", "best_method"] == "A"          # higher is better


# --------------------------------------------------------------------- registry / manifest
def test_tool_registry_exists_and_parses():
    reg = EXT / "tool_registry.tsv"
    assert reg.exists(), "tool_registry.tsv must exist"
    df = pd.read_csv(reg, sep="\t")
    for col in ("tool", "modality", "installed", "version"):
        assert col in df.columns
    # the validated/known tools appear
    tools = set(df["tool"])
    assert {"TissueResolve_hierarchical", "NNLS_baseline", "MuSiC", "cell2location"} <= tools


def test_run_manifest_records_skipped_and_executed_with_reasons():
    man = EXT / "run_manifest.tsv"
    assert man.exists(), "run_manifest.tsv must exist (created by the bulk benchmark)"
    df = pd.read_csv(man, sep="\t")
    assert {"tool", "executed_this_session", "reason"} <= set(df.columns)
    # every row has a non-empty reason (no silent skips)
    assert df["reason"].astype(str).str.len().gt(0).all()
    # at least one fresh internal execution and at least one external tool documented
    assert bool(df["executed_this_session"].astype(str).str.lower().eq("true").any())
    assert len(df) >= 5


def test_label_mapping_and_input_metadata_present():
    for name in ("label_mapping.tsv", "input_manifest.tsv",
                 "normalization_summary.tsv", "gene_overlap_summary.tsv"):
        p = EXT / name
        assert p.exists(), f"{name} must exist"
        assert pd.read_csv(p, sep="\t").shape[0] >= 1


def test_no_raw_data_paths_in_versioned_metadata():
    # small versioned metadata must not embed raw data file paths (.h5ad/.h5/.zarr)
    for p in EXT.glob("*.tsv"):
        text = p.read_text().lower()
        assert ".h5ad" not in text and ".zarr" not in text
        assert "/data/" not in text


def test_external_tools_inputs_gitignore_excludes_large_files():
    gi = EXT / "inputs" / ".gitignore"
    assert gi.exists()
    body = gi.read_text()
    for pat in (".h5ad", ".zarr"):
        assert pat in body
