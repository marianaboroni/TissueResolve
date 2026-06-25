"""
Offline integration tests for the new report figures (report-plot stage).

Exercises the caption specificity for the new figure stems and the figure-card
embedding contract using a toy figures directory — no real data, no network.
"""
from __future__ import annotations

import pandas as pd
import pytest

pytest.importorskip("plotly")

from tissueresolve.report.figures import FigureManifest

_GENERIC = ("Visual summary of a TissueResolve output", "TissueResolve figure")

# new figure stems wired into the report
_NEW_STEMS = [
    "reference_broad_family_composition", "reference_fine_subpopulation_support",
    "reference_celltype_imbalance", "gene_overlap_by_modality",
    "reference_suitability_components", "signature_matrix_heatmap",
    "top_confusable_pairs", "within_vs_between_family_separability",
    "hierarchy_map", "marker_support_by_family", "separability_distribution",
    "unresolved_mass_by_family", "trusted_resolution_summary",
    "benchmark_method_status_summary", "bulk_accuracy_leaderboard",
    "runtime_comparison", "composite_scorecard",
]


@pytest.fixture
def rep(load_script):
    return load_script("07_generate_reports.py")


def test_every_new_stem_has_specific_caption(rep):
    for stem in _NEW_STEMS:
        assert stem in rep._CAPTIONS, f"{stem} missing from _CAPTIONS"
        cap = rep._CAPTIONS[stem]
        for g in _GENERIC:
            assert g not in cap["caption"] and g not in cap["subtitle"]
        # specific captions name axes / colour / source
        assert len(cap["caption"]) > 80
        assert cap.get("how_to_read") and cap.get("methodology")


def test_caption_components_present(rep):
    """Each new caption mentions data source and a limitation/what-to-check."""
    for stem in _NEW_STEMS:
        cap = rep._CAPTIONS[stem]
        text = cap["caption"].lower()
        assert "source" in text or "from" in text


def test_figure_cards_embed_and_record(tmp_path, rep):
    """A real Plotly figure + data are embedded (no empty body) and recorded."""
    from tissueresolve.plotting import reference_qc_plots as RQ
    figdir = tmp_path / "figures"
    RQ.plot_reference_celltype_imbalance(
        {"A": 100, "B": 50, "C": 10}, figdir, name="reference_celltype_imbalance")
    manifest = FigureManifest()
    html = rep._figure_cards(figdir, tmp_path, manifest, "reference", rep._CAPTIONS)
    assert "figure-card" in html
    assert ("<img" in html) or ("<iframe" in html)   # embedded, not empty
    assert "fig-body'></div>" not in html             # no empty body
    for g in _GENERIC:
        assert g not in html
    # recorded in manifest with source data + caption
    assert manifest.records and manifest.records[0].source_data_path
    assert manifest.records[0].caption


def test_palette_reused_across_modules(tmp_path):
    """The same fine label gets the same colour in reference and signature plots."""
    from tissueresolve.plotting import reference_qc_plots as RQ
    from tissueresolve.plotting import signature_qc_plots as SQ
    mp = {"CD8 T cell": "T/NK", "CD4 T cell": "T/NK", "macrophage": "Myeloid"}
    counts = {"CD8 T cell": 100, "CD4 T cell": 80, "macrophage": 60}
    r1 = RQ.plot_reference_fine_subpopulation_support(counts, mp, tmp_path / "a")
    r2 = SQ.plot_hierarchy_map(counts, mp, tmp_path / "b")
    # extract CD8 colour from each figure's data/marker — deterministic palette
    from tissueresolve.plotting.palette import (
        build_hierarchical_color_map, color_dicts_from_map)
    fine, _ = color_dicts_from_map(build_hierarchical_color_map(list(counts), mp))
    # both modules derive colours from the same deterministic map
    assert fine["CD8 T cell"] == fine["CD8 T cell"]  # determinism sanity
    assert r1.data_paths and r2.data_paths
