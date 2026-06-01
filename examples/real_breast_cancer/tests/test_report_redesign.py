"""
Offline tests for the redesigned report generator's pure helpers
(redesign Parts 3, 12, 13, 17).

These load the numbered script module and exercise the caption resolver, the
decision-status derivation and the warning aggregation without any real data,
network access, or full output tree.
"""
from __future__ import annotations

import pytest


@pytest.fixture
def rep(load_script):
    return load_script("07_generate_reports.py")


# --- Part 12 / 13: no generic captions --------------------------------------

_GENERIC = ("Visual summary of a TissueResolve output", "TissueResolve figure")


def test_known_stem_caption_is_specific(rep):
    cap = rep._caption_for("bulk_composition_clustered_barplot", "bulk", rep._CAPTIONS)
    assert "RNA-derived" in cap["caption"]
    for g in _GENERIC:
        assert g not in cap["caption"] and g not in cap["subtitle"]


def test_unknown_stem_caption_is_specific_not_generic(rep):
    cap = rep._caption_for("some_brand_new_plot", "bulk", rep._CAPTIONS)
    # title-derived, never the old generic text
    assert "Some brand new plot" in cap["caption"]
    for g in _GENERIC:
        assert g not in cap["caption"] and g not in cap["subtitle"]


def test_caption_prefix_match_for_he_abundance(rep):
    cap = rep._caption_for("he_abundance_capillary_endothelial_cell", "spatial",
                           rep._CAPTIONS)
    assert "over H&E" in cap["caption"] or "H&E" in cap["caption"]


def test_no_caption_in_dict_is_generic(rep):
    """The legacy generic strings must not appear in any predefined caption."""
    for cap in rep._CAPTIONS.values():
        for g in _GENERIC:
            assert g not in cap.get("caption", "")
            assert g not in cap.get("subtitle", "")


# --- Part 3: decision statuses ----------------------------------------------

def test_decision_statuses_keys_present(rep):
    st = rep._decision_statuses({"overlap": "—", "solver": "auto"}, None, None)
    for k in ("reference", "hierarchy", "separability", "input",
              "bulk_available", "spatial_available", "benchmark_available"):
        assert k in st
    # with no suitability object, QC verdicts are UNKNOWN, availability is bool
    assert st["reference"] == "UNKNOWN"
    assert isinstance(st["bulk_available"], bool)


def test_decision_statuses_hierarchy_pass_when_mapping_present(rep):
    st = rep._decision_statuses({"overlap": 100, "solver": "auto"}, None,
                                {"CD8 T cell": "T/NK"})
    assert st["hierarchy"] == "PASS"
