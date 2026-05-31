"""
Tests for the family-aware palette and label helpers (offline).
"""
from __future__ import annotations

import pandas as pd
import pytest

from tissueresolve.plotting.palette import (
    OTHER_COLOR, UNRESOLVED_COLOR, assign_family_palette,
    infer_cell_type_family, save_color_map, shorten_cell_type_label, wrap_label,
)


@pytest.mark.parametrize("label,family", [
    ("CD8-positive, alpha-beta T cell", "T/NK"),
    ("regulatory T cell", "T/NK"),
    ("natural killer cell", "T/NK"),
    ("macrophage", "myeloid"),
    ("conventional dendritic cell", "myeloid"),
    ("capillary endothelial cell", "endothelial"),
    ("luminal epithelial cell of mammary gland", "epithelial"),
    ("basal cell", "epithelial"),
    ("fibroblast", "stromal"),
    ("pericyte", "mural"),
    ("vascular associated smooth muscle cell", "mural"),
    ("plasma cell", "B/plasma"),
    ("adipocyte", "adipocyte"),
    ("some unknown thing", "other"),
])
def test_infer_cell_type_family(label, family):
    assert infer_cell_type_family(label) == family


def test_palette_consistent_and_grouped():
    cts = ["CD4 T cell", "CD8 T cell", "macrophage", "monocyte",
           "luminal epithelial cell", "fibroblast", "Other"]
    a = assign_family_palette(cts)
    b = assign_family_palette(cts)
    assert a == b                                  # stable
    assert a["Other"] == OTHER_COLOR               # light grey
    # T/NK members share the T/NK ramp (purple family), distinct from myeloid.
    from tissueresolve.plotting.palette import FAMILY_RAMPS
    assert a["CD4 T cell"] in FAMILY_RAMPS["T/NK"]
    assert a["CD8 T cell"] in FAMILY_RAMPS["T/NK"]
    assert a["CD4 T cell"] != a["CD8 T cell"]      # related but distinct shades
    assert a["macrophage"] in FAMILY_RAMPS["myeloid"]


def test_unresolved_is_dark_grey():
    cm = assign_family_palette(["unresolved_T cell", "myeloid cells"])
    assert cm["unresolved_T cell"] == UNRESOLVED_COLOR
    assert cm["myeloid cells"] == UNRESOLVED_COLOR  # family-level label


def test_shorten_and_wrap():
    long = "CD8-positive, alpha-beta memory T cell"
    short = shorten_cell_type_label(long)
    assert len(short) < len(long)
    wrapped = wrap_label("vascular associated smooth muscle cell", width=18)
    assert "<br>" in wrapped


def test_save_color_map(tmp_path):
    cm = assign_family_palette(["CD8 T cell", "macrophage", "Other"])
    p = save_color_map(cm, tmp_path)
    assert p.exists()
    df = pd.read_csv(p, sep="\t")
    assert {"cell_type", "family", "short_label", "color"}.issubset(df.columns)
