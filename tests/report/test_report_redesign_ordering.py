"""
Offline tests for the redesigned report's decision-workflow ordering and the
hierarchical-palette guarantees (redesign Parts 2-3, 6, 17).

These exercise the report-builder *contract* with synthetic sections (no real
data, no network), plus the palette relationships.
"""
from __future__ import annotations

import pandas as pd

from tissueresolve.report import components as C
from tissueresolve.report.unified import Section, build_unified_report


# anchors in the redesigned decision-workflow order
_ORDER = [
    ("summary", "1. Executive decision summary"),
    ("reference", "2. Reference quality"),
    ("signature", "3. Signature quality & hierarchy"),
    ("input", "4. Input data quality"),
    ("bulk", "5. Bulk deconvolution"),
    ("spatial", "6. Spatial deconvolution"),
    ("hierarchical", "7. Hierarchical broad→fine deconvolution"),
    ("resolution", "8. Resolution, separability & spillover"),
    ("benchmark", "9. Benchmark comparison"),
    ("warnings", "10. Warnings & limitations"),
]


def _report(tmp_path):
    secs = [Section(a, t, C.metric_grid({"x": 1})) for a, t in _ORDER]
    return build_unified_report(tmp_path / "report.html", secs).read_text()


def test_qc_sections_appear_before_results(tmp_path):
    """Reference, signature and input QC must precede bulk/spatial results."""
    html = _report(tmp_path)
    for qc in ("reference", "signature", "input"):
        for result in ("bulk", "spatial"):
            assert html.index(f"id='{qc}'") < html.index(f"id='{result}'"), \
                f"{qc} must come before {result}"


def test_signature_quality_before_results(tmp_path):
    html = _report(tmp_path)
    assert html.index("id='signature'") < html.index("id='bulk'")
    assert html.index("id='signature'") < html.index("id='spatial'")


def test_warnings_section_present_and_navigable(tmp_path):
    html = _report(tmp_path)
    assert "href='#warnings'" in html and "id='warnings'" in html


# --- hierarchical palette (Part 6) ------------------------------------------

def _color_map():
    from tissueresolve.plotting import palette as p
    fine = ["CD8 T cell", "CD4 T cell", "Treg", "NK cell",
            "macrophage", "monocyte", "Other"]
    mapping = {"CD8 T cell": "T/NK", "CD4 T cell": "T/NK", "Treg": "T/NK",
               "NK cell": "T/NK", "macrophage": "Myeloid", "monocyte": "Myeloid"}
    return p.build_hierarchical_color_map(fine, mapping)


def _as_dict(df):
    return dict(zip(df["fine_cell_type"], df["color_hex"]))


def _hex_to_rgb(h):
    h = h.lstrip("#")
    return tuple(int(h[i:i + 2], 16) for i in (0, 2, 4))


def _dist(a, b):
    return sum((x - y) ** 2 for x, y in zip(_hex_to_rgb(a), _hex_to_rgb(b))) ** 0.5


def test_broad_families_get_distinct_base_colors():
    df = _color_map()
    broad = df[df["color_role"] == "broad"]
    base = dict(zip(broad["broad_cell_type"], broad["color_hex"]))
    # the two broad family base colors must be visibly distinct
    assert _dist(base["T/NK"], base["Myeloid"]) > 30


def test_fine_labels_share_family_color_region():
    cm = _as_dict(_color_map())
    # within-family spread should be smaller than between-family spread
    within = _dist(cm["CD8 T cell"], cm["CD4 T cell"])
    between = _dist(cm["CD8 T cell"], cm["macrophage"])
    assert within < between


def test_palette_deterministic():
    assert _as_dict(_color_map()) == _as_dict(_color_map())


def test_other_is_neutral_grey():
    cm = _as_dict(_color_map())
    r, g, b = _hex_to_rgb(cm["Other"])
    # grey => channels close together
    assert max(r, g, b) - min(r, g, b) < 25
