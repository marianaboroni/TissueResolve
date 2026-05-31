"""
Family-aware colour palette and label helpers for TissueResolve figures.

Goals
-----
* The **same cell type keeps the same colour** across every figure.
* Fine subtypes of the same broad family use **related shades**.
* ``Other`` is light grey; ``unresolved`` / family-level labels are dark grey.
* Long fine labels get short display forms (axes) and wrapped forms (hover).

The colour map can be persisted to ``figures/cell_type_color_map.tsv`` so it is
reproducible and shared across all panels of a report.
"""
from __future__ import annotations

import re
from pathlib import Path
from typing import Optional

import pandas as pd

__all__ = [
    "FAMILY_RAMPS",
    "OTHER_COLOR",
    "UNRESOLVED_COLOR",
    "infer_cell_type_family",
    "assign_family_palette",
    "shorten_cell_type_label",
    "wrap_label",
    "save_color_map",
]

OTHER_COLOR = "#d9d9d9"        # light grey
UNRESOLVED_COLOR = "#636363"   # dark grey (unresolved / family-level)

# Broad biological family → ordered shade ramp (dark → light).
FAMILY_RAMPS: dict[str, list[str]] = {
    "epithelial": ["#08519c", "#3182bd", "#6baed6", "#9ecae1", "#c6dbef"],
    "stromal":    ["#8c510a", "#bf812d", "#d8b365", "#dfc27d", "#f6e8c3"],
    "endothelial": ["#01665e", "#35978f", "#80cdc1", "#b2e2da", "#c7eae5"],
    "myeloid":    ["#006d2c", "#31a354", "#74c476", "#a1d99b", "#c7e9c0"],
    "T/NK":       ["#54278f", "#756bb1", "#9e9ac8", "#bcbddc", "#dadaeb"],
    "B/plasma":   ["#a50f15", "#de2d26", "#fb6a4a", "#fc9272", "#fcbba1"],
    "mural":      ["#808000", "#9a9a1f", "#b3b34d", "#c2c260", "#d4d48a"],
    "adipocyte":  ["#e6a817", "#f0c64b", "#f7d978"],
    "other":      [OTHER_COLOR],
}

# Keyword → family, checked in order (first match wins).  More specific first.
_FAMILY_KEYWORDS: list[tuple[str, str]] = [
    ("plasma", "B/plasma"), ("b cell", "B/plasma"), ("b-cell", "B/plasma"),
    ("regulatory t", "T/NK"), ("cd4", "T/NK"), ("cd8", "T/NK"),
    ("t cell", "T/NK"), ("t-cell", "T/NK"), ("nkt", "T/NK"),
    ("natural killer", "T/NK"), ("nk cell", "T/NK"), ("lymphocyte", "T/NK"),
    ("macrophage", "myeloid"), ("monocyte", "myeloid"), ("dendritic", "myeloid"),
    ("myeloid", "myeloid"), ("mast", "myeloid"), ("neutrophil", "myeloid"),
    ("granulocyte", "myeloid"),
    ("pericyte", "mural"), ("smooth muscle", "mural"), ("mural", "mural"),
    ("endothel", "endothelial"), ("vascular endothelial", "endothelial"),
    ("lymphatic", "endothelial"),
    ("fibroblast", "stromal"), ("stromal", "stromal"), ("caf", "stromal"),
    ("adipocyte", "adipocyte"), ("adipose", "adipocyte"),
    ("luminal", "epithelial"), ("basal", "epithelial"),
    ("myoepithelial", "epithelial"), ("epithelial", "epithelial"),
    ("carcinoma", "epithelial"), ("tumor", "epithelial"), ("tumour", "epithelial"),
]


def infer_cell_type_family(label: str) -> str:
    """Map a (fine) cell-type label to a broad biological family.

    Returns one of the keys of :data:`FAMILY_RAMPS` (``"other"`` if unknown).
    """
    low = str(label).lower()
    if low.startswith("unresolved") or low.startswith("family:"):
        return "other"
    for kw, fam in _FAMILY_KEYWORDS:
        if kw in low:
            return fam
    return "other"


def _is_family_level(label: str) -> bool:
    low = str(label).lower()
    return (low.startswith("unresolved") or " family" in low
            or low.endswith("cells") and any(
                k in low for k in ("lymphocytes", "myeloid cells",
                                   "endothelial cells", "epithelial cells",
                                   "mural cells")))


def assign_family_palette(cell_types: list[str]) -> dict[str, str]:
    """Assign a stable hex colour per cell type, grouped by biological family.

    Within a family, fine types get successive shades of that family's ramp
    (deterministic by sorted order).  ``Other`` → light grey; family-level /
    ``unresolved`` labels → dark grey.
    """
    by_family: dict[str, list[str]] = {}
    colors: dict[str, str] = {}
    for ct in cell_types:
        if str(ct).strip().lower() == "other":
            colors[ct] = OTHER_COLOR
            continue
        if _is_family_level(ct):
            colors[ct] = UNRESOLVED_COLOR
            continue
        fam = infer_cell_type_family(ct)
        by_family.setdefault(fam, []).append(ct)

    for fam, members in by_family.items():
        ramp = FAMILY_RAMPS.get(fam, FAMILY_RAMPS["other"])
        for i, ct in enumerate(sorted(members, key=str)):
            colors[ct] = ramp[i % len(ramp)]
    return colors


_STRIP_PATTERNS = [
    (r",?\s*alpha-beta", ""), (r"-positive", "+"), (r"-negative", "-"),
    (r"\bof mammary gland\b", ""), (r"\bof artery\b", "(artery)"),
    (r"\bof lymphatic vessel\b", "(lymphatic)"),
    (r"\bconventional\b", "conv."), (r"\bclassical\b", "class."),
    (r"\bnon-classical\b", "non-class."), (r"\beffector memory\b", "EM"),
    (r"\bcentral memory\b", "CM"), (r"\bregulatory\b", "reg."),
    (r"\bnatural killer\b", "NK"), (r"\bdendritic cell\b", "DC"),
    (r"\bmacrophage\b", "macro."), (r"\bvascular associated\b", "vasc."),
]


def shorten_cell_type_label(label: str, *, max_len: int = 28) -> str:
    """Return a compact display label for axes/legends (hover keeps the full one)."""
    s = str(label)
    for pat, rep in _STRIP_PATTERNS:
        s = re.sub(pat, rep, s, flags=re.IGNORECASE)
    s = re.sub(r"\s+", " ", s).strip().strip(",").strip()
    if len(s) > max_len:
        s = s[: max_len - 1].rstrip() + "…"
    return s or str(label)


def wrap_label(label: str, width: int = 18) -> str:
    """Insert ``<br>`` at word boundaries so long labels wrap in hover/captions."""
    words = str(label).split()
    lines, cur = [], ""
    for w in words:
        if cur and len(cur) + 1 + len(w) > width:
            lines.append(cur)
            cur = w
        else:
            cur = f"{cur} {w}".strip()
    if cur:
        lines.append(cur)
    return "<br>".join(lines)


def save_color_map(color_map: dict[str, str], out_dir, *,
                   name: str = "cell_type_color_map.tsv") -> Path:
    """Persist the cell-type → colour mapping (with family + short label)."""
    out_dir = Path(out_dir)
    out_dir.mkdir(parents=True, exist_ok=True)
    rows = [{
        "cell_type": ct,
        "family": infer_cell_type_family(ct),
        "short_label": shorten_cell_type_label(ct),
        "color": color,
    } for ct, color in color_map.items()]
    df = pd.DataFrame(rows)
    path = out_dir / name
    df.to_csv(path, sep="\t", index=False)
    return path
