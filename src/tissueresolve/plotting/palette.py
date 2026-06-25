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

import json
import re
from pathlib import Path
from typing import Optional

import pandas as pd

__all__ = [
    "FAMILY_RAMPS",
    "OTHER_COLOR",
    "UNRESOLVED_COLOR",
    "QUALITATIVE_BASE",
    "infer_cell_type_family",
    "assign_family_palette",
    "shorten_cell_type_label",
    "wrap_label",
    "save_color_map",
    "build_hierarchical_color_map",
    "save_hierarchical_color_map",
    "load_hierarchical_color_map",
    "color_dicts_from_map",
]

# Deterministic qualitative base palette for families with no known ramp.
QUALITATIVE_BASE = [
    "#4e79a7", "#f28e2b", "#59a14f", "#e15759", "#76b7b2",
    "#edc948", "#b07aa1", "#ff9da7", "#9c755f", "#bab0ac",
]

# Map canonical broad-family names (as produced by hierarchy.infer_broad_cell_type_family
# and by explicit annotations) → FAMILY_RAMPS key.
_FAMILY_TO_RAMP = {
    "t/nk": "T/NK", "b/plasma": "B/plasma", "myeloid": "myeloid",
    "endothelial": "endothelial", "epithelial": "epithelial",
    "stromal/fibroblast": "stromal", "stromal": "stromal",
    "fibroblast": "stromal", "mural": "mural", "adipocyte": "adipocyte",
    "other": "other",
}

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


# ---------------------------------------------------------------------------
# Reproducible, family-aware colour map (broad/fine annotation aware)
# ---------------------------------------------------------------------------


def _hex_to_rgb(h: str) -> tuple[int, int, int]:
    h = h.lstrip("#")
    return tuple(int(h[i:i + 2], 16) for i in (0, 2, 4))  # type: ignore[return-value]


def _rgb_to_hex(rgb) -> str:
    return "#{:02x}{:02x}{:02x}".format(*(max(0, min(255, int(round(c)))) for c in rgb))


def _blend(hex_color: str, frac: float, toward: str = "#ffffff") -> str:
    """Blend *hex_color* a *frac* of the way toward *toward* (frac∈[0,1])."""
    a, b = _hex_to_rgb(hex_color), _hex_to_rgb(toward)
    return _rgb_to_hex(tuple(a[i] + (b[i] - a[i]) * frac for i in range(3)))


def _family_shades(base: str, n: int, ramp_key: Optional[str]) -> list[str]:
    """Return *n* visually-related shades for a family's fine subtypes."""
    if n <= 0:
        return []
    ramp = FAMILY_RAMPS.get(ramp_key or "", [])
    if len(ramp) >= n:
        return list(ramp[:n])
    if n == 1:
        return [base]
    # interpolate from the base (dark) toward a light tint
    return [_blend(base, 0.65 * i / (n - 1)) for i in range(n)]


def build_hierarchical_color_map(
    fine_types: list[str],
    mapping: Optional[dict[str, str]] = None,
    *,
    family_order: Optional[list[str]] = None,
    existing: Optional[pd.DataFrame] = None,
) -> pd.DataFrame:
    """Build a deterministic, family-aware colour map for a run.

    Broad families each get a visually distinct base colour; the fine subtypes
    within a family use related shades (lighter tints) of that colour.
    ``Other`` is light grey; ``unresolved_*`` / broad-level labels are dark
    grey.  When a previous *existing* colour map is supplied, colours already
    assigned to a cell type are **preserved** (new types get new colours), so
    figures stay consistent across re-runs.

    Parameters
    ----------
    fine_types:
        Fine cell-type labels appearing in the results (subtype columns).
    mapping:
        ``{fine_cell_type: broad_family}``.  When ``None`` or empty, colours
        are assigned to fine types from the qualitative palette (no
        family-aware shading) and a note is recorded in ``palette_source``.
    family_order:
        Optional explicit broad-family ordering (e.g. by abundance).  Defaults
        to sorted family names.
    existing:
        Optional previously-saved colour-map DataFrame to extend without
        changing existing colours.

    Returns
    -------
    pd.DataFrame
        Columns: ``broad_cell_type, fine_cell_type, color_hex, display_label,
        palette_source, color_role``.  One row per broad family
        (``color_role="broad"``) and one per fine/unresolved/other label.
    """
    fine_types = [str(c) for c in fine_types]
    mapping = {str(k): str(v) for k, v in (mapping or {}).items()}
    have_families = bool(mapping)

    locked: dict[tuple[str, str], str] = {}
    if existing is not None and not existing.empty:
        for _, r in existing.iterrows():
            locked[(str(r.get("broad_cell_type", "")),
                    str(r.get("fine_cell_type", "")))] = str(r["color_hex"])

    # group fine types by family (separating special roles)
    fam_members: dict[str, list[str]] = {}
    specials: list[tuple[str, str, str]] = []  # (label, role, color)
    for ct in fine_types:
        low = ct.strip().lower()
        if low == "other":
            specials.append((ct, "other", OTHER_COLOR))
            continue
        if low.startswith("unresolved"):
            specials.append((ct, "unresolved", UNRESOLVED_COLOR))
            continue
        fam = mapping.get(ct, ct if not have_families else "Other")
        fam_members.setdefault(fam, []).append(ct)

    families = family_order or sorted(fam_members)
    # base colour per family (known ramp first, else qualitative cycle)
    base_for: dict[str, tuple[str, Optional[str]]] = {}
    qi = 0
    for fam in families:
        ramp_key = _FAMILY_TO_RAMP.get(fam.strip().lower())
        if ramp_key and ramp_key != "other":
            base = FAMILY_RAMPS[ramp_key][0]
        else:
            base = QUALITATIVE_BASE[qi % len(QUALITATIVE_BASE)]
            qi += 1
            ramp_key = None
        base_for[fam] = (base, ramp_key)

    rows: list[dict] = []

    def _color(broad: str, fine: str, default: str) -> tuple[str, str]:
        key = (broad, fine)
        if key in locked:
            return locked[key], "reused_from_existing"
        return default, ("family_ramp" if have_families else "qualitative")

    for fam in families:
        base, ramp_key = base_for[fam]
        # broad-family row (used for broad-level plots)
        col, src = _color(fam, "", base if have_families else OTHER_COLOR)
        if have_families:
            rows.append({
                "broad_cell_type": fam, "fine_cell_type": "",
                "color_hex": col, "display_label": shorten_cell_type_label(fam),
                "palette_source": src, "color_role": "broad",
            })
        members = sorted(fam_members[fam], key=str)
        shades = _family_shades(base, len(members), ramp_key)
        for ct, shade in zip(members, shades):
            col, src = _color(fam, ct, shade)
            rows.append({
                "broad_cell_type": fam if have_families else "",
                "fine_cell_type": ct,
                "color_hex": col,
                "display_label": shorten_cell_type_label(ct),
                "palette_source": src if have_families else "qualitative_no_family",
                "color_role": "fine",
            })

    for label, role, default in specials:
        col, src = _color("", label, default)
        rows.append({
            "broad_cell_type": "", "fine_cell_type": label,
            "color_hex": col, "display_label": shorten_cell_type_label(label),
            "palette_source": "reused_from_existing" if (("", label) in locked) else "neutral_grey",
            "color_role": role,
        })

    return pd.DataFrame(
        rows,
        columns=["broad_cell_type", "fine_cell_type", "color_hex",
                 "display_label", "palette_source", "color_role"],
    )


def save_hierarchical_color_map(
    color_map: pd.DataFrame, out_dir, *,
    tsv_name: str = "cell_type_color_map.tsv",
    json_name: str = "color_map.json",
) -> dict[str, Path]:
    """Persist a hierarchical colour map as both TSV and JSON.

    The JSON form is a flat ``{label: color_hex}`` dict (fine labels preferred,
    falling back to the broad family name) for quick programmatic reuse.
    """
    out_dir = Path(out_dir)
    out_dir.mkdir(parents=True, exist_ok=True)
    tsv_path = out_dir / tsv_name
    color_map.to_csv(tsv_path, sep="\t", index=False)

    flat: dict[str, str] = {}
    for _, r in color_map.iterrows():
        label = str(r["fine_cell_type"]) or str(r["broad_cell_type"])
        if label:
            flat.setdefault(label, str(r["color_hex"]))
    json_path = out_dir / json_name
    json_path.write_text(json.dumps(flat, indent=2), encoding="utf-8")
    return {"tsv": tsv_path, "json": json_path}


def load_hierarchical_color_map(path) -> pd.DataFrame:
    """Load a colour map saved by :func:`save_hierarchical_color_map` (TSV)."""
    return pd.read_csv(Path(path), sep="\t").fillna("")


def color_dicts_from_map(color_map: pd.DataFrame) -> tuple[dict[str, str], dict[str, str]]:
    """Split a hierarchical colour-map frame into ``(fine→hex, broad→hex)`` dicts.

    Convenience for plotting code that needs to colour either fine cell types or
    broad families from the single deterministic map.  Fine/special rows
    (``color_role`` in ``fine/other/unresolved``) populate the fine dict keyed by
    ``fine_cell_type``; broad rows populate the broad dict keyed by
    ``broad_cell_type``.
    """
    fine: dict[str, str] = {}
    broad: dict[str, str] = {}
    for _, r in color_map.iterrows():
        role = str(r.get("color_role", ""))
        hexc = str(r["color_hex"])
        ft = str(r.get("fine_cell_type", "") or "")
        bt = str(r.get("broad_cell_type", "") or "")
        if role == "broad":
            if bt:
                broad[bt] = hexc
        else:
            if ft:
                fine[ft] = hexc
    return fine, broad
