"""
Visium H&E histology integration for spatial figures.

Overlays spots / predictions on the tissue image when available, with a clean
coordinate-only fallback (and a recorded warning) when it is not.  Any axis
inversion or coordinate transform is **recorded in metadata**, never applied
silently.

The plotting functions take an explicit image array (or ``None``) and a
coordinate DataFrame, so they are fully unit-testable offline with a toy image.
``load_visium_histology_image`` / ``extract_visium_coordinates`` /
``get_spatial_scale_factors`` pull these out of a Visium AnnData.
"""
from __future__ import annotations

from pathlib import Path
from typing import Any, Optional

import numpy as np
import pandas as pd

from tissueresolve.plotting.export import FigureResult, export_figure, require_plotly

__all__ = [
    "get_spatial_scale_factors",
    "load_visium_histology_image",
    "extract_visium_coordinates",
    "plot_he_with_spots",
    "plot_dominant_cell_type_on_he",
    "plot_abundance_on_he",
    "plot_spot_pies_on_he",
]


# ---------------------------------------------------------------------------
# AnnData extraction
# ---------------------------------------------------------------------------


def get_spatial_scale_factors(adata, library_id: Optional[str] = None) -> dict:
    """Return Visium scale factors from ``adata.uns['spatial']`` (or defaults)."""
    sp = getattr(adata, "uns", {}).get("spatial", {}) if adata is not None else {}
    if not sp:
        return {"tissue_hires_scalef": 1.0, "tissue_lowres_scalef": 1.0,
                "spot_diameter_fullres": 1.0, "library_id": None}
    lib = library_id or next(iter(sp))
    sf = sp.get(lib, {}).get("scalefactors", {})
    sf = dict(sf)
    sf["library_id"] = lib
    sf.setdefault("tissue_hires_scalef", 1.0)
    sf.setdefault("tissue_lowres_scalef", 1.0)
    return sf


def load_visium_histology_image(adata, *, resolution: str = "hires",
                                library_id: Optional[str] = None):
    """Return ``(image_array, scalefactors, warnings)`` — prefer hires, else lowres.

    ``image_array`` is ``None`` when no image is present (caller should fall
    back to coordinate-only plots).
    """
    warnings: list[str] = []
    sp = getattr(adata, "uns", {}).get("spatial", {}) if adata is not None else {}
    sf = get_spatial_scale_factors(adata, library_id)
    if not sp:
        warnings.append("No H&E image found in adata.uns['spatial']; using "
                        "coordinate-only spatial plots.")
        return None, sf, warnings
    lib = sf.get("library_id") or next(iter(sp))
    images = sp.get(lib, {}).get("images", {})
    order = [resolution, "hires", "lowres"]
    for res in order:
        if res in images and images[res] is not None:
            img = np.asarray(images[res])
            if img.dtype.kind == "f" and img.max() <= 1.0:
                img = (img * 255).astype(np.uint8)
            sf["used_resolution"] = res
            if res != resolution:
                warnings.append(f"Requested {resolution!r} image unavailable; "
                                f"using {res!r}.")
            return img, sf, warnings
    warnings.append("No usable H&E image array; using coordinate-only plots.")
    return None, sf, warnings


def extract_visium_coordinates(adata, *, scalefactors: Optional[dict] = None,
                               resolution: str = "hires") -> pd.DataFrame:
    """Build a coordinate table: spot, array_row, array_col, x, y, image_x/y.

    ``image_*`` are pixel coordinates scaled by the chosen resolution's scale
    factor.  No axis flips are applied here — the figure records any inversion.
    """
    sf = scalefactors or get_spatial_scale_factors(adata)
    scalef = sf.get(f"tissue_{resolution}_scalef", 1.0)
    obs = adata.obs
    spatial = (adata.obsm["spatial"] if "spatial" in getattr(adata, "obsm", {})
               else None)
    n = adata.n_obs
    if spatial is not None:
        px, py = np.asarray(spatial)[:, 0], np.asarray(spatial)[:, 1]
    else:
        px = obs.get("array_col", pd.Series(range(n))).to_numpy(dtype=float)
        py = obs.get("array_row", pd.Series(range(n))).to_numpy(dtype=float)
    df = pd.DataFrame({
        "spot": list(adata.obs_names),
        "array_row": obs.get("array_row", pd.Series([np.nan] * n)).to_numpy(),
        "array_col": obs.get("array_col", pd.Series([np.nan] * n)).to_numpy(),
        "x": px, "y": py,
        "image_x": px * scalef, "image_y": py * scalef,
        "tissue_status": obs.get("in_tissue", pd.Series([1] * n)).to_numpy(),
    })
    return df


# ---------------------------------------------------------------------------
# Plotting (explicit image + coords → testable)
# ---------------------------------------------------------------------------

_SP_SUB = "spot-level RNA-derived composition, not cell counts"


def _base_figure(go, image: Optional[np.ndarray]):
    fig = go.Figure()
    inverted = False
    if image is not None:
        fig.add_trace(go.Image(z=image))
        inverted = True  # image origin is top-left → y increases downward
    return fig, inverted


def _coord_comment(inverted: bool, image: Optional[np.ndarray]) -> list[str]:
    return [
        "estimate_type: spot_rna_composition",
        f"has_he_image: {image is not None}",
        f"y_axis_inverted: {inverted} (image origin top-left)",
        "columns: spot,x,y,array_row,array_col,image_x,image_y,tissue_status",
    ]


def plot_he_with_spots(image, coords: pd.DataFrame, output_dir, *,
                       name: str = "he_spots_check",
                       color: Optional[str] = "#888888") -> FigureResult:
    """Visual check: H&E (if available) with spots overlaid in grey."""
    go = require_plotly()
    fig, inverted = _base_figure(go, image)
    x = coords.get("image_x", coords.get("x"))
    y = coords.get("image_y", coords.get("y"))
    fig.add_trace(go.Scatter(x=x, y=y, mode="markers",
                             marker={"color": color, "size": 5, "opacity": 0.6},
                             name="spots", hovertext=coords["spot"].astype(str)))
    if image is None:
        fig.update_yaxes(autorange="reversed")
    fig.update_layout(title={"text": "H&E with spots (visual check)", "x": 0.5},
                      template="plotly_white", height=560,
                      xaxis={"visible": False}, yaxis={"visible": False})
    warns = [] if image is not None else ["No H&E image; spots shown on blank canvas."]
    res = export_figure(fig, output_dir, name, data={"data": coords},
                        caption="H&E + spot positions (visual check). " + _SP_SUB,
                        data_comment=_coord_comment(inverted, image), formats=())
    res.warnings.extend(warns)
    return res


def plot_dominant_cell_type_on_he(image, coords: pd.DataFrame,
                                  proportions: pd.DataFrame, color_map: dict,
                                  output_dir, *,
                                  name: str = "he_dominant_cell_type",
                                  low_conf_threshold: float = 0.5) -> FigureResult:
    """Colour each spot by its dominant cell type over the H&E image."""
    go = require_plotly()
    from tissueresolve.plotting.palette import shorten_cell_type_label

    fig, inverted = _base_figure(go, image)
    dom = proportions.idxmax(axis=1).astype(str)
    domfrac = proportions.max(axis=1).to_numpy()
    x = coords.get("image_x", coords.get("x")).to_numpy()
    y = coords.get("image_y", coords.get("y")).to_numpy()
    for ct in dom.unique():
        mask = (dom == ct).to_numpy()
        conf = domfrac[mask] >= low_conf_threshold
        fig.add_trace(go.Scatter(
            x=x[mask], y=y[mask], mode="markers", name=shorten_cell_type_label(ct),
            marker={"color": color_map.get(ct, "#888888"), "size": 6,
                    "opacity": np.where(conf, 0.9, 0.35).tolist()},
            hovertext=[f"{ct}: {f:.2f}" for f in domfrac[mask]]))
    if image is None:
        fig.update_yaxes(autorange="reversed")
    fig.update_layout(title={"text": "Dominant cell type on H&E", "x": 0.5},
                      template="plotly_white", height=560,
                      xaxis={"visible": False}, yaxis={"visible": False})
    data = coords.copy()
    data["dominant_type"] = dom.to_numpy()
    data["dominant_fraction"] = domfrac
    return export_figure(fig, output_dir, name, data={"data": data},
                         caption="Dominant cell type per spot on H&E. " + _SP_SUB,
                         data_comment=_coord_comment(inverted, image), formats=())


def plot_abundance_on_he(image, coords: pd.DataFrame, abundance: pd.Series,
                         output_dir, *, cell_type: str = "abundance",
                         name: Optional[str] = None) -> FigureResult:
    """Continuous abundance of one cell type over the H&E image."""
    go = require_plotly()
    fig, inverted = _base_figure(go, image)
    x = coords.get("image_x", coords.get("x"))
    y = coords.get("image_y", coords.get("y"))
    fig.add_trace(go.Scatter(
        x=x, y=y, mode="markers",
        marker={"color": abundance.to_numpy(), "colorscale": "Viridis",
                "cmin": 0, "cmax": 1, "size": 6,
                "colorbar": {"title": "composition"}}))
    if image is None:
        fig.update_yaxes(autorange="reversed")
    fig.update_layout(title={"text": f"Abundance: {cell_type}", "x": 0.5},
                      template="plotly_white", height=560,
                      xaxis={"visible": False}, yaxis={"visible": False})
    data = coords.copy()
    data[str(cell_type)] = abundance.to_numpy()
    return export_figure(fig, output_dir, name or f"he_abundance_{cell_type}",
                         data={"data": data},
                         caption=f"{cell_type} abundance on H&E. " + _SP_SUB,
                         data_comment=_coord_comment(inverted, image), formats=())


def plot_spot_pies_on_he(image, coords: pd.DataFrame, proportions: pd.DataFrame,
                         color_map: dict, output_dir, *, max_spots: int = 200,
                         top_n: int = 5, name: str = "he_spot_pies") -> FigureResult:
    """Per-spot pie charts over the H&E image (exploratory; capped + warned)."""
    go = require_plotly()
    warns: list[str] = []
    n = len(proportions)
    idx = np.arange(n)
    if n > max_spots:
        idx = idx[:: int(np.ceil(n / max_spots))]
        warns.append(f"{n} spots > max_spots={max_spots}; showing "
                     f"{len(idx)} sampled spots.")
    # top_n + Other
    order = proportions.mean(0).sort_values(ascending=False).index.tolist()
    keep = order[:top_n]
    props = proportions[keep].copy()
    if len(order) > top_n:
        props["Other"] = proportions[order[top_n:]].sum(axis=1)
    cts = list(props.columns)
    colours = [color_map.get(c, "#d9d9d9") for c in cts]
    x = coords.get("image_x", coords.get("x")).to_numpy()
    y = coords.get("image_y", coords.get("y")).to_numpy()
    xn = (x - x.min()) / (x.max() - x.min() + 1e-9)
    yn = (y - y.min()) / (y.max() - y.min() + 1e-9)
    fig = go.Figure()
    half = 0.02
    for s in idx:
        vals = props.iloc[int(s)].to_numpy(dtype=float)
        if vals.sum() <= 0:
            continue
        fig.add_trace(go.Pie(labels=cts, values=vals, sort=False, textinfo="none",
                             marker={"colors": colours}, showlegend=bool(s == idx[0]),
                             domain={"x": [max(0, xn[s] - half), min(1, xn[s] + half)],
                                     "y": [max(0, yn[s] - half), min(1, yn[s] + half)]}))
    fig.update_layout(title={"text": "Per-spot pies (exploratory)", "x": 0.5},
                      template="plotly_white", height=620)
    long = []
    for s in idx:
        for ct in cts:
            long.append({"spot": str(props.index[int(s)]), "image_x": x[s],
                         "image_y": y[s], "cell_type": ct,
                         "fraction": float(props.iloc[int(s)][ct])})
    res = export_figure(fig, output_dir, name, data={"data": pd.DataFrame(long)},
                        caption="Exploratory per-spot pies on H&E. " + _SP_SUB,
                        data_comment=["estimate_type: spot_rna_composition"], formats=())
    res.warnings.extend(warns)
    return res
