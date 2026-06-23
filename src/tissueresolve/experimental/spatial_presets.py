"""Experimental spatial smoothing presets (opt-in; default behaviour untouched).

A *preset* only sets ``cfg.spatial_solver.lambda_spatial`` to a documented value
and returns provenance metadata. It does NOT add a solver, change soft gating,
touch unresolved mass, or alter bulk behaviour. The package default
(``lambda_spatial = 0.1``) is unchanged: the ``"default"`` preset is a no-op.

Presets
-------
default          λ = 0.1   (package default, NOT experimental)
weak_smoothing   λ = 0.02  (EXPERIMENTAL — less spatial smoothing; intended for
                            boundary / rare-niche analyses where the default
                            oversmooths)
no_smoothing     λ = 0.0   (EXPERIMENTAL — pure NB-MAP per spot; diagnostic control)
"""
from __future__ import annotations

from dataclasses import dataclass

FEATURE_STATUS = "experimental"

# Package default spatial lambda (kept in sync with SpatialSolverConfig.lambda_spatial).
DEFAULT_LAMBDA_SPATIAL = 0.1

SPATIAL_PRESETS: dict[str, dict] = {
    "default": {"lambda_spatial": DEFAULT_LAMBDA_SPATIAL, "experimental": False,
                "description": "package default CAR smoothing"},
    "weak_smoothing": {"lambda_spatial": 0.02, "experimental": True,
                       "description": "experimental: reduced CAR smoothing for "
                                      "boundary / rare-niche analyses"},
    "no_smoothing": {"lambda_spatial": 0.0, "experimental": True,
                     "description": "experimental: no spatial smoothing (NB-MAP control)"},
}


@dataclass
class SpatialPresetInfo:
    preset: str
    lambda_spatial: float
    experimental: bool
    smoothing_used: bool
    is_default: bool
    description: str

    def to_metadata(self) -> dict:
        return {
            "spatial_preset": self.preset,
            "lambda_spatial": self.lambda_spatial,
            "smoothing_used": self.smoothing_used,
            "is_default_spatial_behaviour": self.is_default,
            "spatial_preset_experimental": self.experimental,
            "spatial_preset_description": self.description,
        }


def available_presets() -> list[str]:
    return list(SPATIAL_PRESETS)


def apply_spatial_preset(cfg, preset: str) -> SpatialPresetInfo:
    """Set ``cfg.spatial_solver.lambda_spatial`` from a named preset, in place.

    Returns provenance metadata. ``"default"`` is a no-op that leaves the
    configured (package-default) lambda untouched. Raises on an unknown preset.
    """
    if preset not in SPATIAL_PRESETS:
        raise ValueError(
            f"Unknown spatial preset {preset!r}. "
            f"Available: {', '.join(SPATIAL_PRESETS)}.")
    spec = SPATIAL_PRESETS[preset]
    if preset != "default":
        cfg.spatial_solver.lambda_spatial = float(spec["lambda_spatial"])
    lam = float(cfg.spatial_solver.lambda_spatial)
    return SpatialPresetInfo(
        preset=preset,
        lambda_spatial=lam,
        experimental=bool(spec["experimental"]),
        smoothing_used=lam > 0.0,
        is_default=(preset == "default"),
        description=str(spec["description"]),
    )
