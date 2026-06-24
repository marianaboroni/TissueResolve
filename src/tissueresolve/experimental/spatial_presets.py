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
    "edge_aware_smoothing": {"lambda_spatial": 0.05, "experimental": True,
                             "edge_aware": True,
                             "description": "experimental: in-solver edge-aware smoothing "
                                            "(edge-weighted spot graph; less smoothing across "
                                            "likely boundaries)"},
    "combined_weak_edge_smoothing": {"lambda_spatial": 0.02, "experimental": True,
                                     "edge_aware": True, "min_edge_weight": 0.05,
                                     "combined": True,
                                     "description": "experimental: combined weak (low lambda) + "
                                                    "in-solver edge-aware smoothing"},
    # Experimental state-similarity / Redeconve-inspired presets (opt-in). These
    # leave lambda_spatial at the package default and instead enable the post-fit,
    # within-family, mass-conserving state refinement (no solver change). They are
    # no-ops unless a fine->broad family_map is supplied to deconv_spatial.
    "state_regularized_experimental": {
        "lambda_spatial": DEFAULT_LAMBDA_SPATIAL, "experimental": True,
        "state_regularization": {"enabled": True, "mode": "state_regularized",
                                 "lambda_state": 0.01, "lambda_sparse": 0.001},
        "description": "experimental: within-family state-similarity concentration + "
                       "sparsity refinement (post-fit; not Redeconve)"},
    "sparsity_state_regularized_experimental": {
        "lambda_spatial": DEFAULT_LAMBDA_SPATIAL, "experimental": True,
        "state_regularization": {"enabled": True, "mode": "sparsity",
                                 "lambda_state": 0.0, "lambda_sparse": 0.01},
        "description": "experimental: within-family sparsity-aware refinement "
                       "(reduce inflated effective-N; rare states protected)"},
    "adaptive_resolution_experimental": {
        "lambda_spatial": DEFAULT_LAMBDA_SPATIAL, "experimental": True,
        "state_regularization": {"enabled": True, "mode": "adaptive_resolution"},
        "description": "experimental: diagnostic state grouping / adaptive resolution "
                       "(estimates unchanged; reporting only)"},
    # Experimental in-solver state-regularized solver presets (Option A; opt-in).
    # A separate projected-gradient solver minimises a joint objective; the
    # production NB-CAR solver is untouched. No-op for the state term unless a
    # fine->broad family_map is supplied to deconv_spatial. NOT Redeconve.
    "state_regularized_solver_experimental": {
        "lambda_spatial": DEFAULT_LAMBDA_SPATIAL, "experimental": True,
        "state_regularized_solver": {"enabled": True, "state_penalty": "competition",
                                     "lambda_spatial": 0.02, "lambda_state": 0.01,
                                     "lambda_sparse": 0.001},
        "description": "experimental: in-solver state-regularized solver "
                       "(competition penalty; separate projected-gradient path)"},
    "state_regularized_solver_competition": {
        "lambda_spatial": DEFAULT_LAMBDA_SPATIAL, "experimental": True,
        "state_regularized_solver": {"enabled": True, "state_penalty": "competition",
                                     "lambda_spatial": 0.02, "lambda_state": 0.01,
                                     "lambda_sparse": 0.001},
        "description": "experimental: in-solver state-regularized solver, competition penalty"},
    "state_regularized_solver_laplacian": {
        "lambda_spatial": DEFAULT_LAMBDA_SPATIAL, "experimental": True,
        "state_regularized_solver": {"enabled": True, "state_penalty": "laplacian",
                                     "lambda_spatial": 0.02, "lambda_state": 0.01,
                                     "lambda_sparse": 0.001},
        "description": "experimental: in-solver state-regularized solver, laplacian penalty"},
    "state_regularized_solver_weak": {
        "lambda_spatial": DEFAULT_LAMBDA_SPATIAL, "experimental": True,
        "state_regularized_solver": {"enabled": True, "state_penalty": "competition",
                                     "lambda_spatial": 0.02, "lambda_state": 0.005,
                                     "lambda_sparse": 0.0},
        "description": "experimental: in-solver state-regularized solver, weak competition"},
}


@dataclass
class SpatialPresetInfo:
    preset: str
    lambda_spatial: float
    experimental: bool
    smoothing_used: bool
    is_default: bool
    description: str
    edge_aware: bool = False
    combined: bool = False
    min_edge_weight: float | None = None
    state_regularization: bool = False
    state_regularization_mode: str | None = None
    state_regularized_solver: bool = False
    state_solver_penalty: str | None = None

    def to_metadata(self) -> dict:
        return {
            "spatial_preset": self.preset,
            "lambda_spatial": self.lambda_spatial,
            "smoothing_used": self.smoothing_used,
            "edge_aware_smoothing_used": self.edge_aware,
            "combined_preset": self.combined,
            "min_edge_weight": self.min_edge_weight,
            "state_regularization_used": self.state_regularization,
            "state_regularization_mode": self.state_regularization_mode,
            "state_regularized_solver_used": self.state_regularized_solver,
            "state_solver_penalty": self.state_solver_penalty,
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
    # edge-aware flag (default False; only edge-aware/combined presets turn it on)
    if hasattr(cfg.spatial_solver, "edge_aware"):
        cfg.spatial_solver.edge_aware = bool(spec.get("edge_aware", False))
    # optional min edge weight override (combined preset)
    min_ew = spec.get("min_edge_weight")
    if min_ew is not None and hasattr(cfg.spatial_solver, "edge_aware_min_weight"):
        cfg.spatial_solver.edge_aware_min_weight = float(min_ew)
    # optional state-similarity regularization (experimental presets)
    sr_spec = spec.get("state_regularization")
    if sr_spec is not None and hasattr(cfg, "state_regularization"):
        for key, val in sr_spec.items():
            if hasattr(cfg.state_regularization, key):
                setattr(cfg.state_regularization, key, val)
    elif hasattr(cfg, "state_regularization"):
        # presets without a state_regularization block must leave it disabled
        cfg.state_regularization.enabled = False
    # optional in-solver state-regularized solver (experimental presets)
    srs_spec = spec.get("state_regularized_solver")
    if srs_spec is not None and hasattr(cfg, "state_regularized_solver"):
        for key, val in srs_spec.items():
            if hasattr(cfg.state_regularized_solver, key):
                setattr(cfg.state_regularized_solver, key, val)
    elif hasattr(cfg, "state_regularized_solver"):
        cfg.state_regularized_solver.enabled = False
    lam = float(cfg.spatial_solver.lambda_spatial)
    info = SpatialPresetInfo(
        preset=preset,
        lambda_spatial=lam,
        experimental=bool(spec["experimental"]),
        smoothing_used=lam > 0.0,
        is_default=(preset == "default"),
        description=str(spec["description"]),
    )
    info.edge_aware = bool(spec.get("edge_aware", False))
    info.combined = bool(spec.get("combined", False))
    info.min_edge_weight = float(min_ew) if min_ew is not None else None
    info.state_regularization = bool(sr_spec is not None and sr_spec.get("enabled", False))
    info.state_regularization_mode = sr_spec.get("mode") if sr_spec is not None else None
    info.state_regularized_solver = bool(srs_spec is not None and srs_spec.get("enabled", False))
    info.state_solver_penalty = srs_spec.get("state_penalty") if srs_spec is not None else None
    return info
