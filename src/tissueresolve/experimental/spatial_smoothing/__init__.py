"""Experimental level-specific spatial smoothing (Task D).

EXPERIMENTAL — not imported by default, not wired into any pipeline. Applies
DIFFERENT CAR smoothing strengths to the broad (λ_broad) and fine (λ_fine) levels
by composing the EXISTING spatial solver at two λ values (no new model; no NB-CAR
VI). Hypothesis: broad families tolerate stronger smoothing than fine subtypes.
"""
from .level_specific import (
    fit_level_specific_spatial, LevelSpecificSpatialResult, FEATURE_STATUS, ALGORITHM_VERSION,
)

__all__ = ["fit_level_specific_spatial", "LevelSpecificSpatialResult",
           "FEATURE_STATUS", "ALGORITHM_VERSION"]
