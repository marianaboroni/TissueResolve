"""Config + result dataclasses for experimental partial-confidence gating.

EXPERIMENTAL (Phase 2A). Not imported by default; not wired into any pipeline.
"""
from __future__ import annotations

from dataclasses import dataclass, field
from typing import Any, Optional

import pandas as pd

FEATURE_STATUS = "experimental"
ALGORITHM_VERSION = "soft_gating-0.1.0"


@dataclass
class SoftGatingConfig:
    """Settings for partial-confidence-weighted unresolved mass."""
    min_confidence: float = 0.0          # floor on calibrated confidence
    max_confidence: float = 1.0          # ceiling
    mass_tol: float = 1e-6               # allowed family/total mass drift
    clip_negative: bool = True
    feature_status: str = FEATURE_STATUS
    algorithm_version: str = ALGORITHM_VERSION


@dataclass
class SoftGatingResult:
    """Output of :func:`apply_partial_confidence_gating`."""
    raw_fine_estimates: pd.DataFrame            # samples × subtypes (pre-gating absolute)
    confidence_by_subtype: pd.DataFrame         # samples × subtypes, in [0,1]
    resolved_fine_estimates: pd.DataFrame       # samples × subtypes (after weighting)
    unresolved_by_family: pd.DataFrame          # samples × unresolved_<family>
    mass_conservation_error: float              # max |row_total_after - row_total_before|
    confidence_features: Optional[pd.DataFrame] = None
    warnings: list = field(default_factory=list)
    metadata: dict[str, Any] = field(default_factory=dict)

    @property
    def combined(self) -> pd.DataFrame:
        """Resolved subtypes + ``unresolved_<family>`` columns (rows preserve mass)."""
        return pd.concat([self.resolved_fine_estimates, self.unresolved_by_family], axis=1)
