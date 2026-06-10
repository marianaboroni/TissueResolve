"""Experimental partial-confidence-weighted unresolved mass (Phase 2A).

EXPERIMENTAL — not imported by default, not wired into any pipeline. Replaces the
binary hierarchical subtype gate with a continuous, calibrated confidence weight
while conserving mass.  See docs/PARTIAL_CONFIDENCE_GATING_AUDIT.md.
"""
from .config import SoftGatingConfig, SoftGatingResult, FEATURE_STATUS, ALGORITHM_VERSION
from .partial_gating import apply_partial_confidence_gating
from .confidence import (
    compute_reference_confidence_features, add_run_confidence_features,
    ConfidenceModel, fit_confidence_model, calibration_metrics,
)
from .joint_solver import (
    JointHierarchyResult, fit_alternating_soft_hierarchy, fit_joint_soft_hierarchy,
    reconcile_broad_fine, compute_hierarchy_consistency, compute_error_decomposition,
)

__all__ = [
    "SoftGatingConfig", "SoftGatingResult", "FEATURE_STATUS", "ALGORITHM_VERSION",
    "apply_partial_confidence_gating",
    "compute_reference_confidence_features", "add_run_confidence_features",
    "ConfidenceModel", "fit_confidence_model", "calibration_metrics",
    "JointHierarchyResult", "fit_alternating_soft_hierarchy", "fit_joint_soft_hierarchy",
    "reconcile_broad_fine", "compute_hierarchy_consistency", "compute_error_decomposition",
]
