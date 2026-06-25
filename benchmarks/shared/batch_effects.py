"""
Back-compat shim. Batch-effect diagnostics moved into the core package at
``tissueresolve.diagnostics.batch_effects`` so the shipped library does not
depend on this development benchmark tree (see docs/V0_1_PRUNING_PLAN.md).

This module re-exports the core implementation; new code should import from
``tissueresolve.diagnostics.batch_effects`` directly.
"""
from __future__ import annotations

from tissueresolve.diagnostics.batch_effects import (  # noqa: F401
    __all__,
    BATCH_COL_CANDIDATES,
    detect_batch_columns,
    compute_celltype_batch_confounding,
    compute_marker_batch_stability,
    compute_batch_mixing_score,
    batch_aware_reference_summary,
)
