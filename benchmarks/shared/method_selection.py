"""
Method-selection / compatibility advisor.

Given a scenario (normalization, protocol, reference library type) and a list of
methods, decide which can run, which need conversion, and which should be
skipped — recording the rationale.  This drives the protocol/library-aware
comparison tables.
"""
from __future__ import annotations

import pandas as pd

from . import normalization as _norm

__all__ = ["build_compatibility_table"]


def build_compatibility_table(methods, scenario: dict) -> pd.DataFrame:
    """One row per method: availability, normalization fit, compatibility notes."""
    status = scenario.get("normalization_status", "unknown")
    rows = []
    for m in methods:
        avail = m.is_available()
        norm_warns = _norm.validate_method_normalization_requirements(m, status)
        compat = m.check_input_compatibility(scenario)
        rows.append({
            "method": m.name,
            "modality": m.modality,
            "available": avail,
            "external": getattr(m, "external", False),
            "requires_raw_counts": m.requires_raw_counts,
            "supports_snrna": m.supports_single_nucleus_reference,
            "supports_mixed": m.supports_mixed_sc_sn_reference,
            "supports_hierarchical": m.supports_hierarchical_reference,
            "normalization_fit": "ok" if not norm_warns else "; ".join(norm_warns),
            "compatibility_notes": "; ".join(compat) if compat else "ok",
            "install_hint": m.install_hint() if not avail else "",
        })
    return pd.DataFrame(rows).set_index("method")
