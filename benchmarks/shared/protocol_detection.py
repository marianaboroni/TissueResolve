"""
Protocol detection for reference and query, with confidence and conservative
defaults.  Wraps tissueresolve.protocol.detect when available and adds
benchmark-level helpers.  Never guesses silently — unknown stays unknown with a
CAUTION.
"""
from __future__ import annotations

from typing import Optional

import pandas as pd

from . import library_type as _lt

__all__ = ["detect_reference_protocol", "detect_query_protocol",
           "protocol_compatibility"]


def detect_reference_protocol(adata, library_col: Optional[str] = None) -> dict:
    """Detect reference protocol/library type from obs/uns metadata."""
    info = {"protocol": "unknown", "library_type": "unknown",
            "confidence": 0.0, "source": "none", "caution": True}
    try:
        lt = _lt.detect_library_type(adata.obs, library_col)
        info["library_type"] = lt["overall"]
        info["confidence"] = lt["confidence"]
        info["source"] = lt["column"] or "none"
    except Exception:
        pass
    # try the package detector if present
    try:
        from tissueresolve.protocol import detect as pdet
        if hasattr(pdet, "detect_reference_protocol"):
            pkg = pdet.detect_reference_protocol(adata)
            if isinstance(pkg, dict):
                info.update({k: v for k, v in pkg.items() if v})
    except Exception:
        pass
    info["caution"] = info["library_type"] in ("unknown",) or info["confidence"] < 0.5
    return info


def detect_query_protocol(kind: str, meta: Optional[dict] = None) -> dict:
    """Detect query protocol from a hint ('bulk'/'visium'/...) and optional meta."""
    meta = meta or {}
    label = "unknown"
    k = str(kind).lower()
    if "bulk" in k:
        label = "bulk_rnaseq"
    elif "visium" in k or "spatial" in k:
        ffpe = str(meta.get("preservation", "")).lower()
        label = "visium_ffpe" if "ffpe" in ffpe else "visium"
    return {"protocol": label, "confidence": 0.5 if label != "unknown" else 0.0,
            "caution": label == "unknown"}


def protocol_compatibility(ref_protocol: dict, query_protocol: dict) -> dict:
    """Coarse compatibility note between a reference and a query protocol."""
    ref_lib = ref_protocol.get("library_type", "unknown")
    q = query_protocol.get("protocol", "unknown")
    risk = "low"
    notes = []
    if ref_lib == "single_nucleus" and q.startswith("bulk"):
        risk, n = "medium", ("snRNA reference vs bulk (polyA) — intronic/3' "
                             "composition differs; protocol-aware weighting advised.")
        notes.append(n)
    if ref_lib == "unknown" or q == "unknown":
        risk = "unknown"
        notes.append("Protocol unknown for reference and/or query; using "
                     "conservative defaults (CAUTION).")
    if ref_lib == "mixed":
        notes.append("Mixed sc/sn reference — check library-type confounding "
                     "by cell type.")
    return {"risk": risk, "notes": notes,
            "reference_library_type": ref_lib, "query_protocol": q}
