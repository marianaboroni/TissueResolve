from __future__ import annotations

from typing import Dict, Any

import anndata as ad


def detect_protocol_from_anndata(adata: ad.AnnData) -> Dict[str, Any]:
    # Heuristic checks using uns/obs/var
    out: Dict[str, Any] = {"protocol": "unknown", "confidence": 0.0, "notes": []}
    platform = None
    if "platform" in adata.uns:
        platform = str(adata.uns.get("platform"))
    if platform is None and "platform" in adata.obs.columns:
        platform = str(adata.obs["platform"].iloc[0])
    if platform:
        p = platform.lower()
        if "10x" in p or "10x" in platform:
            out["protocol"] = "scrna_10x_3prime"
            out["confidence"] = 0.6
            out["notes"].append(f"platform={platform}")
    # library type
    lib = adata.uns.get("library_type") or adata.obs.columns[0:1]
    if isinstance(lib, str) and "smart" in lib.lower():
        out["protocol"] = "smartseq"
        out["confidence"] = max(out["confidence"], 0.5)
    # feature types
    ft = adata.var.get("feature_types") if hasattr(adata, "var") else None
    if ft is not None:
        out["notes"].append("var.feature_types present")
    return out


def detect_reference_protocol(adata: ad.AnnData) -> Dict[str, Any]:
    return detect_protocol_from_anndata(adata)


def detect_query_protocol(input_path: str) -> Dict[str, Any]:
    # Best-effort: return unknown for now.
    return {"protocol": "unknown", "confidence": 0.0, "notes": [f"path={input_path}"]}


def protocol_confidence_score(det: Dict[str, Any]) -> float:
    return float(det.get("confidence", 0.0))


def protocol_detection_summary(det: Dict[str, Any]) -> str:
    return f"protocol={det.get('protocol')} confidence={det.get('confidence')} notes={det.get('notes')}"
