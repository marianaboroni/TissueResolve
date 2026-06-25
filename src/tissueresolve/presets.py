from __future__ import annotations

from typing import Dict, Any

PRESETS: Dict[str, Dict[str, Any]] = {
    "quick": {
        "plots": "minimal",
        "bootstrap": False,
        "auto_tune": False,
        "spatial": {"lambda": "auto", "n_neighbors": 6},
        # permissive: split subtypes readily, fewer unresolved families
        "hierarchical": {"min_discriminating_genes": 5,
                         "within_family_spillover_threshold": 0.40,
                         "unresolved_threshold": 0.05},
    },
    "standard": {
        "plots": "standard",
        "bootstrap": False,
        "auto_tune": False,
        "spatial": {"lambda": "auto", "n_neighbors": 8},
        "hierarchical": {"min_discriminating_genes": 10,
                         "within_family_spillover_threshold": 0.30,
                         "unresolved_threshold": 0.10},
    },
    "publication": {
        "plots": "publication",
        "bootstrap": True,
        "n_bootstrap": 100,
        "auto_tune": True,
        "spatial": {"lambda": "auto", "n_neighbors": 8},
        # stricter: demand strong within-family evidence before splitting
        "hierarchical": {"min_discriminating_genes": 30,
                         "within_family_spillover_threshold": 0.25,
                         "unresolved_threshold": 0.12},
    },
    "diagnostic": {
        "plots": "all",
        "bootstrap": True,
        "n_bootstrap": 200,
        "auto_tune": True,
        "spatial": {"lambda": "auto", "n_neighbors": 12},
        "hierarchical": {"min_discriminating_genes": 5,
                         "within_family_spillover_threshold": 0.35,
                         "unresolved_threshold": 0.08},
    },
}


def get_preset(name: str) -> Dict[str, Any]:
    return PRESETS.get(name, PRESETS["standard"])  # default
