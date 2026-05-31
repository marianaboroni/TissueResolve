from __future__ import annotations

from typing import Dict, Any

PRESETS: Dict[str, Dict[str, Any]] = {
    "quick": {
        "plots": "minimal",
        "bootstrap": False,
        "auto_tune": False,
        "spatial": {"lambda": "auto", "n_neighbors": 6},
    },
    "standard": {
        "plots": "standard",
        "bootstrap": False,
        "auto_tune": False,
        "spatial": {"lambda": "auto", "n_neighbors": 8},
    },
    "publication": {
        "plots": "publication",
        "bootstrap": True,
        "n_bootstrap": 100,
        "auto_tune": True,
        "spatial": {"lambda": "auto", "n_neighbors": 8},
    },
    "diagnostic": {
        "plots": "all",
        "bootstrap": True,
        "n_bootstrap": 200,
        "auto_tune": True,
        "spatial": {"lambda": "auto", "n_neighbors": 12},
    },
}


def get_preset(name: str) -> Dict[str, Any]:
    return PRESETS.get(name, PRESETS["standard"])  # default
