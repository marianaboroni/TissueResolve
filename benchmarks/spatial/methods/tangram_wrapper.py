"""Tangram spatial deconvolution wrapper (optional)."""
from __future__ import annotations

import pandas as pd

from benchmarks.shared.base import BenchmarkMethod
from benchmarks.shared import environment as _env
from benchmarks.spatial.methods._export import export_spatial_inputs


class TangramWrapper(BenchmarkMethod):
    name = "Tangram"
    modality = "spatial"
    requires_raw_counts = True
    supports_single_nucleus_reference = True
    external = True

    def is_available(self) -> bool:
        return (_env.python_module_available("tangram")
                or _env.r_package_available("Tangram"))

    def install_hint(self) -> str:
        return _env.INSTALL_HINTS.get("tangram", "see method documentation")

    def _run(self, scenario: dict) -> pd.DataFrame:
        export_spatial_inputs(scenario, self.name)
        raise NotImplementedError(
            "Tangram is installed but automatic execution is not wired in this "
            "harness; spatial inputs were exported for a manual run.")
