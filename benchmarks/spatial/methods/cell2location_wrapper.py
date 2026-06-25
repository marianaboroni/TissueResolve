"""cell2location spatial deconvolution wrapper (optional)."""
from __future__ import annotations

import pandas as pd

from benchmarks.shared.base import BenchmarkMethod
from benchmarks.shared import environment as _env
from benchmarks.spatial.methods._export import export_spatial_inputs


class cell2locationWrapper(BenchmarkMethod):
    name = "cell2location"
    modality = "spatial"
    requires_raw_counts = True
    supports_single_nucleus_reference = True
    external = True

    def is_available(self) -> bool:
        return (_env.python_module_available("cell2location")
                or _env.r_package_available("cell2location"))

    def install_hint(self) -> str:
        return _env.INSTALL_HINTS.get("cell2location", "see method documentation")

    def _run(self, scenario: dict) -> pd.DataFrame:
        export_spatial_inputs(scenario, self.name)
        raise NotImplementedError(
            "cell2location is installed but automatic execution is not wired in this "
            "harness; spatial inputs were exported for a manual run.")
