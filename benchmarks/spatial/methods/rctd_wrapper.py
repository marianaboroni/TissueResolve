"""RCTD spatial deconvolution wrapper (optional)."""
from __future__ import annotations

import pandas as pd

from benchmarks.shared.base import BenchmarkMethod
from benchmarks.shared import environment as _env
from benchmarks.spatial.methods._export import export_spatial_inputs


class RCTDWrapper(BenchmarkMethod):
    name = "RCTD"
    modality = "spatial"
    requires_raw_counts = True
    supports_single_nucleus_reference = True
    external = True

    def is_available(self) -> bool:
        return (_env.python_module_available("rctd")
                or _env.r_package_available("RCTD"))

    def install_hint(self) -> str:
        return _env.INSTALL_HINTS.get("rctd", "see method documentation")

    def _run(self, scenario: dict) -> pd.DataFrame:
        export_spatial_inputs(scenario, self.name)
        raise NotImplementedError(
            "RCTD is installed but automatic execution is not wired in this "
            "harness; spatial inputs were exported for a manual run.")
