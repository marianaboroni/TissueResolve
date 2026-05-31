"""BisqueRNA bulk deconvolution wrapper (optional)."""
from __future__ import annotations

import pandas as pd

from benchmarks.shared.base import BenchmarkMethod
from benchmarks.shared import environment as _env
from benchmarks.bulk.methods.music_wrapper import _export_bulk_inputs


class BisqueWrapper(BenchmarkMethod):
    name = "Bisque"
    modality = "bulk"
    requires_raw_counts = True
    supports_single_nucleus_reference = True
    external = True

    def is_available(self) -> bool:
        return _env.python_module_available("bisque") or _env.r_package_available("BisqueRNA")

    def install_hint(self) -> str:
        return _env.INSTALL_HINTS["bisque"]

    def _run(self, scenario: dict) -> pd.DataFrame:
        _export_bulk_inputs(scenario, self.name)
        raise NotImplementedError(
            "Bisque is installed but automatic execution is not wired; inputs exported.")
