"""DWLS bulk deconvolution wrapper (optional, via rpy2)."""
from __future__ import annotations

import pandas as pd

from benchmarks.shared.base import BenchmarkMethod
from benchmarks.shared import environment as _env
from benchmarks.bulk.methods.music_wrapper import _export_bulk_inputs


class DWLSWrapper(BenchmarkMethod):
    name = "DWLS"
    modality = "bulk"
    requires_raw_counts = True
    external = True

    def is_available(self) -> bool:
        return _env.r_package_available("DWLS")

    def install_hint(self) -> str:
        return _env.INSTALL_HINTS["dwls"]

    def _run(self, scenario: dict) -> pd.DataFrame:
        _export_bulk_inputs(scenario, self.name)
        raise NotImplementedError(
            "DWLS is installed but automatic rpy2 execution is not wired; inputs exported.")
