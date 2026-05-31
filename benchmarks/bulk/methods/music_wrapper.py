"""MuSiC bulk deconvolution wrapper (optional, via rpy2 + Bioconductor)."""
from __future__ import annotations

import pandas as pd

from benchmarks.shared.base import BenchmarkMethod
from benchmarks.shared import environment as _env


class MuSiCWrapper(BenchmarkMethod):
    name = "MuSiC"
    modality = "bulk"
    requires_raw_counts = True
    supports_single_nucleus_reference = True
    supports_mixed_sc_sn_reference = False
    external = True

    def is_available(self) -> bool:
        return _env.r_package_available("MuSiC")

    def install_hint(self) -> str:
        return _env.INSTALL_HINTS["music"]

    def _run(self, scenario: dict) -> pd.DataFrame:
        # Available but automatic rpy2 wiring is environment-specific; export
        # inputs and document manual steps rather than guessing the bridge.
        _export_bulk_inputs(scenario, self.name)
        raise NotImplementedError(
            "MuSiC is installed but automatic rpy2 execution is not wired in this "
            "harness; reference + bulk inputs were exported for a manual run.")


def _export_bulk_inputs(scenario: dict, method: str) -> None:
    from benchmarks.shared.io import write_tsv, OUTPUTS_DIR
    out = OUTPUTS_DIR / "bulk" / "exports" / method
    bulk = scenario.get("bulk")
    if bulk is not None:
        write_tsv(bulk, out / "bulk_counts.tsv")
    ref = scenario.get("reference")
    if ref is not None:
        import pandas as pd
        sig = pd.DataFrame(ref.as_R_cpm().T, index=list(ref.gene_names),
                           columns=list(ref.cell_types))
        write_tsv(sig, out / "reference_signature_cpm.tsv")
