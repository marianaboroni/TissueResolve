"""CIBERSORTx-compatible input export (web tool — never run automatically)."""
from __future__ import annotations

import pandas as pd

from benchmarks.shared.base import BenchmarkMethod, MethodResult, STATUS_EXPORTED
from benchmarks.shared import environment as _env
from benchmarks.shared.io import write_tsv, OUTPUTS_DIR


class CIBERSORTxExport(BenchmarkMethod):
    name = "CIBERSORTx_export"
    modality = "bulk"
    requires_raw_counts = True
    external = True

    def is_available(self) -> bool:
        # Export is always possible; it never executes the web tool.
        return True

    def install_hint(self) -> str:
        return _env.INSTALL_HINTS["cibersortx"]

    def run(self, scenario: dict) -> MethodResult:  # override: always export, never run
        out = OUTPUTS_DIR / "bulk" / "exports" / self.name
        bulk = scenario.get("bulk")
        ref = scenario.get("reference")
        files = []
        if bulk is not None:
            # CIBERSORTx "mixture" file: genes × samples, tab-delimited
            files.append(str(write_tsv(bulk, out / "mixture.tsv")))
        if ref is not None:
            sig = pd.DataFrame(ref.as_R_cpm().T, index=list(ref.gene_names),
                               columns=list(ref.cell_types))
            files.append(str(write_tsv(sig, out / "signature_matrix.tsv")))
        return MethodResult(
            method=self.name, modality=self.modality, status=STATUS_EXPORTED,
            skip_reason="CIBERSORTx runs on the web; inputs exported for manual upload",
            install_hint=self.install_hint(),
            metadata={"exported_files": files})
