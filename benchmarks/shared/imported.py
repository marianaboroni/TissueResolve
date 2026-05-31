"""
Imported external-method results.

External tools (RCTD, cell2location, MuSiC, …) that cannot run automatically in
this environment can still be *fairly compared*: run them yourself, then drop
their predictions into ``benchmarks/outputs/<modality>/external/<Method>.tsv``
(obs × cell_type) via ``import_external_results.py``.  The benchmark then treats
them as **executed (imported)** — clearly distinguished from export-only tools.
"""
from __future__ import annotations

from pathlib import Path

import pandas as pd

from benchmarks.shared.base import BenchmarkMethod, MethodResult, STATUS_SUCCESS
from benchmarks.shared.io import OUTPUTS_DIR


class ImportedMethod(BenchmarkMethod):
    """Wrap a TSV of predictions produced by an external tool run elsewhere."""

    external = True

    def __init__(self, name: str, modality: str, path: Path):
        self.name = name
        self.modality = modality
        self._path = Path(path)

    def is_available(self) -> bool:
        return self._path.exists()

    def install_hint(self) -> str:
        return f"import results: python benchmarks/import_external_results.py --method {self.name} ..."

    def run(self, scenario: dict) -> MethodResult:
        if not self._path.exists():
            return MethodResult(method=self.name, modality=self.modality,
                                status="skipped",
                                skip_reason=f"no imported file at {self._path}",
                                install_hint=self.install_hint())
        df = pd.read_csv(self._path, sep="\t", index_col=0, comment="#")
        # renormalise rows to sum to 1 defensively
        s = df.sum(axis=1).replace(0, 1.0)
        df = df.div(s, axis=0)
        return MethodResult(method=self.name, modality=self.modality,
                            status=STATUS_SUCCESS, predictions=df, runtime_s=0.0,
                            metadata={"executed_or_exported": "executed_imported",
                                      "source_file": str(self._path)},
                            warnings=["results imported from an external run "
                                      "(not executed inside TissueResolve)"])


def discover_imported(modality: str) -> list[ImportedMethod]:
    """Find imported external prediction TSVs for a modality."""
    ext_dir = OUTPUTS_DIR / modality / "external"
    if not ext_dir.exists():
        return []
    return [ImportedMethod(p.stem, modality, p) for p in sorted(ext_dir.glob("*.tsv"))
            if not p.name.endswith(".data.tsv")]
