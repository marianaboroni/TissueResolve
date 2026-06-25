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


class ExecutedExternalMethod(BenchmarkMethod):
    """An external tool that was **executed locally** (its R/Python runner ran in
    this environment and wrote a predictions TSV + executed metadata).

    Distinct from ImportedMethod (results produced elsewhere): this is marked
    ``executed`` (not imported)."""

    external = True

    def __init__(self, name: str, modality: str, path: Path, runtime: float = 0.0):
        self.name = name
        self.modality = modality
        self._path = Path(path)
        self._runtime = runtime

    def is_available(self) -> bool:
        return self._path.exists()

    def run(self, scenario: dict) -> MethodResult:
        if not self._path.exists():
            return MethodResult(method=self.name, modality=self.modality,
                                status="skipped",
                                skip_reason=f"no prediction at {self._path}")
        df = pd.read_csv(self._path, sep="\t", index_col=0, comment="#")
        s = df.sum(axis=1).replace(0, 1.0)
        df = df.div(s, axis=0)
        return MethodResult(method=self.name, modality=self.modality,
                            status=STATUS_SUCCESS, predictions=df,
                            runtime_s=self._runtime,
                            metadata={"executed_or_exported": "executed",
                                      "external_tool": True,
                                      "source_file": str(self._path)},
                            warnings=["external tool executed locally via its runner"])


def discover_executed_external(modality: str) -> list:
    """Find external tools that executed locally (predictions/ + executed metadata)."""
    import json
    base = OUTPUTS_DIR / modality
    pred_dir = base / "predictions"
    meta_dir = base / "method_metadata"
    if not pred_dir.exists():
        return []
    out = []
    for p in sorted(pred_dir.glob("*.tsv")):
        if p.name.endswith(".data.tsv"):
            continue
        runtime = 0.0
        mp = meta_dir / f"{p.stem}.json"
        executed = True
        if mp.exists():
            try:
                m = json.loads(mp.read_text())
                executed = bool(m.get("executed", True))
                runtime = float(m.get("runtime_seconds", 0) or 0)
            except Exception:
                pass
        if executed:
            out.append(ExecutedExternalMethod(p.stem, modality, p, runtime))
    return out


def discover_imported(modality: str) -> list[ImportedMethod]:
    """Find imported external prediction TSVs for a modality.

    Looks in both ``<modality>/imported/`` (the documented destination) and
    ``<modality>/external/`` (legacy)."""
    found = {}
    for sub in ("imported", "external"):
        d = OUTPUTS_DIR / modality / sub
        if not d.exists():
            continue
        for p in sorted(d.glob("*.tsv")):
            if p.name.endswith(".data.tsv"):
                continue
            found.setdefault(p.stem, ImportedMethod(p.stem, modality, p))
    return list(found.values())
