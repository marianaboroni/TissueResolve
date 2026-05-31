"""
Benchmark method interface.

Every benchmarked method (internal or external) implements :class:`BenchmarkMethod`.
The contract is deliberately tolerant: a method that is not installed, or that
cannot run automatically, returns a structured ``MethodResult`` with a status
other than ``"success"`` instead of raising — so one missing tool never stops
the whole benchmark.
"""
from __future__ import annotations

import time
import traceback
from dataclasses import dataclass, field
from typing import Any, Optional

import pandas as pd

STATUS_SUCCESS = "success"
STATUS_SKIPPED = "skipped"          # tool not installed / not applicable
STATUS_FAILED = "failed"            # tool installed but errored
STATUS_EXPORTED = "exported_not_run"  # inputs written for manual/web execution


@dataclass
class MethodResult:
    """Outcome of running one method on one scenario."""

    method: str
    modality: str
    status: str
    predictions: Optional[pd.DataFrame] = None   # obs × cell_type proportions
    runtime_s: float = 0.0
    peak_memory_mb: Optional[float] = None
    warnings: list[str] = field(default_factory=list)
    metadata: dict[str, Any] = field(default_factory=dict)
    skip_reason: Optional[str] = None
    install_hint: Optional[str] = None

    def to_row(self) -> dict[str, Any]:
        return {
            "method": self.method,
            "modality": self.modality,
            "status": self.status,
            "runtime_s": round(self.runtime_s, 3),
            "peak_memory_mb": self.peak_memory_mb,
            "n_warnings": len(self.warnings),
            "skip_reason": self.skip_reason or "",
        }


class BenchmarkMethod:
    """Base class for all benchmarked deconvolution methods.

    Subclasses set the capability attributes and implement :meth:`run`.
    Capability flags let the runner decide compatibility (normalization,
    protocol, reference type) *before* attempting a run, and record the
    decision transparently.
    """

    name: str = "base"
    modality: str = "bulk"  # "bulk" | "spatial"
    supported_input_types: tuple = ("counts",)
    supported_normalizations: tuple = ("counts",)
    supported_protocols: tuple = ("any",)
    requires_raw_counts: bool = True
    requires_normalized_input: bool = False
    requires_single_cell_reference: bool = True
    supports_single_nucleus_reference: bool = True
    supports_mixed_sc_sn_reference: bool = True
    supports_batch_covariates: bool = False
    supports_hierarchical_reference: bool = False
    external: bool = False

    # --- availability ----------------------------------------------------
    def is_available(self) -> bool:
        """Whether the method can run in this environment.  Internal → True."""
        return True

    def install_hint(self) -> str:
        return ""

    # --- compatibility ----------------------------------------------------
    def check_input_compatibility(self, scenario: dict) -> list[str]:
        """Return a list of compatibility warnings (empty = fully compatible)."""
        warns: list[str] = []
        norm = scenario.get("normalization_status", "unknown")
        if self.requires_raw_counts and norm not in ("counts", "raw", "unknown"):
            warns.append(
                f"{self.name} expects raw counts but input looks '{norm}'.")
        ref_kind = scenario.get("reference_library_type", "unknown")
        if ref_kind == "single_nucleus" and not self.supports_single_nucleus_reference:
            warns.append(f"{self.name} does not officially support snRNA references.")
        if ref_kind == "mixed" and not self.supports_mixed_sc_sn_reference:
            warns.append(f"{self.name} does not officially support mixed sc/sn references.")
        return warns

    # --- execution --------------------------------------------------------
    def prepare_inputs(self, scenario: dict) -> dict:
        """Hook to transform the scenario inputs for this method.  Default: pass-through."""
        return scenario

    def _run(self, scenario: dict) -> pd.DataFrame:
        """Actual deconvolution; subclasses implement.  Return obs × cell_type."""
        raise NotImplementedError

    def run(self, scenario: dict) -> MethodResult:
        """Run the method, capturing status/runtime/warnings.  Never raises."""
        if not self.is_available():
            return MethodResult(
                method=self.name, modality=self.modality, status=STATUS_SKIPPED,
                skip_reason="not installed", install_hint=self.install_hint())
        warns = self.check_input_compatibility(scenario)
        t0 = time.perf_counter()
        try:
            prepared = self.prepare_inputs(scenario)
            preds = self._run(prepared)
            return MethodResult(
                method=self.name, modality=self.modality, status=STATUS_SUCCESS,
                predictions=preds, runtime_s=time.perf_counter() - t0,
                warnings=warns,
                metadata={"n_cell_types": int(preds.shape[1]),
                          "n_obs": int(preds.shape[0])})
        except NotImplementedError:
            return MethodResult(
                method=self.name, modality=self.modality, status=STATUS_EXPORTED,
                runtime_s=time.perf_counter() - t0, warnings=warns,
                skip_reason="requires manual/external execution",
                install_hint=self.install_hint())
        except Exception as exc:  # noqa: BLE001 — one method must not stop the suite
            return MethodResult(
                method=self.name, modality=self.modality, status=STATUS_FAILED,
                runtime_s=time.perf_counter() - t0,
                warnings=warns + [f"{type(exc).__name__}: {exc}"],
                skip_reason=traceback.format_exc(limit=3))

    def summarize_warnings(self, result: MethodResult) -> str:
        return f"{result.method}: {len(result.warnings)} warning(s)"
