"""TissueResolve bulk methods — flat and hierarchical broad→fine."""
from __future__ import annotations

import warnings

import pandas as pd

from benchmarks.shared.base import BenchmarkMethod


class TissueResolveBulkFlat(BenchmarkMethod):
    name = "TissueResolve_flat"
    modality = "bulk"
    requires_raw_counts = True
    supports_single_nucleus_reference = True
    supports_mixed_sc_sn_reference = True
    supports_hierarchical_reference = False
    external = False

    def _run(self, scenario: dict) -> pd.DataFrame:
        import tissueresolve as tr
        with warnings.catch_warnings():
            warnings.simplefilter("ignore")
            res = tr.deconv_bulk(scenario["bulk"], scenario["reference"],
                                 resolution_mode="flat", n_bootstrap=0)
        return res.deconv.proportions


class TissueResolveBulkAuto(BenchmarkMethod):
    """Improved TissueResolve: solver=auto backbone chosen by gene-masking CV."""
    name = "TissueResolve_auto"
    modality = "bulk"
    requires_raw_counts = True
    supports_single_nucleus_reference = True
    supports_mixed_sc_sn_reference = True
    supports_hierarchical_reference = False
    external = False

    def _run(self, scenario: dict):
        import warnings
        from tissueresolve.solver import AutoSolver
        from benchmarks.shared.io import OUTPUTS_DIR, write_tsv, write_json
        import pandas as pd
        auto = AutoSolver(n_splits=2)
        with warnings.catch_warnings():
            warnings.simplefilter("ignore")
            best, comp, reason = auto.select(scenario["bulk"], scenario["reference"])
            res = best.solve(scenario["bulk"], scenario["reference"])
        out = OUTPUTS_DIR / "bulk"
        write_tsv(comp, out / "solver_comparison.tsv")
        write_json({"selected_solver": best.name, "reason": reason},
                   out / "selected_solver.json")
        return res.proportions


class TissueResolveBulkHierarchical(BenchmarkMethod):
    name = "TissueResolve_hierarchical"
    modality = "bulk"
    requires_raw_counts = True
    supports_single_nucleus_reference = True
    supports_mixed_sc_sn_reference = True
    supports_hierarchical_reference = True
    external = False

    def _run(self, scenario: dict) -> pd.DataFrame:
        import tissueresolve as tr
        mapping = scenario.get("hierarchy_mapping")
        if not mapping:
            raise ValueError(
                "TissueResolve_hierarchical needs a hierarchy_mapping in the scenario")
        with warnings.catch_warnings():
            warnings.simplefilter("ignore")
            res = tr.deconv_bulk(scenario["bulk"], scenario["reference"],
                                 resolution_mode="hierarchical",
                                 hierarchy_mapping=mapping, n_bootstrap=0)
        # combined fine (resolved subtypes + unresolved_<family>); rows sum to 1
        return res.deconv.proportions


class TissueResolveBulkStateAware(BenchmarkMethod):
    """EXPERIMENTAL: broad→cell-type→state hierarchical deconvolution.

    Distinct from ``TissueResolve_hierarchical`` (the standard broad→fine path):
    this routes through the state-aware solver.  On a reference with no state
    labels it runs the two-level broad→cell-type fallback (so its proportions are
    cell types + ``unresolved_<broad>``, comparable to the hierarchical method);
    state-level estimates appear only when state labels exist.
    """
    name = "TissueResolve_state_aware"
    modality = "bulk"
    requires_raw_counts = True
    supports_single_nucleus_reference = True
    supports_mixed_sc_sn_reference = True
    supports_hierarchical_reference = True
    external = False

    def _run(self, scenario: dict) -> pd.DataFrame:
        import tissueresolve as tr
        mapping = scenario.get("hierarchy_mapping")
        if not mapping:
            raise ValueError(
                "TissueResolve_state_aware needs a hierarchy_mapping in the scenario")
        with warnings.catch_warnings():
            warnings.simplefilter("ignore")
            res = tr.deconv_bulk(
                scenario["bulk"], scenario["reference"],
                resolution_mode="hierarchical", hierarchy_mapping=mapping,
                state_aware=True,
                state_to_celltype=scenario.get("state_to_celltype"),
                n_bootstrap=0)
        # state proportions when states exist, else the two-level fallback
        # (cell types + unresolved_<broad>); rows sum to 1 either way
        return res.state_proportions if res.state_proportions is not None \
            else res.cell_type_proportions
