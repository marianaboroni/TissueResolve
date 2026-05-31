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
