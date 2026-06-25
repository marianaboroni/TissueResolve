"""TissueResolve spatial methods — flat and hierarchical broad→fine."""
from __future__ import annotations

import warnings

import pandas as pd

from benchmarks.shared.base import BenchmarkMethod


def _run_spatial(scenario, resolution_mode, mapping=None):
    import tissueresolve as tr
    with warnings.catch_warnings():
        warnings.simplefilter("ignore")
        res = tr.deconv_spatial(
            scenario["Y"], scenario["reference"], scenario["array_row"],
            scenario["array_col"], scenario["lib_sizes"], scenario["gene_names"],
            spot_ids=scenario.get("spot_ids"),
            resolution_mode=resolution_mode, hierarchy_mapping=mapping,
            run_neighbourhood=False)
    return res.deconv.proportions


class TissueResolveSpatialFlat(BenchmarkMethod):
    name = "TissueResolve_flat"
    modality = "spatial"
    requires_raw_counts = True
    supports_hierarchical_reference = False
    external = False

    def _run(self, scenario: dict) -> pd.DataFrame:
        return _run_spatial(scenario, "flat")


class TissueResolveSpatialHierarchical(BenchmarkMethod):
    name = "TissueResolve_hierarchical"
    modality = "spatial"
    requires_raw_counts = True
    supports_hierarchical_reference = True
    external = False

    def _run(self, scenario: dict) -> pd.DataFrame:
        mapping = scenario.get("hierarchy_mapping")
        if not mapping:
            raise ValueError(
                "TissueResolve_hierarchical needs a hierarchy_mapping in the scenario")
        return _run_spatial(scenario, "hierarchical", mapping)
