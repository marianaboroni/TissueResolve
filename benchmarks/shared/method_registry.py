"""Registry of benchmark methods per modality (internal + optional external)."""
from __future__ import annotations

from benchmarks.shared.base import BenchmarkMethod

__all__ = ["bulk_methods", "spatial_methods", "all_methods"]


def bulk_methods(include_external: bool = True) -> list[BenchmarkMethod]:
    from benchmarks.bulk.methods.nnls_baseline import NNLSBaseline
    from benchmarks.bulk.methods.extra_baselines import (
        WeightedNNLSBaseline, MarkerOnlyNNLSBaseline,
    )
    from benchmarks.bulk.methods.tissueresolve import (
        TissueResolveBulkFlat, TissueResolveBulkHierarchical,
    )
    methods: list[BenchmarkMethod] = [
        TissueResolveBulkHierarchical(),  # recommended/default first
        TissueResolveBulkFlat(),
        NNLSBaseline(),
        WeightedNNLSBaseline(),
        MarkerOnlyNNLSBaseline(),
    ]
    if include_external:
        from benchmarks.bulk.methods.music_wrapper import MuSiCWrapper
        from benchmarks.bulk.methods.bisque_wrapper import BisqueWrapper
        from benchmarks.bulk.methods.dwls_wrapper import DWLSWrapper
        from benchmarks.bulk.methods.cibersortx_export import CIBERSORTxExport
        methods += [MuSiCWrapper(), BisqueWrapper(), DWLSWrapper(),
                    CIBERSORTxExport()]
    return methods


def spatial_methods(include_external: bool = True) -> list[BenchmarkMethod]:
    from benchmarks.spatial.methods.nnls_spot_baseline import NNLSSpotBaseline
    from benchmarks.spatial.methods.tissueresolve import (
        TissueResolveSpatialFlat, TissueResolveSpatialHierarchical,
    )
    methods: list[BenchmarkMethod] = [
        TissueResolveSpatialHierarchical(),
        TissueResolveSpatialFlat(),
        NNLSSpotBaseline(),
    ]
    if include_external:
        from benchmarks.spatial.methods.rctd_wrapper import RCTDWrapper
        from benchmarks.spatial.methods.cell2location_wrapper import cell2locationWrapper
        from benchmarks.spatial.methods.stereoscope_wrapper import stereoscopeWrapper
        from benchmarks.spatial.methods.spotlight_wrapper import SPOTlightWrapper
        from benchmarks.spatial.methods.tangram_wrapper import TangramWrapper
        methods += [RCTDWrapper(), cell2locationWrapper(), stereoscopeWrapper(),
                    SPOTlightWrapper(), TangramWrapper()]
    return methods


def all_methods(modality: str, include_external: bool = True):
    return (bulk_methods(include_external) if modality == "bulk"
            else spatial_methods(include_external))
