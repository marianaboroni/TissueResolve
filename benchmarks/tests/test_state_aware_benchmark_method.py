"""The state-aware bulk method is registered as a distinct, non-external row."""
from __future__ import annotations


def test_state_aware_registered_and_distinct():
    from benchmarks.shared.method_registry import bulk_methods
    methods = {m.name: m for m in bulk_methods(include_external=False)}
    assert "TissueResolve_state_aware" in methods
    assert "TissueResolve_hierarchical" in methods           # standard kept
    assert "TissueResolve_auto" in methods
    # distinct rows (no name collision)
    names = [m.name for m in bulk_methods(include_external=False)]
    assert len(set(names)) == len(names)


def test_state_aware_method_is_internal_and_hierarchical():
    from benchmarks.bulk.methods.tissueresolve import TissueResolveBulkStateAware
    m = TissueResolveBulkStateAware()
    assert m.external is False
    assert m.modality == "bulk"
    assert m.supports_hierarchical_reference is True
