"""
Tests: package import and public API surface.
"""
import pytest


def test_package_imports_cleanly():
    import tissueresolve
    assert hasattr(tissueresolve, "__version__")
    assert isinstance(tissueresolve.__version__, str)
    assert len(tissueresolve.__version__) > 0


def test_result_classes_importable():
    from tissueresolve.results import (
        BenchmarkResult,
        BulkDeconvResult,
        PairSeparability,
        QCReport,
        ReferenceSignature,
        SeparabilityReport,
        SpatialDeconvResult,
    )
    assert ReferenceSignature is not None
    assert BulkDeconvResult is not None
    assert SpatialDeconvResult is not None
    assert QCReport is not None
    assert PairSeparability is not None
    assert SeparabilityReport is not None
    assert BenchmarkResult is not None


def test_config_importable():
    from tissueresolve.config import (
        BootstrapConfig,
        BulkQCConfig,
        BulkSolverConfig,
        DiscordanceConfig,
        GeneConfig,
        ReferenceConfig,
        SpatialQCConfig,
        SpatialSolverConfig,
        TissueResolveConfig,
    )
    cfg = TissueResolveConfig()
    assert cfg is not None
    assert cfg.reference.genome == "hg38"


def test_utils_importable():
    from tissueresolve.utils import (
        MemoryTracker,
        bh_fdr,
        dense_subset,
        nb_loglik_batch,
        project_simplex_batch,
        safe_log,
        set_random_state,
        setup_logging,
        timer,
    )
    assert callable(bh_fdr)
    assert callable(setup_logging)


def test_subpackage_imports():
    import tissueresolve.bulk
    import tissueresolve.io
    import tissueresolve.plotting
    import tissueresolve.protocol
    import tissueresolve.reference
    import tissueresolve.report
    import tissueresolve.spatial
    import tissueresolve.uncertainty


def test_cli_importable():
    from tissueresolve.cli import cli
    assert cli is not None


def test_all_exports_declared():
    """Verify __all__ is defined and non-empty for core modules."""
    import tissueresolve.results as results_mod
    import tissueresolve.utils as utils_mod
    import tissueresolve.config as config_mod
    assert hasattr(results_mod, "__all__"), "results.py must declare __all__"
    assert hasattr(utils_mod, "__all__"), "utils.py must declare __all__"
    assert hasattr(config_mod, "__all__"), "config.py must declare __all__"
    assert len(results_mod.__all__) > 0
    assert len(utils_mod.__all__) > 0
    assert len(config_mod.__all__) > 0
