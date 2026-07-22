"""
Environment probing — which optional benchmark dependencies are installed.

Used so the benchmark can run with only internal baselines and gracefully skip
external tools, recording the reason and an install hint for each.
"""
from __future__ import annotations

import importlib.util
import shutil
import subprocess
from pathlib import Path

__all__ = ["python_module_available", "r_package_available", "package_versions",
           "INSTALL_HINTS", "rectangle_env_python", "rectangle_available"]

#: Isolated Python >=3.10 env that holds rectanglepy (see
#: benchmarks/envs/install_rectangle.sh). Rectangle cannot import into the
#: project 3.9 venv, so it is invoked out-of-process — like cell2location.
_RECTANGLE_ENV_PY = (Path(__file__).resolve().parents[1]
                     / "envs" / "rectangle_py311" / "bin" / "python")

INSTALL_HINTS = {
    "rectangle": ("Isolated Python>=3.10 env: bash benchmarks/envs/"
                  "install_rectangle.sh (installs rectanglepy in "
                  "benchmarks/envs/rectangle_py311; run via subprocess)."),
    # bulk
    "music": "R/Bioconductor: BiocManager::install('MuSiC') (via rpy2)",
    "bisque": "pip install BisqueRNA  (or R: install.packages('BisqueRNA'))",
    "dwls": "R: DWLS from github.com/dtsoucas/DWLS (via rpy2)",
    "scdc": "R: SCDC from github.com/meichendong/SCDC (via rpy2)",
    "cibersortx": "Web tool — TissueResolve exports CIBERSORTx-compatible inputs only.",
    "bayesprism": "R: BayesPrism from github.com/Danko-Lab/BayesPrism (export-only here).",
    # spatial
    "rctd": "R: spacexr (RCTD) from github.com/dmcable/spacexr (via rpy2)",
    "cell2location": "pip install cell2location  (needs scvi-tools, GPU recommended)",
    "stereoscope": "pip install scvi-tools  (stereoscope model)",
    "spotlight": "R/Bioconductor: BiocManager::install('SPOTlight') (via rpy2)",
    "tangram": "pip install tangram-sc",
    "destvi": "pip install scvi-tools  (DestVI model)",
    "card": "R: CARD from github.com/YingMa0107/CARD (via rpy2)",
}


def python_module_available(name: str) -> bool:
    try:
        return importlib.util.find_spec(name) is not None
    except (ImportError, ValueError):
        return False


def r_package_available(pkg: str) -> bool:
    """Best-effort: needs rpy2 + the R package.  Never raises."""
    if not python_module_available("rpy2"):
        return False
    if shutil.which("R") is None and shutil.which("Rscript") is None:
        return False
    try:
        from rpy2.robjects.packages import isinstalled  # type: ignore
        return bool(isinstalled(pkg))
    except Exception:
        return False


def rectangle_env_python() -> "Path | None":
    """Path to the isolated Rectangle interpreter, or None if not provisioned."""
    return _RECTANGLE_ENV_PY if _RECTANGLE_ENV_PY.exists() else None


def rectangle_available() -> bool:
    """True iff the isolated env exists and can import rectanglepy. Never raises."""
    py = rectangle_env_python()
    if py is None:
        return False
    try:
        r = subprocess.run(
            [str(py), "-c", "import rectanglepy"],
            capture_output=True, timeout=60)
        return r.returncode == 0
    except (OSError, subprocess.SubprocessError):
        return False


def package_versions() -> dict:
    import platform
    out = {"python": platform.python_version()}
    for mod in ("numpy", "scipy", "pandas", "anndata", "tissueresolve",
                "scvi", "cell2location", "tangram", "rpy2"):
        try:
            import importlib.metadata as m
            out[mod] = m.version(mod)
        except Exception:
            out[mod] = "not installed"
    return out
