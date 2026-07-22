#!/usr/bin/env bash
# Provision an isolated Python >=3.10 environment for Rectangle (rectanglepy).
#
# Rectangle is a benchmark comparison tool only. It requires Python >=3.10, but
# the TissueResolve project venv (.venv) and the cell2location env are Python
# 3.9, so Rectangle gets its own isolated environment invoked via subprocess —
# the same pattern used for cell2location (benchmarks/envs/c2l_py39).
#
# This script is documented for reproducibility; it is NOT run automatically by
# the test suite or the default benchmark. Re-running it is idempotent.
#
# Verified populated state on this machine (macOS/arm64):
#   Homebrew python@3.11 (3.11.15) -> benchmarks/envs/rectangle_py311
#   rectanglepy 1.5.0 (pulls anndata 0.10.8, pydeseq2 0.4.11, scikit-learn 1.9,
#   osqp, statsmodels, matplotlib) — no TissueResolve code is installed here.
set -euo pipefail

REPO="$(cd "$(dirname "${BASH_SOURCE[0]}")/../.." && pwd)"
ENV_DIR="${REPO}/benchmarks/envs/rectangle_py311"

# 1. Ensure a Python >=3.10 interpreter. Homebrew is the route used here.
PY311="$(command -v python3.11 || true)"
if [[ -z "${PY311}" && -x /opt/homebrew/bin/python3.11 ]]; then
  PY311=/opt/homebrew/bin/python3.11
fi
if [[ -z "${PY311}" ]]; then
  echo "Python 3.11 not found. Install it first, e.g.:" >&2
  echo "  brew install python@3.11        # macOS/Homebrew" >&2
  echo "  # or use pyenv/uv to provide a >=3.10 interpreter" >&2
  exit 1
fi
echo "Using interpreter: ${PY311} ($(${PY311} --version))"

# 2. Isolated venv.
if [[ ! -d "${ENV_DIR}" ]]; then
  "${PY311}" -m venv "${ENV_DIR}"
fi
"${ENV_DIR}/bin/python" -m pip install --upgrade pip >/dev/null

# 3. Rectangle. Pin for reproducibility; bump deliberately.
"${ENV_DIR}/bin/pip" install "rectanglepy==1.5.0"

# 4. Verify.
"${ENV_DIR}/bin/python" - <<'PYEOF'
import rectanglepy as rp
print("rectanglepy", rp.version("rectanglepy"), "OK")
PYEOF
echo "Rectangle env ready: ${ENV_DIR}"
