#!/usr/bin/env bash
# Attempt to install external deconvolution tools, logging each attempt.
# Safe to re-run; never touches the main TissueResolve .venv.
# Usage:  bash benchmarks/envs/install_external_tools.sh
set -u
LOGDIR="$(dirname "$0")/../logs"
mkdir -p "$LOGDIR"

r_install () {  # name  R-expression
  local name="$1"; local expr="$2"
  echo ">> installing $name ..."
  if command -v Rscript >/dev/null 2>&1; then
    Rscript -e "$expr" >"$LOGDIR/install_${name}.log" 2>&1 \
      && echo "   $name: ok" || echo "   $name: FAILED (see logs/install_${name}.log)"
  else
    echo "Rscript not found" >"$LOGDIR/install_${name}.log"
    echo "   $name: skipped (no R)"
  fi
}

r_install music      'if(!requireNamespace("BiocManager",quietly=TRUE)) install.packages("BiocManager",repos="https://cloud.r-project.org"); BiocManager::install("MuSiC", ask=FALSE, update=FALSE)'
r_install bayesprism 'if(!requireNamespace("devtools",quietly=TRUE)) install.packages("devtools",repos="https://cloud.r-project.org"); devtools::install_github("Danko-Lab/BayesPrism/BayesPrism")'
r_install bisque     'install.packages("BisqueRNA", repos="https://cloud.r-project.org")'
r_install card       'if(!requireNamespace("devtools",quietly=TRUE)) install.packages("devtools",repos="https://cloud.r-project.org"); devtools::install_github("YingMa0107/CARD")'
r_install spotlight  'if(!requireNamespace("BiocManager",quietly=TRUE)) install.packages("BiocManager",repos="https://cloud.r-project.org"); BiocManager::install("SPOTlight", ask=FALSE, update=FALSE)'

# Python spatial tool (separate env recommended; pip into current env as a fallback)
echo ">> installing cell2location (pip) ..."
python -m pip install cell2location >"$LOGDIR/install_cell2location.log" 2>&1 \
  && echo "   cell2location: ok" || echo "   cell2location: FAILED (see logs/install_cell2location.log)"
