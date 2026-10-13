#!/usr/bin/env bash
set -euo pipefail
# Each background worker uses one numerical-library thread.
export OMP_NUM_THREADS=1
export MKL_NUM_THREADS=1
export OPENBLAS_NUM_THREADS=1
export NUMEXPR_NUM_THREADS=1

repo_dir="$(cd "$(dirname "${BASH_SOURCE[0]}")/.." && pwd)"
# Usage: bash workflows/run_netflow3d.sh config/example_config.yaml
#    or: bash workflows/run_netflow3d.sh --job-name NAME --gene-score FILE ...
if [[ $# -eq 1 && "$1" != -* ]]; then
    python3 "${repo_dir}/scripts/run_netflow_from_config.py" "$1"
else
    python3 "${repo_dir}/src/netflow3d/run_netflow3d.py" "$@"
fi
