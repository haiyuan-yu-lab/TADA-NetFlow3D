#!/usr/bin/env bash
set -euo pipefail
# Each background worker uses one numerical-library thread.
export OMP_NUM_THREADS=1
export MKL_NUM_THREADS=1
export OPENBLAS_NUM_THREADS=1
export NUMEXPR_NUM_THREADS=1

repo_dir="$(cd "$(dirname "${BASH_SOURCE[0]}")/.." && pwd)"
python3 "${repo_dir}/src/netflow3d/run_netflow3d.py" "$@"
