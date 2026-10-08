#!/usr/bin/env bash
set -euo pipefail
repo_dir="$(cd "$(dirname "${BASH_SOURCE[0]}")/.." && pwd)"
# Run both complete TADA cohort examples. NetFlow3D example data are not bundled yet.
bash "${repo_dir}/workflows/run_tada.sh" "${repo_dir}/config/example_config.yaml"
