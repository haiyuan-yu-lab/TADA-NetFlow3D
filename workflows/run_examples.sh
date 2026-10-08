#!/usr/bin/env bash
set -euo pipefail
repo_dir="$(cd "$(dirname "${BASH_SOURCE[0]}")/.." && pwd)"
# Run the two TADA examples and the prepared NetFlow3D subset independently.
bash "${repo_dir}/workflows/run_tada.sh" "${repo_dir}/config/example_config.yaml"
bash "${repo_dir}/workflows/run_netflow3d.sh" "${repo_dir}/config/example_config.yaml"
