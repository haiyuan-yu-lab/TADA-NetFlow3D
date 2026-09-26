#!/usr/bin/env bash
set -euo pipefail
repo_dir="$(cd "$(dirname "${BASH_SOURCE[0]}")/.." && pwd)"
python3 "${repo_dir}/src/netflow3d/run_netflow3d.py" "$@"
