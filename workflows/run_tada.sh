#!/usr/bin/env bash
set -euo pipefail
# Usage: bash workflows/run_tada.sh config/example_config.yaml [--cohort SPARK]

repo_dir="$(cd "$(dirname "${BASH_SOURCE[0]}")/.." && pwd)"
python3 "${repo_dir}/scripts/run_tada_from_config.py" "$@"
