#!/usr/bin/env bash
set -euo pipefail

if [[ $# -lt 6 ]]; then
  echo "Usage: $0 CODING NONCODING_OR_NA INFO PRIORS SAMPLE_SIZES OUTPUT [extra TADA flags]" >&2
  exit 2
fi
repo_dir="$(cd "$(dirname "${BASH_SOURCE[0]}")/.." && pwd)"
coding="$1"; noncoding="$2"; info="$3"; priors="$4"; samples="$5"; output="$6"
shift 6
args=(--coding-file "$coding" --info-file "$info" --prior-file "$priors" --sample-file "$samples" --output-file "$output")
if [[ "$noncoding" != "NA" ]]; then args+=(--noncoding-file "$noncoding"); fi
Rscript "${repo_dir}/src/tada/run_tada.R" "${args[@]}" "$@"
