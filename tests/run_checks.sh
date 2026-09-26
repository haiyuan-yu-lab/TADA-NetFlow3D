#!/usr/bin/env bash
set -euo pipefail
repo_dir="$(cd "$(dirname "${BASH_SOURCE[0]}")/.." && pwd)"
export PYTHONPYCACHEPREFIX="${repo_dir}/.pycache"

python3 -m py_compile \
  "${repo_dir}/scripts/tada_to_netflow3d.py" \
  "${repo_dir}/scripts/validate_inputs.py" \
  "${repo_dir}/scripts/download_reference_data.py" \
  "${repo_dir}/src/netflow3d/run_netflow3d.py" \
  "${repo_dir}/src/netflow3d/network_functions.py"
python3 "${repo_dir}/tests/test_converter.py"
python3 "${repo_dir}/tests/test_references.py"
python3 "${repo_dir}/scripts/validate_inputs.py" \
  --info "${repo_dir}/examples/input/gene_information.tsv" \
  --coding "${repo_dir}/examples/input/coding_counts.tsv" \
  --noncoding "${repo_dir}/examples/input/noncoding_counts.tsv" \
  --priors "${repo_dir}/examples/input/priors.tsv" \
  --sample-sizes "${repo_dir}/examples/input/sample_sizes.tsv"

for script in "${repo_dir}"/workflows/*.sh "${repo_dir}"/workflows/slurm/*.sbatch; do
  bash -n "$script"
done

echo "Structural checks passed."
