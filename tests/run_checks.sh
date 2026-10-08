#!/usr/bin/env bash
set -euo pipefail
repo_dir="$(cd "$(dirname "${BASH_SOURCE[0]}")/.." && pwd)"
export PYTHONPYCACHEPREFIX="${repo_dir}/.pycache"

python3 -m py_compile \
  "${repo_dir}/scripts/tada_to_netflow3d.py" \
  "${repo_dir}/scripts/validate_inputs.py" \
  "${repo_dir}/scripts/run_tada_from_config.py" \
  "${repo_dir}/scripts/run_netflow_from_config.py" \
  "${repo_dir}/src/netflow3d/run_netflow3d.py" \
  "${repo_dir}/src/netflow3d/diffusion.py"
python3 "${repo_dir}/tests/test_converter.py"
python3 "${repo_dir}/tests/test_references.py"
for cohort in SSC_KOR_STAR SPARK; do
  python3 "${repo_dir}/scripts/validate_inputs.py" \
    --info "${repo_dir}/examples/input/tada/shared/gene_information.tsv" \
    --coding "${repo_dir}/examples/input/tada/${cohort}/coding_counts.tsv" \
    --noncoding "${repo_dir}/examples/input/tada/${cohort}/noncoding_counts.tsv" \
    --priors "${repo_dir}/examples/input/tada/shared/priors.tsv" \
    --sample-sizes "${repo_dir}/examples/input/tada/${cohort}/sample_sizes.tsv" \
    --mut-col mut.em
done

python3 "${repo_dir}/scripts/validate_inputs.py" \
  --gene-score "${repo_dir}/examples/input/netflow3d/gene_scores.tsv" \
  --disrupted-interactions "${repo_dir}/examples/input/netflow3d/disrupted_interactions.tsv"

for script in "${repo_dir}"/workflows/*.sh; do
  bash -n "$script"
done

echo "Structural checks passed."
