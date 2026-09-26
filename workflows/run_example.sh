#!/usr/bin/env bash
set -euo pipefail

repo_dir="$(cd "$(dirname "${BASH_SOURCE[0]}")/.." && pwd)"
result_dir="${repo_dir}/results/example"
mkdir -p "${result_dir}/netflow3d"

python3 "${repo_dir}/scripts/validate_inputs.py" \
  --info "${repo_dir}/examples/input/gene_information.tsv" \
  --coding "${repo_dir}/examples/input/coding_counts.tsv" \
  --noncoding "${repo_dir}/examples/input/noncoding_counts.tsv" \
  --priors "${repo_dir}/examples/input/priors.tsv" \
  --sample-sizes "${repo_dir}/examples/input/sample_sizes.tsv"

Rscript "${repo_dir}/src/tada/run_tada.R" \
  --coding-file "${repo_dir}/examples/input/coding_counts.tsv" \
  --noncoding-file "${repo_dir}/examples/input/noncoding_counts.tsv" \
  --info-file "${repo_dir}/examples/input/gene_information.tsv" \
  --prior-file "${repo_dir}/examples/input/priors.tsv" \
  --sample-file "${repo_dir}/examples/input/sample_sizes.tsv" \
  --dn-em --in-em \
  --output-file "${result_dir}/tada_results.tsv"

python3 "${repo_dir}/scripts/tada_to_netflow3d.py" \
  --tada-results "${result_dir}/tada_results.tsv" \
  --output "${result_dir}/netflow3d_gene_scores.tsv"

python3 "${repo_dir}/src/netflow3d/run_netflow3d.py" \
  --job_name example \
  --gene_score "${result_dir}/netflow3d_gene_scores.tsv" \
  --binary_interactome "${repo_dir}/examples/input/interactome.tsv" \
  --id_mapping "${repo_dir}/examples/input/id_mapping.tsv" \
  --disrupted_interaction "${repo_dir}/examples/input/disrupted_interactions.tsv" \
  --threads 1 --delta_trial 2 --max_subnetwork_size 10 \
  --output_path "${result_dir}/netflow3d"

echo "Example completed: ${result_dir}"
