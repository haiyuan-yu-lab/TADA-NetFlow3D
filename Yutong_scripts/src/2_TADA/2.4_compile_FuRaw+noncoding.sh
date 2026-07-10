#!/bin/bash
#SBATCH --job-name=compile_fig3_TADA
#SBATCH -o logs/compile_fig3_TADA.%A.out
#SBATCH -e logs/compile_fig3_TADA.%A.err
#SBATCH --mail-type=ALL
#SBATCH --mail-user=yz2676@cornell.edu
#SBATCH -N 1
#SBATCH --ntasks-per-node=1
#SBATCH --mem=8G
#SBATCH --time=04:00:00

set -euo pipefail

analysis_tag="InhAll_DNVAll"
n_rep="1000"
output_prefix="Fig3_downsampled"
downsample_dir=""
noncoding_file=""
netflow3d_input_dir="data/netflow3d/input"

while [[ $# -gt 0 ]]; do
  case "$1" in
    --analysis-tag)           analysis_tag="$2"; shift 2;;
    --n-rep)                  n_rep="$2"; shift 2;;
    --output-prefix)          output_prefix="$2"; shift 2;;
    --downsample-dir)         downsample_dir="$2"; shift 2;;
    --noncoding-file)         noncoding_file="$2"; shift 2;;
    --netflow3d-input-dir)    netflow3d_input_dir="$2"; shift 2;;
    *) echo "Unknown arg: $1" >&2; exit 1;;
  esac
done

if [[ -z "${downsample_dir}" ]]; then
  downsample_dir="data/TADA/downsampled/${analysis_tag}"
fi
if [[ -z "${noncoding_file}" ]]; then
  noncoding_file="data/TADA/${analysis_tag}/ALL/TADA_results/TADA_NoncodingOnly_Combined_Simple_Roulette_MR_rd_control.txt"
fi

mkdir -p "${downsample_dir}" "${netflow3d_input_dir}"

echo "Compiling Fig3 downsampled TADA:"
echo "  downsample_dir=${downsample_dir}"
echo "  noncoding_file=${noncoding_file}"
echo "  netflow3d_input_dir=${netflow3d_input_dir}"
echo "  n_rep=${n_rep}"
echo "  output_prefix=${output_prefix}"

python "src/2_TADA/compile_fig3_tada.py" \
  --downsample-dir "${downsample_dir}" \
  --noncoding-file "${noncoding_file}" \
  --netflow3d-input-dir "${netflow3d_input_dir}" \
  --n-rep "${n_rep}" \
  --output-prefix "${output_prefix}"
