#!/bin/bash
#SBATCH --job-name=downsampled_TADA
#SBATCH -o logs/downsampled_TADA.%A_%a.out
#SBATCH -e logs/downsampled_TADA.%A_%a.err
#SBATCH --mail-type=ALL
#SBATCH --mail-user=yz2676@cornell.edu
#SBATCH -N 1
#SBATCH --ntasks-per-node=1
#SBATCH --mem=8G
#SBATCH --time=48:00:00
#SBATCH --nodelist=cbsuhy02

set -euo pipefail

# Examples:
#   sbatch --array=0-19 src/2_TADA/2.3_run_downsampled_TADA.sh \
#     --cohort SSC_KOR_STAR --n-rep 1000 --reps-per-task 50
#
#   sbatch --array=0-19 src/2_TADA/2.3_run_downsampled_TADA.sh \
#     --cohort SPARK --n-rep 1000 --reps-per-task 50

analysis_tag="InhAll_DNVAll"
cohort="SSC_KOR_STAR"
coding="FuRaw"
mut_col="Roulette_MR"
n_rep="1000"
reps_per_task="50"
rep_start=""
rep_end=""
seed="2026"
output_dir=""
coding_file=""
sample_file=""
prior_file="data/TADA/prior.txt"
info_file="data/TADA/info.txt"

while [[ $# -gt 0 ]]; do
  case "$1" in
    --analysis-tag)     analysis_tag="$2"; shift 2;;
    --cohort)           cohort="$2"; shift 2;;
    --coding)           coding="$2"; shift 2;;
    --mut-col)          mut_col="$2"; shift 2;;
    --n-rep)            n_rep="$2"; shift 2;;
    --reps-per-task)    reps_per_task="$2"; shift 2;;
    --rep-start)        rep_start="$2"; shift 2;;
    --rep-end)          rep_end="$2"; shift 2;;
    --seed)             seed="$2"; shift 2;;
    --output-dir)       output_dir="$2"; shift 2;;
    --coding-file)      coding_file="$2"; shift 2;;
    --sample-file)      sample_file="$2"; shift 2;;
    --prior-file)       prior_file="$2"; shift 2;;
    --info-file)        info_file="$2"; shift 2;;
    *) echo "Unknown arg: $1" >&2; exit 1;;
  esac
done

if [[ -z "${coding_file}" ]]; then
  coding_file="data/TADA/TADA_input_coding_${coding}_${cohort}.txt"
fi
if [[ -z "${sample_file}" ]]; then
  sample_file="data/TADA/sample_${coding}_${cohort}_${analysis_tag}.txt"
fi
if [[ -z "${output_dir}" ]]; then
  output_dir="data/TADA/downsampled/${analysis_tag}/${cohort}/coding_only"
fi

if [[ -n "${SLURM_ARRAY_TASK_ID:-}" ]]; then
  if (( reps_per_task < 1 )); then
    echo "--reps-per-task must be at least 1" >&2
    exit 1
  fi
  rep_start=$((SLURM_ARRAY_TASK_ID * reps_per_task))
  rep_end=$((rep_start + reps_per_task))
  if (( rep_end > n_rep )); then
    rep_end="${n_rep}"
  fi
fi

if [[ -z "${rep_start}" ]]; then
  rep_start="0"
fi
if [[ -z "${rep_end}" ]]; then
  rep_end="${n_rep}"
fi

mkdir -p "${output_dir}"
input_dir="${output_dir}/inputs"
result_dir="${output_dir}/TADA_results"
run_sample_file="${output_dir}/downsampled_sample.txt"
mkdir -p "${input_dir}" "${result_dir}"

echo "Generating downsampled TADA inputs:"
echo "  cohort=${cohort}"
echo "  replicate range=${rep_start}-${rep_end} of ${n_rep}"
echo "  coding_file=${coding_file}"
echo "  sample_file=${sample_file}"
echo "  output_dir=${output_dir}"

python "src/2_TADA/downsample_tada.py" \
  --coding-file "${coding_file}" \
  --sample-file "${sample_file}" \
  --input-dir "${input_dir}" \
  --run-sample-file "${run_sample_file}" \
  --n-rep "${n_rep}" \
  --rep-start "${rep_start}" \
  --rep-end "${rep_end}" \
  --seed "${seed}"

for (( replicate = rep_start; replicate < rep_end; replicate++ )); do
    printf -v input_file "%s/Downsampled_rep%04d_TADA_input.txt" "${input_dir}" "${replicate}"
    printf -v output_file "%s/Downsampled_rep%04d_TADA_output.txt" "${result_dir}" "${replicate}"

    echo "Running coding-only TADA replicate ${replicate}: ${output_file}"

    Rscript src/2_TADA/run_TADA.r \
        --coding_file "${input_file}" \
        --prior_file "${prior_file}" \
        --info_file "${info_file}" \
        --sample_file "${run_sample_file}" \
        --output_file "${output_file}" \
        --mut_col "${mut_col}" \
        --dn_coding \
        --in_coding

done
