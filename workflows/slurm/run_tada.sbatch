#!/usr/bin/env bash
#SBATCH --job-name=tada
#SBATCH --cpus-per-task=1
#SBATCH --mem=8G
#SBATCH --time=04:00:00
#SBATCH --output=logs/tada_%j.out
set -euo pipefail

# Submit from the repository root after activating the documented environment:
# sbatch workflows/slurm/run_tada.sbatch CODING NONCODING_OR_NA INFO PRIORS SAMPLE OUTPUT [flags]
mkdir -p logs
bash workflows/run_tada.sh "$@"
