#!/usr/bin/env bash
#SBATCH --job-name=netflow3d
#SBATCH --cpus-per-task=5
#SBATCH --mem=64G
#SBATCH --time=1-00:00:00
#SBATCH --output=logs/netflow3d_%j.out
set -euo pipefail

# Submit from the repository root and pass normal NetFlow3D CLI arguments.
mkdir -p logs
bash workflows/run_netflow3d.sh "$@"
