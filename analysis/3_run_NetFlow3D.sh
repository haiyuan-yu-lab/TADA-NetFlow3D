#!/bin/bash
#SBATCH --time=1-00:00:00
#SBATCH --mem=835G
#SBATCH --cpus-per-task=5
#SBATCH --job-name=netflow3d_%j
#SBATCH --mail-user=yz8667@princeton.edu
#SBATCH --mail-type=ALL
#SBATCH --output=logs/netflow3d_%j.out
#SBATCH --array=0-4

### Example:
# sbatchlog 3_run_NetFlow3D.sh /Genomics/argo/users/yz8667/Autism_network_integration/data/TADA/output/InhAll_DNVAll/ TADA_coding_UCSF_rare_Combined_Simple_Roulette_MR_ScaleByLOEUF.txt 
#  $1 = tada_output_path; $2 = tada_output_suffix; $3 = delta_permu [optional]; $4 = restart_prob [optional]; $5 = disrupt_edges [optional]; $6 = binary_interactome [optional]; $7 = expressed_genes [optional];

set -euo pipefail

eval "$(conda shell.bash hook)"
conda activate jupyter

# ---- Paths and parameters ----
netflow3d_input_path="/Genomics/argo/users/yz8667/Autism_network_integration/data/NetFlow3D/input/"
netflow3d_output_path="/Genomics/argo/users/yz8667/Autism_network_integration/data/NetFlow3D/output/"

# ---- Setup batch jobs ----
tada_output_path=$1
tada_output_suffix=$2
delta_permu=${3:-100}
restart_prob=${4:-0.5}
interface_dn=${5:-None}
binary_interactome=${6:-/Genomics/argo/users/yz8667/Autism_network_integration/scripts/NetFlow3D/metadata/HomoSapiens_binary_HINThq.txt}
expressed_genes=${7:-/Genomics/argo/users/yz8667/Autism_network_integration/data/ENCODE_GTEx/fetal_frontal_cortex/fetal_frontal_cortex_expressed_genes_TPM0.1.txt}

# ---- Prepare tags ----
interactome_tag="$(basename "$binary_interactome" .txt)"
if [[ -n "$interface_dn" && "$interface_dn" != "None" ]]; then
  disrupt_edges_tag="$(basename "$interface_dn" .txt)"
  disrupt_command="-e $interface_dn"
else
  disrupt_edges_tag="NoDisruptEdges"
  disrupt_command=""
fi
if [[ -n "$expressed_genes" ]]; then
  expressed_genes_tag="$(basename "$expressed_genes" .txt)"
  express_command="-X $expressed_genes" 
else
  expressed_genes_tag="NoExprFilter"
  express_command=""
fi

mkdir -p "$netflow3d_input_path" "$netflow3d_output_path" logs

# ---- Collect TADA files (sorted, deterministic) ----
# Note: uses grep -E for a clean filter.
echo "Collecting TADA files from: $tada_output_path with suffix: $tada_output_suffix"
mapfile -t TADA_FILES < <(
  find "$tada_output_path" -maxdepth 1 -type f -name "$tada_output_suffix" \
)

# ---- Select TADA file for this array task ----
n="${#TADA_FILES[@]}"
if [[ "$n" -eq 0 ]]; then
  echo "No selected TADA files found in: $tada_output_path" >&2
  exit 1
fi

if [[ "${SLURM_ARRAY_TASK_COUNT}" -lt "$n" ]]; then
  echo "ERROR: Found $n TADA files, but SLURM array has ${SLURM_ARRAY_TASK_COUNT} tasks." >&2
  echo "       Please adjust --array to 0-$((n-1))." >&2
  exit 1
fi

task_id="${SLURM_ARRAY_TASK_ID}"
tada_output="${TADA_FILES[$task_id]}"
tada_file_tag="$(basename "$tada_output_path")"_"$(basename "$tada_output" .txt)"
infile="${netflow3d_input_path}/${tada_file_tag}_NetFlow3D_Input.txt"

echo "Task ${SLURM_ARRAY_JOB_ID}_${task_id} running: $tada_output"
echo "Using disrupt edges file: $interface_dn"
echo "Using binary interactome file: $binary_interactome"
echo "Using expressed flag: ${express_command}"

# ---- Build NetFlow3D input (Gene, Score) ----
python - <<PY
import pandas as pd
import numpy as np

tada_output = "${tada_output}"
infile = "${infile}"

df = pd.read_csv(tada_output, sep="\\t").rename(columns={"gene_id": "Gene"})
df["Score"] = df["bf"].apply(lambda x: np.log(x))
df[["Gene", "Score"]].to_csv(infile, sep="\\t", header=True, index=False)
PY

# ---- Run NetFlow3D ----
python /Genomics/argo/users/yz8667/Autism_network_integration/scripts/NetFlow3D/NetFlow3D_ASD_version.py \
  -I "${tada_file_tag}_${disrupt_edges_tag}_${interactome_tag}_${expressed_genes_tag}_rp${restart_prob}_dp${delta_permu}" \
  -m "$infile" \
  -o "$netflow3d_output_path" \
  -t "${SLURM_CPUS_PER_TASK}" \
  -r "$restart_prob" \
  -D "$delta_permu" \
  -n "${binary_interactome}" \
  ${disrupt_command} \
  ${express_command}

echo "Done."