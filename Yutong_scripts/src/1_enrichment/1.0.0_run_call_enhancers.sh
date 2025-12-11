#!/bin/bash
#SBATCH --job-name=call_enhancers
#SBATCH -o logs/call_enhancers.%A_%a.out
#SBATCH -e logs/call_enhancers.%A_%a.err

set -euo pipefail

TSV="${1:?TSV path is required}"
OUTPUT_DIR="${2:?OUTPUT_DIR is required}"
DIV_EXT_LEFT="${3:-60}"
DIV_EXT_RIGHT="${4:-60}"
SINGLE_EXT_LEFT="${5:-60}"
SINGLE_EXT_RIGHT="${6:-60}"
PROMOTER_FILE="${7:-}"
ELEMENT_NAME="${8:-all_peaks_element}"
shift 8
EXTRA_ARGS=("$@") 

mkdir -p "$OUTPUT_DIR"

# Use SLURM_ARRAY_TASK_ID to pick the TSV line.
# We assume the TSV has a header on line 1, so the array will start at 2.
line_raw="$(sed -n "${SLURM_ARRAY_TASK_ID}p" "$TSV" )"

# Parse columns: cell_type  bw_pl  bw_mn
IFS=$'\t' read -r CELL_TYPE BW_PL BW_MN <<< "$line_raw"

# Basic checks
if [[ -z "$CELL_TYPE" || -z "$BW_PL" || -z "$BW_MN" ]]; then
    echo "ERROR: Missing fields in TSV line $SLURM_ARRAY_TASK_ID"
    echo "Line content: <$line_raw>"
    exit 1
fi
if [[ ! -f "$BW_PL" ]]; then echo "ERROR: bw_pl not found: $BW_PL"; exit 1; fi
if [[ ! -f "$BW_MN" ]]; then echo "ERROR: bw_mn not found: $BW_MN"; exit 1; fi

log(){ echo "[$(date '+%F %T')] $*" >&2; }

log "Running cell_type: $CELL_TYPE"
log "bw_pl: $BW_PL"
log "bw_mn: $BW_MN"

ENH_DIR="$OUTPUT_DIR/$CELL_TYPE"
mkdir -p "$ENH_DIR"
srun python src/1_enrichment/call_enhancers.py \
    --bw-pl "$BW_PL" \
    --bw-mn "$BW_MN" \
    --output-path "$ENH_DIR" \
    --job-name "$CELL_TYPE" \
    --element-name "$ELEMENT_NAME" \
    --div-ext-left "$DIV_EXT_LEFT" \
    --div-ext-right "$DIV_EXT_RIGHT" \
    --single-ext-left "$SINGLE_EXT_LEFT" \
    --single-ext-right "$SINGLE_EXT_RIGHT" \
    ${PROMOTER_FILE:+--promoter-file "$PROMOTER_FILE"} \
    "${EXTRA_ARGS[@]}"

log "Finished cell_type: $CELL_TYPE"