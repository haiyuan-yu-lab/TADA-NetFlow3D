#!/bin/bash
#SBATCH --job-name=untrans_chunk
#SBATCH -o logs/untrans_chunk.%A_%a.out
#SBATCH -e logs/untrans_chunk.%A_%a.err
#SBATCH -N 1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=1
#SBATCH --mem=36G
#SBATCH --time=48:00:00
#SBATCH --nodelist=cbsuhy02
#SBATCH --array=0-24

set -euo pipefail

log(){ echo "[$(date '+%F %T')] $*" >&2; }

# ---- CONFIG ----
ASD_DIR="/fs/cbsuhy01/storage/yz2676/projects/1-ASD/NetFlow_ASD/Yutong_scripts/data/asdinfo"
CHUNK_DIR="$ASD_DIR/trio_chunks"

# These should match what your main pipeline produced
FINAL_PARENT_BCF="$ASD_DIR/UCSF_parents_WGS_hg38.includesimplex.codingonly.filtered.vep.with_asdinfo.bcf"
INH_FOR_UNTRANS_BCF="$ASD_DIR/All_inherited_WGS_hg38.includesimplex.codingonly.filtered.vep.with_asdinfo.bcf"

# TSV index paths created by build_sample_indexed_tsv()
TX_TSV="$ASD_DIR/All_inherited_WGS_hg38.includesimplex.codingonly.filtered.vep.with_asdinfo.by_child.tsv.gz"
PR_TSV="$ASD_DIR/UCSF_parents_WGS_hg38.includesimplex.codingonly.filtered.vep.with_asdinfo.by_parent.tsv.gz"

OUT_PREFIX="$ASD_DIR/UCSF_parents_WGS_hg38.includesimplex.codingonly"
OUT_BODY_DIR="$ASD_DIR/untrans_bodies"
mkdir -p "$OUT_BODY_DIR"
# -----------------------------

CHUNK=$(printf "%s/trios.%04d" "$CHUNK_DIR" "$SLURM_ARRAY_TASK_ID")
OUT_BODY=$(printf "%s/%s.untransmitted.body.%04d.vcf" \
  "$OUT_BODY_DIR" "$(basename "$OUT_PREFIX")" "$SLURM_ARRAY_TASK_ID")

[[ -s "$CHUNK" ]] || { log "Missing chunk: $CHUNK"; exit 1; }
[[ -s "$TX_TSV" ]] || { log "Missing transmitted TSV: $TX_TSV"; exit 1; }
[[ -s "$PR_TSV" ]] || { log "Missing parent TSV: $PR_TSV"; exit 1; }

log "Chunk: $CHUNK"
log "Writing: $OUT_BODY"

python3 src/0_process_data/make_untransmitted.py \
  --trios "$CHUNK" \
  --transmitted-tsv "$TX_TSV" \
  --parent-tsv "$PR_TSV" \
  --out-body "$OUT_BODY"

log "Done chunk $SLURM_ARRAY_TASK_ID"
