#!/usr/bin/env bash
set -euo pipefail

TSV="manifests/procap_data.tsv"
OUTPUT_DIR="/fs/cbsuhy01/storage/yz2676/projects/1-ASD/NetFlow_ASD/Yutong_scripts/data/enhancer"
DIV_EXT_LEFT=200
DIV_EXT_RIGHT=200
SINGLE_EXT_LEFT=200
SINGLE_EXT_RIGHT=200
PROMOTER_FILE="/fs/cbsuhy01/storage/NetFlow_ASD/Week8_429/Promoters/promoters_1kb_tss_centered.coding.v43.bed"

# Count lines; assume header on line 1
LINES=$(wc -l < "$TSV")

# Array indexes: start at 2 to skip header; end at LINES
ARRAY_SPEC="2-${LINES}"

echo "Submitting array ${ARRAY_SPEC} for TSV: $TSV"
sbatch --array="${ARRAY_SPEC}" \
       src/1_enrichment/1.0.0_run_call_enhancers.sh \
       "$TSV" "$OUTPUT_DIR" \
       "$DIV_EXT_LEFT" "$DIV_EXT_RIGHT" "$SINGLE_EXT_LEFT" "$SINGLE_EXT_RIGHT" "$PROMOTER_FILE"