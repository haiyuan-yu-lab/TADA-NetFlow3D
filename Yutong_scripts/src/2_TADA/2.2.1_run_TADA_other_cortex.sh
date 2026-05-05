#!/bin/bash
# prepare the file for heatmap

BASE="/fs/cbsuhy02/storage/yz2676/ASD/NetFlow_ASD"
TADA_DIR="data/TADA/InhAll_DNVAll/ALL"
TMP="/tmp/TADA_Other_Cerebral_Cortex_Elements.no_overlap.bed"

bedtools intersect -v \
  -a "$TADA_DIR/dnv/Cerebral Cortex/enh.merged.bed" \
  -b "$TADA_DIR/dnv/Deconv Microglia/enh.merged.bed" \
     "$TADA_DIR/dnv/Deconv Astrocytes/enh.merged.bed" \
     "$TADA_DIR/dnv/Deconv Excitatory Neurons/enh.merged.bed" \
  > "$TMP"

awk 'NR==FNR {ids[$4]=1; next} FNR==1 || ($1 in ids)' \
  "$TMP" \
  "$TADA_DIR/input/TADA_Cerebral_Cortex_roulette_rd_all.txt" \
  > "$TADA_DIR/input/TADA_Other_Cerebral_Cortex_roulette_rd_all.txt"

awk 'NR==FNR {ids[$4]=1; next} FNR==1 || ($1 in ids)' \
  "$TMP" \
  "$TADA_DIR/input/TADA_Cerebral_Cortex_roulette_rd_control.txt" \
  > "$TADA_DIR/input/TADA_Other_Cerebral_Cortex_roulette_rd_control.txt"

Rscript src/2_TADA/run_TADA.r \
  --coding_file $BASE/Yutong_scripts/data/TADA/TADA_input_coding_FuRaw_ALL.txt \
  --prior_file $BASE/prior.txt \
  --info_file $BASE/info.txt \
  --sample_file $BASE/Yutong_scripts/data/TADA/sample_ALL_InhAll_DNVAll.txt \
  --output_file $TADA_DIR/TADA_results/TADA_FuRaw_Other_Cerebral_Cortex_Roulette_MR_rd_all.txt \
  --noncoding_file $TADA_DIR/input/TADA_Other_Cerebral_Cortex_roulette_rd_all.txt \
  --mut_col Roulette_MR \
  --dn_coding --in_coding --dn_noncoding --in_noncoding

Rscript src/2_TADA/run_TADA.r \
  --coding_file $BASE/Yutong_scripts/data/TADA/TADA_input_coding_FuRaw_ALL.txt \
  --prior_file $BASE/prior.txt \
  --info_file $BASE/info.txt \
  --sample_file $BASE/Yutong_scripts/data/TADA/sample_ALL_InhAll_DNVAll.txt \
  --output_file $TADA_DIR/TADA_results/TADA_FuRaw_Other_Cerebral_Cortex_Roulette_MR_rd_control.txt \
  --noncoding_file $TADA_DIR/input/TADA_Other_Cerebral_Cortex_roulette_rd_control.txt \
  --mut_col Roulette_MR \
  --dn_coding --in_coding --dn_noncoding --in_noncoding

rm "$TMP"