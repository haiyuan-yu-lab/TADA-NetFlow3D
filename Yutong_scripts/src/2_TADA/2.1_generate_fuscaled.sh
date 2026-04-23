#!/bin/bash

set -euo pipefail

base_dir=/fs/cbsuhy02/storage/yz2676/ASD/NetFlow_ASD
fu_raw="${base_dir}/TADA_input_coding_FuRaw.txt"
sample_tag="InhAll_DNVAll"
cohorts=("SSC_KOR_STAR" "SPARK")

python src/2_TADA/scale_fuscaled.py \
    --input "${fu_raw}" \
    --sample-tag "${sample_tag}" \
    --cohort "${cohorts[@]}" \
    --output-dir "${base_dir}"
