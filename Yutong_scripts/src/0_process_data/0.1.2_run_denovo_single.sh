#!/bin/bash

# For a small de novo file, a single pass is simpler and avoids array overhead.
# We still dedup → unique sites to drastically shrink the VEP work.

source src/0_process_data/vep.env
source src/0_process_data/vep_common.sh

echo "== De novo full-file annotation =="
out_vcf="$DNM_OUT_DIR/UCSF_20211001_dnms_WGS_hg38.no_match.vep.vcf.gz"

# Smaller fork/buffer keeps memory/IO tame for tiny inputs
annotate_full_vcf_as_sites "$DNM_VCF" "$out_vcf" "$TMP_DIR" 4 50000

echo "Done: $out_vcf"
