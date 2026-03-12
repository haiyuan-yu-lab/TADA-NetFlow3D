#!/bin/bash
#SBATCH --job-name=vep_denovo
#SBATCH -o logs/vep_denovo.%A.out
#SBATCH -e logs/vep_denovo.%A.err
#SBATCH --mail-type=ALL
#SBATCH --mail-user=yz2676@cornell.edu

#SBATCH -N 1
#SBATCH --cpus-per-task=4
#SBATCH --mem=64G
#SBATCH --time=24:00:00
#SBATCH --nodelist=cbsuhy02

# For a small de novo file, a single pass is simpler and avoids array overhead.
# We still dedup → unique sites to drastically shrink the VEP work.

source src/0_process_data/vep.env
source src/0_process_data/vep_common.sh

echo "== De novo full-file annotation =="
out_vcf="$DNM_OUT_DIR/UCSF_20211001_dnms_WGS_hg38.vep.vcf.gz"

# Smaller fork/buffer keeps memory/IO tame for tiny inputs
annotate_full_vcf_as_sites "$DNM_VCF" "$out_vcf" "$TMP_DIR" 4 50000

echo "Done: $out_vcf"
