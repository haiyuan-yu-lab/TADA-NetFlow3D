#!/bin/bash
#SBATCH --job-name=vep_inh_chr
#SBATCH -o logs/vep_inh_chr.%A_%a.out
#SBATCH -e logs/vep_inh_chr.%A_%a.err
#SBATCH --array=1-24    # 1..22,X,Y

source src/0_process_data/vep.env
source src/0_process_data/vep_common.sh

idx=$((SLURM_ARRAY_TASK_ID-1))
CHR="${CHROMS[$idx]}"
CHR_NC="${CHROMS_NO_CHR[$idx]}"

echo "== Inherited per-chrom: $CHR =="
sites_vcf="$TMP_DIR/inherited.${CHR}.unique.sites.vcf.gz"
out_vcf="$INH_OUT_DIR/All_inherited_WGS_hg38.${CHR}.vep.vcf.gz"

# Build per-chrom unique sites from the full inherited VCF
make_unique_sites_vcf "$INH_VCF" "$sites_vcf" "$CHR_NC"

# VEP (bigger file → keep higher fork/buffer)
run_vep_sites "$sites_vcf" "$out_vcf" 16 80000

echo "Done: $out_vcf"
