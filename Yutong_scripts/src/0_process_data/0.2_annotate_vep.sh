#!/bin/bash
#SBATCH --job-name=annotate_vcf
#SBATCH -o logs/annotate_vcf.%j.out
#SBATCH -e logs/annotate_vcf.%j.err

source src/0_process_data/vep.env
source src/0_process_data/vep_common.sh
set -euo pipefail

FIELDS="INFO/CSQ"

echo "== Annotate de novo vcf with VEP =="
dnm_vep="$DNM_OUT_DIR/UCSF_20211001_dnms_WGS_hg38.vep.vcf.gz"
tabix -f -p vcf "$dnm_vep"

bcftools annotate \
  -a "$dnm_vep" \
  -c "$FIELDS" \
  -Oz -o "$DNM_OUT_DIR/UCSF_20211001_dnms_WGS_hg38.with_vep.vcf.gz" \
  "$DNM_VCF"

tabix -f -p vcf "$DNM_OUT_DIR/UCSF_20211001_dnms_WGS_hg38.with_vep.vcf.gz"
echo "Done: $DNM_OUT_DIR/UCSF_20211001_dnms_WGS_hg38.with_vep.vcf.gz"


echo "== Annotate inherited per chromosome =="
outs=()
for CHR in "${CHROMS[@]}"; do
    inh_vep="$INH_OUT_DIR/All_inherited_WGS_hg38.${CHR}.vep.vcf.gz"
    out_chr="$INH_OUT_DIR/All_inherited_WGS_hg38.${CHR}.with_vep.vcf.gz"
    CHR_NC="${CHR//chr/}"

    echo "  * $CHR"
    tabix -f -p vcf "$inh_vep"

    bcftools annotate \
        -a "$inh_vep" \
        -c "$FIELDS" \
        -r "$CHR_NC" \
        -Oz -o "$out_chr" \
        "$INH_VCF"

    tabix -f -p vcf "$out_chr"
    outs+=("$out_chr")
done

echo "== Concatenate inherited =="
bcftools concat -Oz -o "$INH_OUT_DIR/All_inherited_WGS_hg38.with_vep.vcf.gz" "${outs[@]}"
tabix -f -p vcf "$INH_OUT_DIR/All_inherited_WGS_hg38.with_vep.vcf.gz"
echo "Done: $INH_OUT_DIR/All_inherited_WGS_hg38.with_vep.vcf.gz"
