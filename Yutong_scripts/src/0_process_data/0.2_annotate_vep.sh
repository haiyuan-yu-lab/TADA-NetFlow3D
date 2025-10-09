#!/bin/bash
#SBATCH --job-name=annotate_vcf
#SBATCH -o logs/annotate_vcf.%j.out
#SBATCH -e logs/annotate_vcf.%j.err

source src/0_process_data/vep.env
source src/0_process_data/vep_common.sh
set -euo pipefail

FIELDS="CHROM,POS,REF,ALT,INFO/CSQ"

echo "== Annotate de novo vcf with VEP =="
dnm_vep="$DNM_OUT_DIR/UCSF_20211001_dnms_WGS_hg38.vep.vcf.gz"
tabix -f -p vcf "$dnm_vep"

bcftools view -h "$dnm_vep" | grep '^##INFO=<ID=CSQ' > $DNM_OUT_DIR/csq.hdr
bcftools query -f '%CHROM\t%POS\t%REF\t%ALT\t%INFO/CSQ\n' "$dnm_vep" \
| bgzip > $DNM_OUT_DIR/csq.4key.tsv.gz
tabix -s1 -b2 -e2 $DNM_OUT_DIR/csq.4key.tsv.gz

bcftools annotate \
  -a $DNM_OUT_DIR/csq.4key.tsv.gz \
  -c $FIELDS \
  -h $DNM_OUT_DIR/csq.hdr \
  -Oz -o $DNM_OUT_DIR/UCSF_20211001_dnms_WGS_hg38.with_vep.vcf.gz \
  "$DNM_VCF"

tabix -f -p vcf "$DNM_OUT_DIR/UCSF_20211001_dnms_WGS_hg38.with_vep.vcf.gz"
echo "Done: $DNM_OUT_DIR/UCSF_20211001_dnms_WGS_hg38.with_vep.vcf.gz"

rm $DNM_OUT_DIR/csq.hdr $DNM_OUT_DIR/csq.4key.tsv.gz $DNM_OUT_DIR/csq.4key.tsv.gz.tbi


echo "== Annotate inherited per chromosome =="
outs=()
tmpfiles=()
csq_hdr="$INH_OUT_DIR/csq.hdr"
bcftools view -h "$INH_OUT_DIR/All_inherited_WGS_hg38.${CHROMS[0]}.vep.vcf.gz" \
| grep '^##INFO=<ID=CSQ' > "$csq_hdr"
for CHR in "${CHROMS[@]}"; do
    inh_vep="$INH_OUT_DIR/All_inherited_WGS_hg38.${CHR}.vep.vcf.gz"
    out_chr="$INH_OUT_DIR/All_inherited_WGS_hg38.${CHR}.with_vep.vcf.gz"
    CHR_NC="${CHR//chr/}"

    echo "  * $CHR"
    tabix -f -p vcf "$inh_vep"

    csq_tsv="$INH_OUT_DIR/${CHR}.csq.4key.tsv.gz"
    bcftools query -f '%CHROM\t%POS\t%REF\t%ALT\t%INFO/CSQ\n' "$inh_vep" \
        | bgzip > "$csq_tsv"
    tabix -s1 -b2 -e2 "$csq_tsv"

    bcftools annotate \
        -a "$csq_tsv" \
        -c "$FIELDS" \
        -h "$csq_hdr" \
        -r "$CHR_NC" \
        -Oz -o "$out_chr" \
        "$INH_VCF"

    tabix -f -p vcf "$out_chr"
    outs+=("$out_chr")
    tmpfiles+=("$csq_hdr" "$csq_tsv" "$csq_tsv.tbi")
done

echo "== Concatenate inherited =="
bcftools concat -Oz -o "$INH_OUT_DIR/All_inherited_WGS_hg38.with_vep.vcf.gz" "${outs[@]}"
tabix -f -p vcf "$INH_OUT_DIR/All_inherited_WGS_hg38.with_vep.vcf.gz"
echo "Done: $INH_OUT_DIR/All_inherited_WGS_hg38.with_vep.vcf.gz"

# Clean up temps
rm -f "${tmpfiles[@]}"