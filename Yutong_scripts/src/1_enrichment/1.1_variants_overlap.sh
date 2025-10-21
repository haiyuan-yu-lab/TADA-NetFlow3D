#!/bin/bash
#SBATCH --job-name=variants_enrichment
#SBATCH -o logs/variants_enrichment.%A.out
#SBATCH -e logs/variants_enrichment.%A.err


python src/1_enrichment/enrichment.py \
    --vcf /fs/cbsuhy01/storage/yz2676/projects/1-ASD/NetFlow_ASD/Yutong_scripts/data/asdinfo/UCSF_20211001_dnms_WGS_hg38.filtered.vep.with_asdinfo.bcf \
    --output data/enrichment/UCSF_20211001_dnms_WGS_hg38.filtered.vep.with_asdinfo.vcf.gz
