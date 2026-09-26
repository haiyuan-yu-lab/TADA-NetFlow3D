#!/bin/bash
#SBATCH --job-name=vep_anno
#SBATCH --partition=pritykinlab,main
#SBATCH --cpus-per-task=12
#SBATCH --mem=64G
#SBATCH --time=3-00:00:00
#SBATCH --output=logs/vep_anno.%j.out
#SBATCH --error=logs/vep_anno.%j.err
#SBATCH --mail-user=yz8667@princeton.edu
#SBATCH --mail-type=ALL

### Example:
# sbatchlog 0_run_VEP.sh /Genomics/argo/users/yz8667/Autism_network_integration/data/PPI_interface/UCSF_Fu_de_novo_variants.vcf /Genomics/argo/users/yz8667/Autism_network_integration/data/PPI_interface/
####
#  $1 = input_vcf; $2 = output_dir;


set -euo pipefail

module load vep/110

mkdir -p logs

ALPHAMISSENSE=/Genomics/argo/users/yz8667/.vep/Plugins/data/AlphaMissense_hg38.tsv.gz
REVEL=/Genomics/argo/users/yz8667/.vep/Plugins/data/new_tabbed_revel_grch38.tsv.gz
GNOMAD=/Genomics/argo/users/yz8667/.vep/Plugins/data/gnomad.ch.genomesv3.tabbed.tsv.gz

input_vcf=$1
outdir=$2

vep \
  -i ${input_vcf} \
  -o ${outdir}/$(basename ${input_vcf} .vcf).vep.txt \
  --cache \
  --offline \
  --assembly GRCh38 \
  --fork ${SLURM_CPUS_PER_TASK} \
  --buffer_size 10000 \
  --protein \
  --symbol \
  --canonical \
  --variant_class \
  --plugin AlphaMissense,file=${ALPHAMISSENSE} \
  --plugin REVEL,file=${REVEL} \
  --af_gnomadg \
  --af_gnomade \
  --force_overwrite \
  --no_stats

  # --plugin gnomADc,${GNOMAD} \
  # --hgvs \