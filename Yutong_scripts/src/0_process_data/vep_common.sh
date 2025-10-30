#!/usr/bin/env bash

# Annotate inherited VCF files with VEP
# to get the functional annotations:
# Consequence, CADD, phyloP, phastCons, REVEL, gnomADg_AF, etc.

# Note: first download VEP cache for human GRCh38: 
# vep_install -a cf -s homo_sapiens -y GRCh38 -c data/vep_cache --CONVERT --NO_UPDATE --USE_HTTPS_PROTO

# Plugins: 
# Please refer to https://useast.ensembl.org/info/docs/tools/vep/script/vep_plugins.html
# gnomAD v3.1 genome at /fs/cbsuhyfs1/storage/resources/gnomad/r3.1/hg38/genome
# CADD v1.7 SNV+InDels at data/supplemental/cadd
# REVEL v1.3 at data/supplemental/revel

set -euo pipefail

require_file() {
  local f="$1"
  [[ -f "$f" ]] || { echo "Missing file: $f" >&2; exit 1; }
}

# Create a unique-sites VCF (sorted, bgzipped, indexed) for a given input (optionally a single chrom)
#   $1 = input VCF[.gz]
#   $2 = output sites VCF.gz
#   $3 = (optional) region spec like "chr3" or "3"
make_unique_sites_vcf() {
  local in_vcf="$1" out_sites="$2" region="${3:-}"

  require_file "$in_vcf"
  local region_arg=()
  [[ -n "$region" ]] && region_arg=(-r "$region")

  # Extract CHROM,POS,REF,ALT, deduplicate, rebuild minimal VCF header, bgzip + index
  bcftools view "${region_arg[@]}" "$in_vcf" \
  | bcftools query -f '%CHROM\t%POS\t%REF\t%ALT\n' \
  | LC_ALL=C sort -u -k1,1n -k2,2n -k3,3 -k4,4 \
  | awk 'BEGIN{OFS="\t"}{print $1,$2,".",$3,$4,".",".","."}' \
  | (echo '##fileformat=VCFv4.2'; echo -e '#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO'; cat) \
  | bgzip > "$out_sites"

  tabix -f -p vcf "$out_sites"
}

# Build gnomAD per-chrom path (gnomAD files are per-chr and have chr prefix)
#   $1 = chrom label with "chr" (e.g., chr1)
gnomad_sites_vcf_for_chr() {
  local chr="$1"
  echo "$GNOMAD_DIR/gnomad.genomes.v3.1.2.sites.${chr}.vcf.bgz"
}

# Run VEP with consistent plugins and customs
#   $1 = input VCF.gz (sites)
#   $2 = output VCF.gz
#   $3 = fork (threads)
#   $4 = buffer_size
run_vep_sites() {
  local in_vcf="$1" out_vcf="$2" fork="$3" buf="$4"
  require_file "$in_vcf"
  require_file "$REVEL_TSV"
  require_file "$CADD_SNV"
  require_file "$CADD_INDEL"
  require_file "$PHYLOP_BW"
  require_file "$PHCONS_BW"

  # --- determine chromosomes present in input ---
  # normalize to {1..22,X,Y}, and build "chr*" names
  mapfile -t _chroms < <(
    zcat "$in_vcf" \
      | awk '$1!~/^#/{print $1}' \
      | awk '!seen[$0]++' \
      | awk '$0~/^([1-9]|1[0-9]|2[0-2]|X|Y)$/ {print "chr"$0}'
  )
  # --- build gnomAD --custom args ---
  local -a GNOMAD_CUSTOMS=()
  if ((${#_chroms[@]}==1)); then
    local chrom="${_chroms[0]}"
    local gnomad_vcf="$(gnomad_sites_vcf_for_chr "$chrom")"
    require_file "$gnomad_vcf"
    GNOMAD_CUSTOMS=( --custom "file=$gnomad_vcf,short_name=gnomADg_v3,format=vcf,type=exact,fields=AF" )
  elif ((${#_chroms[@]}>1)); then
    local c v
    for c in "${_chroms[@]}"; do
      v="$(gnomad_sites_vcf_for_chr "$c")"
      require_file "$v"
      GNOMAD_CUSTOMS+=( --custom "file=$v,short_name=gnomADg_v3,format=vcf,type=exact,fields=AF" )
    done
  fi

  vep -i "$in_vcf" \
      --format vcf \
      --offline --cache --dir "$CACHE_DIR" \
      --assembly GRCh38 \
      --vcf --force_overwrite --no_stats \
      --compress_output bgzip \
      --fork "$fork" --buffer_size "$buf" \
      --sift b --polyphen b \
      --plugin REVEL,file="$REVEL_TSV",no_match=1 \
      --plugin CADD,snv="$CADD_SNV",indels="$CADD_INDEL" \
      "${GNOMAD_CUSTOMS[@]}" \
      --custom file="$PHYLOP_BW",short_name=phyloP100way,format=bigwig,type=exact \
      --custom file="$PHCONS_BW",short_name=phastCons100way,format=bigwig,type=exact \
      --fields "$VEP_FIELDS" \
      -o "$out_vcf" \
      --verbose
}

# Convenience wrapper: annotate a full VCF (all variants) in one shot (for small files)
#   $1 = input VCF.gz (full variants, not only sites)
#   $2 = output VCF.gz
#   $3 = tmp dir
#   $4 = fork
#   $5 = buffer_size
annotate_full_vcf_as_sites() {
  local full_vcf="$1" out_vcf="$2" tmp="$3" fork="$4" buf="$5"
  local base="$(basename "$out_vcf" .vcf.gz)"
  local sites="$tmp/${base}.unique.sites.vcf.gz"
  make_unique_sites_vcf "$full_vcf" "$sites"
  run_vep_sites "$sites" "$out_vcf" "$fork" "$buf"
}
