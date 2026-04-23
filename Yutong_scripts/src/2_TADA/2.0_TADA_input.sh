#!/bin/bash
#SBATCH --job-name=generate_TADA_input
#SBATCH -o logs/generate_TADA_input.%A.out
#SBATCH -e logs/generate_TADA_input.%A.err
#SBATCH --mail-type=ALL
#SBATCH --mail-user=yz2676@cornell.edu

#SBATCH -N 1
#SBATCH --ntasks-per-node=1
#SBATCH --mem=32G
#SBATCH --time=48:00:00
#SBATCH --nodelist=cbsuhy02

COHORT="InhAll_DNVAll" # InhMult_DNVSimp / InhAll_DNVAll
samples="SSC_KOR_STAR" # SSC_KOR_STAR / SPARK / ALL
rare_damaging_mode="all" # all / control

while [[ $# -gt 0 ]]; do
  case "$1" in
    --cohort)               COHORT="$2"; shift 2;;
    --samples)              samples="$2"; shift 2;;
    --rare-damaging-mode)   rare_damaging_mode="$2"; shift 2;;
    *) echo "Unknown arg: $1" >&2; exit 1;;
  esac
done

var_prefilter_dir="/fs/cbsuhy01/storage/yz2676/projects/1-ASD/NetFlow_ASD/Yutong_scripts/data/asdinfo"
outbase=/fs/cbsuhy02/storage/yz2676/ASD/NetFlow_ASD/Yutong_scripts/data/TADA
outdir="${outbase}/${COHORT}"
mkdir -p "$outdir"

prefix_dnv="UCSF_dnms.midconseq_damagingmissense.filtered.vep.with_asdinfo"
prefix_dnv_all="UCSF_20211001_dnms_WGS_hg38.with_multiplex.filtered.vep.with_asdinfo"
prefix_inh_mult="All_inherited_WGS_hg38.filtered.vep.with_asdinfo"
# prefix_inh_simp="All_inherited_WGS_hg38.onlysimplex.filtered.vep.with_asdinfo"
prefix_inh_all="All_inherited_WGS_hg38.includesimplex.filtered.vep.with_asdinfo"
# prefix_untrans_mult="UCSF_parents_WGS_hg38.onlymultiplex.untransmitted.dedup"
prefix_untrans_all="UCSF_parents_WGS_hg38.includesimplex.untransmitted.dedup.with_batch"

# de novo
dnv_simp_vcf="${var_prefilter_dir}/${prefix_dnv}.bcf"
dnv_all_vcf="${var_prefilter_dir}/${prefix_dnv_all}.bcf"

# transmitted
inh_mult_vcf="${var_prefilter_dir}/${prefix_inh_mult}.bcf"
# inh_simp_vcf="${var_prefilter_dir}/${prefix_inh_simp}.bcf"
inh_all_vcf="${var_prefilter_dir}/${prefix_inh_all}.bcf"

# untransmitted
untrans_mult_vcf="${var_prefilter_dir}/${prefix_untrans_mult}.bcf"
untrans_all_vcf="${var_prefilter_dir}/${prefix_untrans_all}.bcf"

case "$COHORT" in
  InhMult_DNVSimp)
    dnv_vcf="$dnv_simp_vcf"
    inh_vcf="$inh_mult_vcf"
    untrans_vcf="$untrans_mult_vcf"
    ;;
  InhAll_DNVAll)
    dnv_vcf="$dnv_all_vcf"
    inh_vcf="$inh_all_vcf"
    untrans_vcf="$untrans_all_vcf"
    ;;
  *)
    echo "Unknown COHORT: $COHORT" >&2; exit 1;;
esac


config=/fs/cbsuhy02/storage/yz2676/ASD/NetFlow_ASD/Yutong_scripts/manifests/TADA_input.json

# bedtools slop -i /fs/cbsuhy01/storage/NetFlow_ASD/Week8_429/Promoters/promoters_1kb_tss_centered.coding.v43.bed \
#     -g /fs/cbsuhy01/storage/yz2676/data/DetailedProCapNet/wdir/reference/hg38.22XY.fasta.fai \
#     -l -500 -r -499 \
#     > /fs/cbsuhy02/storage/yz2676/ASD/NetFlow_ASD/Yutong_scripts/data/reference/promoters_tss_centered.coding.v43.bed

coding_tss=/fs/cbsuhy02/storage/yz2676/ASD/NetFlow_ASD/Yutong_scripts/data/reference/promoters_tss_centered.coding.v43.bed
enhancer_dir=/fs/cbsuhy01/storage/yz2676/projects/1-ASD/NetFlow_ASD/Yutong_scripts/data/enhancer

mutation_rate_file="/fs/cbsuhy01/storage/NetFlow_ASD/Mutability/20230919_samocha14_calcBaseLineRate_EG_mapping_allsample.bed"

echo "Generating TADA input files in $outdir for cohort: $COHORT, samples: $samples"
echo "DNV VCF: $dnv_vcf"
echo "Inherited VCF: $inh_vcf"
echo "Untransmitted VCF: $untrans_vcf"

python src/2_TADA/generate_TADA_input.py \
    --outdir "$outdir" \
    --config "$config" \
    build \
    --vcf "$dnv_vcf" \
    --type "dnv" \
    --cohort "$samples" \
    --tss-ref "$coding_tss" \
    --enhancer-dir "$enhancer_dir" \
    --enhancer-prefix "bidirectional_peaks_element_2kb_200bp_e" \
    --gnomad-null-only \
    --cadd-min 12 \
    --phylop-min 2 \
    --phastcons-min 0.5

python src/2_TADA/generate_TADA_input.py \
    --outdir "$outdir" \
    --config "$config" \
    build \
    --vcf "$inh_vcf" \
    --type "inh" \
    --cohort "$samples" \
    --tss-ref "$coding_tss" \
    --enhancer-dir "$enhancer_dir" \
    --enhancer-prefix "bidirectional_peaks_element_2kb_200bp_e" \
    --gnomad-null-only \
    --cadd-min 12 \
    --phylop-min 2 \
    --phastcons-min 0.5

python src/2_TADA/generate_TADA_input.py \
    --outdir "$outdir" \
    --config "$config" \
    build \
    --vcf "$untrans_vcf" \
    --type "untrans" \
    --cohort "$samples" \
    --tss-ref "$coding_tss" \
    --enhancer-dir "$enhancer_dir" \
    --enhancer-prefix "bidirectional_peaks_element_2kb_200bp_e" \
    --gnomad-null-only \
    --cadd-min 12 \
    --phylop-min 2 \
    --phastcons-min 0.5

echo "Combining dnv and inh TADA input files..."
rare_damaging_case=$(bcftools view -H "$outdir/$samples/dnv/case.vcf.gz" | wc -l)
rare_damaging_control=$(bcftools view -H "$outdir/$samples/dnv/ctrl.vcf.gz" | wc -l)

#  Match the denominator to the same cohort split used when generating
# the cohort-specific case/control VCFs above.
case "$samples" in
    SPARK)
        cohort_filter='INFO/BATCH="SPARK_WGS_1" || INFO/BATCH="SPARK_WGS_2" || INFO/BATCH="SPARK_WGS_3" || INFO/BATCH="SPARK"'
        ;;
    SSC_KOR_STAR)
        cohort_filter='INFO/BATCH="SSC_WGS_9209" || INFO/BATCH="KOR_STAR" || INFO/BATCH="KOR-STAR"'
        ;;
    ALL)
        cohort_filter=""
        ;;
    *)
        echo "Unknown samples cohort: $samples" >&2
        exit 1
        ;;
esac

if [[ "$rare_damaging_mode" == "all" ]]; then
    echo "Calculating rare damaging variant rate using all DNVs"
    if [[ -n "$cohort_filter" ]]; then
        dnv_var_count=$(bcftools view -i "$cohort_filter" -H "$dnv_vcf" | wc -l)
    else 
        dnv_var_count=$(bcftools view -H "$dnv_vcf" | wc -l)
    fi
    rare_damaging_rate=$(awk -v num=$((rare_damaging_case + rare_damaging_control)) -v den="$dnv_var_count" \
        'BEGIN{ if(den==0) {print "nan"} else {printf "%.6f", num/den} }')
elif [[ "$rare_damaging_mode" == "control" ]]; then
    echo "Calculating rare damaging variant rate using only control DNVs"
    if [[ -n "$cohort_filter" ]]; then
        dnv_var_count=$(bcftools view -i "($cohort_filter) && INFO/STATUS=\"control\"" -H "$dnv_vcf" | wc -l)
    else
        dnv_var_count=$(bcftools view -i 'INFO/STATUS="control"' -H "$dnv_vcf" | wc -l)
    fi
    rare_damaging_rate=$(awk -v num="$rare_damaging_control" -v den="$dnv_var_count" \
        'BEGIN{ if(den==0) {print "nan"} else {printf "%.6f", num/den} }')
fi

echo "Rare damaging variant rate in DNV control: $rare_damaging_rate"
# providing mut rate files (Samocha et al. baseline rates)
# python src/2_TADA/generate_TADA_input.py \
#     --outdir "$outdir" \
#     --config "$config" \
#     merge \
#     --mut-rate-file "$mutation_rate_file"

# using Roulette mutation rates
roulette_raw_dir="/fs/cbsuhy02/storage/yz2676/ASD/NetFlow_ASD/Yutong_scripts/data/supplemental/roulette"
roulette_generation_factor=1.015e-7
python src/2_TADA/generate_TADA_input.py \
    --outdir "$outdir" \
    --config "$config" \
    merge \
    --cohort "$samples"

mapfile -t keys_sanitized < <(
  jq -r 'keys[] | gsub(" "; "_")' $config
)
for key in "${keys_sanitized[@]}"; do
    regions_file="${outdir}/${samples}/input/TADA_${key}.txt"
    out="${outdir}/${samples}/input/TADA_${key}_roulette_rd_${rare_damaging_mode}.txt"
    echo "Processing region: $key"

    : > "$out".tmp
    tmp_bed="$(mktemp)"
    # Skip header and normalize chr prefix: chr1->1, chrX->X (if present)
    awk -v OFS="\t" 'NR==1{next} {sub(/^chr/,"",$1); print}' "$regions_file" > "$tmp_bed"

    # Generate mutation rates sums using bcftools and bedtools
    for chr in {1..22}; do
        vcf="${roulette_raw_dir}/${chr}_rate_v5.2_TFBS_correction_all.vcf.bgz"

        awk -v c="$chr" '$1==c' "${tmp_bed}" > "${regions_file}.chr${chr}"
        
        [[ ! -s "${regions_file}.chr${chr}" ]] && continue

        bcftools query -R "${regions_file}.chr${chr}" -f '%CHROM\t%POS\t%INFO/MR\t%INFO/MG\n' "$vcf" \
            | awk -v OFS="\t" '{print $1, $2-1, $2, $3, $4}' \
            | sort -k1,1 -k2,2n \
            | bedtools map -a "${regions_file}.chr${chr}" -b - -c 4,5,4,5 -o sum,sum,count,count -null "." \
            | awk -v OFS="\t" -v s="$roulette_generation_factor" -v rd="$rare_damaging_rate" -v chrom="$chr" \
            'BEGIN {miss_mr = miss_mg = 0}
            {
                len = 3*($3 - $2)

                mr_sum = $(NF-3); mg_sum = $(NF-2)
                mr_cnt = $(NF-1); mg_cnt = $NF
                
                if (mr_cnt != len) {mr_out = "."; miss_mr++} else {mr_out = mr_sum * s * rd}
                if (mg_cnt != len) {mg_out = "."; miss_mg++} else {mg_out = mg_sum * s * rd}

                print $4, $5, $6, $7, $8, $9, mr_out, mg_out;
            }
            END {
                printf "SUMMARY: chrom=%s\tmissing_MR=%d\tmissing_MG=%d\n", chrom, miss_mr, miss_mg > "/dev/stderr"
            }
            ' >> "$out".tmp
        rm -f "${regions_file}.chr${chr}"
    done
    rm -f "$tmp_bed"

    # Add header
    header=$(head -n 1 "$regions_file" | cut -f4-)
    echo -e "${header}\tRoulette_MR\tRoulette_MG" > "$out"
    cat "$out".tmp >> "$out"
    rm -f "$regions_file" "$out".tmp
done

