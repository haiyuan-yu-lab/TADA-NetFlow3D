#!/bin/bash
#SBATCH --job-name=generate_UCSF_coding_table
#SBATCH -o logs/generate_UCSF_coding_table.%A.out
#SBATCH -e logs/generate_UCSF_coding_table.%A.err
#SBATCH --mail-type=ALL
#SBATCH --mail-user=yz2676@cornell.edu

#SBATCH -N 1
#SBATCH --ntasks-per-node=1
#SBATCH --mem=32G
#SBATCH --time=48:00:00
#SBATCH --nodelist=cbsuhy02

COHORT="InhAll_DNVAll" # InhMult_DNVSimp / InhAll_DNVAll
while [[ $# -gt 0 ]]; do
  case "$1" in
    --cohort)               COHORT="$2"; shift 2;;
    *) echo "Unknown arg: $1" >&2; exit 1;;
  esac
done

outbase=/fs/cbsuhy02/storage/yz2676/ASD/NetFlow_ASD/Yutong_scripts/data/TADA/UCSF_coding/roulette
outdir="${outbase}/${COHORT}"
mkdir -p "$outdir"

lof_regex="frameshift_variant\|splice_acceptor_variant\|splice_donor_variant\|stop_gained\|stop_lost\|start_lost\|inframe_deletion\|inframe_insertion\|protein_altering_variant"
var_prefilter_dir="/fs/cbsuhy01/storage/yz2676/projects/1-ASD/NetFlow_ASD/Yutong_scripts/data/asdinfo"

prefix_dnv="UCSF_20211001_dnms_WGS_hg38.codingonly.filtered.vep.with_asdinfo"
prefix_dnv_all="UCSF_20211001_dnms_WGS_hg38.includemultiplex.codingonly.filtered.vep.with_asdinfo"
prefix_inh_mult="All_inherited_WGS_hg38.codingonly.filtered.vep.with_asdinfo"
prefix_inh_all="All_inherited_WGS_hg38.includesimplex.codingonly.filtered.vep.with_asdinfo"
prefix_untrans_mult="UCSF_parents_WGS_hg38.onlymultiplex.codingonly.untransmitted.dedup"
prefix_untrans_all="UCSF_parents_WGS_hg38.includesimplex.codingonly.untransmitted.dedup"

# de novo
dnv_simp_vcf="${var_prefilter_dir}/${prefix_dnv}.bcf"
dnv_all_vcf="${var_prefilter_dir}/${prefix_dnv_all}.bcf"

# transmitted
inh_mult_vcf="${var_prefilter_dir}/${prefix_inh_mult}.bcf"
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

echo "Generating TADA input files in $outdir for cohort $COHORT"
echo "DNV VCF: $dnv_vcf"
echo "Inherited VCF: $inh_vcf"
echo "Untransmitted VCF: $untrans_vcf"

python src/2_TADA/generate_TADA_input.py \
    --outdir "$outdir" \
    --config "" \
    build-coding \
    --vcf "$dnv_vcf" \
    --type "dnv" \
    --lof-regex "$lof_regex" \
    --gnomad-null-only

python src/2_TADA/generate_TADA_input.py \
    --outdir "$outdir" \
    --config "" \
    build-coding \
    --vcf "$dnv_vcf" \
    --type "dnv" \
    --misb \
    --gnomad-null-only

python src/2_TADA/generate_TADA_input.py \
    --outdir "$outdir" \
    --config "" \
    build-coding \
    --vcf "$dnv_vcf" \
    --type "dnv" \
    --misa \
    --gnomad-null-only

python src/2_TADA/generate_TADA_input.py \
    --outdir "$outdir" \
    --config "" \
    build-coding \
    --vcf "$inh_vcf" \
    --type "inh" \
    --lof-regex "$lof_regex" \
    --gnomad-null-only

python src/2_TADA/generate_TADA_input.py \
    --outdir "$outdir" \
    --config "" \
    build-coding \
    --vcf "$inh_vcf" \
    --type "inh" \
    --misb \
    --gnomad-null-only

python src/2_TADA/generate_TADA_input.py \
    --outdir "$outdir" \
    --config "" \
    build-coding \
    --vcf "$inh_vcf" \
    --type "inh" \
    --misa \
    --gnomad-null-only

python src/2_TADA/generate_TADA_input.py \
    --outdir "$outdir" \
    --config "" \
    build-coding \
    --vcf "$untrans_vcf" \
    --type "untrans" \
    --lof-regex "$lof_regex" \
    --gnomad-null-only

python src/2_TADA/generate_TADA_input.py \
    --outdir "$outdir" \
    --config "" \
    build-coding \
    --vcf "$untrans_vcf" \
    --type "untrans" \
    --misb \
    --gnomad-null-only

python src/2_TADA/generate_TADA_input.py \
    --outdir "$outdir" \
    --config "" \
    build-coding \
    --vcf "$untrans_vcf" \
    --type "untrans" \
    --misa \
    --gnomad-null-only

echo "Merging coding TADA input files..."
python src/2_TADA/generate_TADA_input.py \
    --outdir "$outdir" \
    --config "" \
    merge-coding \
    --gencode-file "/fs/cbsuhy01/storage/yz2676/projects/1-ASD/NetFlow_ASD/Yutong_scripts/data/enhancer/gencode_v43_coding/gencode.v43.annotation.genes.bed" \
    --constraint-file "/fs/cbsuhy02/storage/yz2676/ASD/NetFlow_ASD/Yutong_scripts/data/supplemental/gnomad_v4_constraint/gnomad.v4.1.constraint_metrics.tsv"

echo "Calculating the mutation rates (ptv/misb/misa) for each gene ..."
ptv_control=$(bcftools view -H "$outdir/dnv/ptv/ctrl.vcf.gz" | wc -l)
misa_control=$(bcftools view -H "$outdir/dnv/misa/ctrl.vcf.gz" | wc -l)
misb_control=$(bcftools view -H "$outdir/dnv/misb/ctrl.vcf.gz" | wc -l)

dnv_var_count=$(bcftools view -i 'INFO/STATUS="control"' -H "$dnv_vcf" | wc -l)

ptv_mut_rate=$(awk -v num="$ptv_control" -v den="$dnv_var_count" \
  'BEGIN{ if(den==0) {print "nan"} else {printf "%.6f", num/den} }')
misa_mut_rate=$(awk -v num="$misa_control" -v den="$dnv_var_count" \
  'BEGIN{ if(den==0) {print "nan"} else {printf "%.6f", num/den} }')
misb_mut_rate=$(awk -v num="$misb_control" -v den="$dnv_var_count" \
  'BEGIN{ if(den==0) {print "nan"} else {printf "%.6f", num/den} }')

echo "PTV mutation rate: $ptv_mut_rate per generation"
echo "MISB mutation rate: $misb_mut_rate per generation"
echo "MISA mutation rate: $misa_mut_rate per generation"


# using Roulette mutation rates
roulette_raw_dir="/fs/cbsuhy02/storage/yz2676/ASD/NetFlow_ASD/Yutong_scripts/data/supplemental/roulette"
roulette_generation_factor=1.015e-7

# all exon file to estimate mut rate
regions_file="/fs/cbsuhy01/storage/yz2676/projects/1-ASD/NetFlow_ASD/Yutong_scripts/data/enhancer/gencode_v43_coding/gencode.v43.annotation.exons.bed"
out="${outdir}/input_coding/exons_roulette.txt"
echo "Processing region: coding regions file: $regions_file"


awk -F'\t' 'BEGIN{OFS="\t"}
{
    key = $4 OFS $5
    seen[key] = 1

    if ($6=="." || $7=="." || $8=="." || $9=="." || $10=="." || $11==".") {
        bad[key] = 1
    } else {
        sum6[key]  += $6
        sum7[key]  += $7
        sum8[key]  += $8
        sum9[key]  += $9
        sum10[key] += $10
        sum11[key] += $11
    }
}
END {
    print "gene_name","Hugo_Symbol","Roulette_MR_ptv","Roulette_MG_ptv","Roulette_MR_misb","Roulette_MG_misb","Roulette_MR_misa","Roulette_MG_misa"
    for (k in seen) {
        split(k, a, OFS)
        if (!bad[k]) {
            print a[1], a[2], sum6[k], sum7[k], sum8[k], sum9[k], sum10[k], sum11[k]
        }
        else {
            print a[1], a[2], ".", ".", ".", ".", ".", "."
        }
    }
}' "$out" > "${outdir}/input_coding/TADA_coding_roulette_genelevel.txt"

awk -F'\t' '
{
    key = $4 FS $5
    seen[key] = 1
    if ($6=="." || $7=="." || $8=="." || $9=="." || $10=="." || $11==".") {
        bad[key] = 1
    }
}
END {
    n = 0
    for (k in seen) {
        if (!(k in bad)) n++
    }
    print "Total genes with non-missing Roulette mutation rates:", n
}' "$out"

coding_file="${outdir}/input_coding/TADA_coding.txt"
roulette_file="${outdir}/input_coding/TADA_coding_roulette_genelevel.txt"
merged_file="${outdir}/input_coding/TADA_coding_roulette.txt"

awk -F'\t' 'BEGIN{OFS="\t"}
FNR==NR {
    if (FNR==1) next
    symbol = $1
    gid = $2
    sub(/\.[0-9]+$/, "", gid)
    rates = $3 OFS $4 OFS $5 OFS $6 OFS $7 OFS $8
    if (!(gid in rid)) rid[gid] = rates
    if (!(symbol in rsym)) rsym[symbol] = rates
    next
}
FNR==1 {
    print $0, "Roulette_MR_ptv","Roulette_MG_ptv","Roulette_MR_misb","Roulette_MG_misb","Roulette_MR_misa","Roulette_MG_misa"
    next
}
{
    gid = $1
    symbol = $2
    sub(/\.[0-9]+$/, "", gid)

    if (gid in rid) {
        print $0, rid[gid]
    } else if (symbol in rsym) {
        print $0, rsym[symbol]
    } else {
        print $0, ".", ".", ".", ".", ".", "."
    }
}' "$roulette_file" "$coding_file" > "$merged_file"