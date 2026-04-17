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

outbase=/fs/cbsuhy02/storage/yz2676/ASD/NetFlow_ASD/Yutong_scripts/data/TADA/UCSF_coding/gnomad
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


echo "Calculating the avg misb/misa mutation proportions ..."
ptv_control=$(bcftools view -H "$outdir/dnv/ptv/ctrl.vcf.gz" | wc -l)
misa_control=$(bcftools view -H "$outdir/dnv/misa/ctrl.vcf.gz" | wc -l)
misb_control=$(bcftools view -H "$outdir/dnv/misb/ctrl.vcf.gz" | wc -l)
echo "Total control PTV variants: $ptv_control"
echo "Total control MISA variants: $misa_control"
echo "Total control MISB variants: $misb_control"

dnv_missense=$(
  bcftools view -i 'INFO/STATUS="control"' -Ou "$dnv_vcf" \
  | bcftools +split-vep -c Consequence -i 'Consequence~"missense_variant"' -Ou - \
  | bcftools query -f '%CHROM\t%POS\t%REF\t%ALT\t%INFO/SAMPLE\n' \
  | wc -l
)
echo "Total DNV missense variants: $dnv_missense"

misa_mut_rate=$(awk -v num="$misa_control" -v den="$dnv_missense" \
  'BEGIN{ if(den==0) {print "nan"} else {printf "%.6f", num/den} }')
misb_mut_rate=$(awk -v num="$misb_control" -v den="$dnv_missense" \
  'BEGIN{ if(den==0) {print "nan"} else {printf "%.6f", num/den} }')

echo "MISB mutation rate among missense variants: $misb_mut_rate"
echo "MISA mutation rate among missense variants: $misa_mut_rate"

awk -v OFS='\t' -v MISB="$misb_mut_rate" -v MISA="$misa_mut_rate" '
NR==1 {print $0, "misb.mu", "misa.mu"}
NR>1 {print $0, MISB*$16, MISA*$16}
' "${outdir}/input_coding/TADA_coding.txt" > "${outdir}/input_coding/TADA_coding.mut.txt"

# scale to match the number of control variants
ptv_col=15
misb_col=18
misa_col=19
n_control_trios=$(bcftools query -f '%INFO/SAMPLE\t%INFO/STATUS\n' "$dnv_vcf" \
  | awk '$2=="control"{print $1}' | sort -u | wc -l)

sum_misb=$(awk -F'\t' -v c="$misb_col" 'NR>1 && $c!="nan" {s+=$c} END{printf "%.6e", s+0}' "${outdir}/input_coding/TADA_coding.mut.txt")
expected_misb_control=$(awk -v N="$n_control_trios" -v S="$sum_misb" 'BEGIN{printf "%.6f", N*S}')

sum_misa=$(awk -F'\t' -v c="$misa_col" 'NR>1 && $c!="nan" {s+=$c} END{printf "%.6e", s+0}' "${outdir}/input_coding/TADA_coding.mut.txt")
expected_misa_control=$(awk -v N="$n_control_trios" -v S="$sum_misa" 'BEGIN{printf "%.6f", N*S}')

sum_ptv=$(awk -F'\t' -v c="$ptv_col" 'NR>1 && $c!="nan" {s+=$c} END{printf "%.6e", s+0}' "${outdir}/input_coding/TADA_coding.mut.txt")
expected_ptv_control=$(awk -v N="$n_control_trios" -v S="$sum_ptv" 'BEGIN{printf "%.6f", N*S}')

scaler_ptv=$(awk -v expc="$expected_ptv_control" -v obs="$ptv_control" \
  'BEGIN{ if(expc==0) {print "1.0"} else {printf "%.6f", obs/expc} }')
scaler_misb=$(awk -v expc="$expected_misb_control" -v obs="$misb_control" \
  'BEGIN{ if(expc==0) {print "1.0"} else {printf "%.6f", obs/expc} }')
scaler_misa=$(awk -v expc="$expected_misa_control" -v obs="$misa_control" \
  'BEGIN{ if(expc==0) {print "1.0"} else {printf "%.6f", obs/expc} }')
echo "PTV scaler: $scaler_ptv"
echo "MISB scaler: $scaler_misb"
echo "MISA scaler: $scaler_misa"

awk -v OFS='\t' -v PTV="$scaler_ptv" -v MISB="$scaler_misb" -v MISA="$scaler_misa" '
NR==1 {print $0, "ptv.mu.scaled", "misb.mu.scaled", "misa.mu.scaled"}
NR>1 {print $0, $15*PTV, $18*MISB, $19*MISA}
' "${outdir}/input_coding/TADA_coding.mut.txt" > "${outdir}/input_coding/TADA_coding.mut.scaled.txt"
