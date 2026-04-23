#!/bin/bash
#SBATCH --job-name=run_TADA
#SBATCH -o logs/run_TADA.%A.out
#SBATCH -e logs/run_TADA.%A.err
#SBATCH --mail-type=ALL
#SBATCH --mail-user=yz2676@cornell.edu

#SBATCH -N 1
#SBATCH --ntasks-per-node=1
#SBATCH --mem=8G
#SBATCH --time=48:00:00
#SBATCH --nodelist=cbsuhy02

MODE="coding+noncoding"          # coding | coding+noncoding
VAR_SIGNAL="both"       # dn | in | both

# ---- decide analysis flags + naming from MODE/VAR_SIGNAL ----
base_dir=/fs/cbsuhy02/storage/yz2676/ASD/NetFlow_ASD
TADA_dir=${base_dir}/Yutong_scripts/data/TADA
samples=("InhAll_DNVAll") # InhMult_DNVSimp
enhancer_set=("Deconv Excitatory Neurons" "Deconv Astrocytes" "Deconv Microglia" "Combined" "Combined Simple")

coding="FuScaled" # FuScaled / FuRaw
mut="Roulette_MR"
cohort="SSC_KOR_STAR" # SSC_KOR_STAR / SPARK / ALL
rare_damaging_mode="all" # all / control

coding_file=${base_dir}/TADA_input_coding_${coding}_${cohort}.txt
prior_file=${base_dir}/prior.txt
info_file=${base_dir}/info.txt

r_flags=( --dn_coding --in_coding )      # always do coding here

while [[ $# -gt 0 ]]; do
  case "$1" in
    --mode)               MODE="$2"; shift 2;;
    --coding)             coding="$2"; shift 2;;
    --mut)                mut="$2"; shift 2;;
    --cohort)             cohort="$2"; shift 2;;
    --rare-damaging-mode) rare_damaging_mode="$2"; shift 2;;
    *) echo "Unknown arg: $1" >&2; exit 1;;
  esac
done

if [[ "$MODE" == "coding+noncoding" ]]; then
    case "$VAR_SIGNAL" in
        dn)   r_flags+=( --dn_noncoding ) ;;
        in)   r_flags+=( --in_noncoding ) ;;
        both) r_flags+=( --dn_noncoding --in_noncoding ) ;;
        *)    echo "ERROR: VAR_SIGNAL must be dn|in|both, got: $VAR_SIGNAL" >&2; exit 1 ;;
    esac
elif [[ "$MODE" != "coding" ]]; then
    echo "ERROR: MODE must be coding|coding+noncoding, got: $MODE" >&2
    exit 1
fi

# calculate the number of individuals once
# bcftools query -f '%INFO/SAMPLE\t%INFO/STATUS\n' UCSF_20211001_dnms_WGS_hg38.with_multiplex.filtered.vep.with_asdinfo.bcf | awk '$2=="case"{print $1}' | sort -u | wc -l
# bcftools query -f '%INFO/SAMPLE\t%INFO/STATUS\n' All_inherited_WGS_hg38.includesimplex.filtered.vep.with_asdinfo.bcf | awk '$2=="control"{print $1}' | sort -u | wc -l

jobs=()
if [[ "$MODE" == "coding" ]]; then
    jobs+=( "codingOnly|" )   # no noncoding file
else
    for enhancers in "${enhancer_set[@]}"; do
        noncoding="${enhancers// /_}"
        input_file=TADA_${noncoding}_roulette_rd_${rare_damaging_mode}.txt
        jobs+=( "${noncoding}|${input_file}" )
    done
fi

# ---- run TADA ----
for sample in "${samples[@]}"; do
    sample_file=${base_dir}/sample_${coding}_${sample}.txt

    for job in "${jobs[@]}"; do
        IFS='|' read -r tag input_file <<< "$job"
        noncoding_file=${TADA_dir}/${sample}/input/${input_file}

        if [[ "$MODE" == "coding" ]]; then
        output_prefix=${TADA_dir}/${sample}/TADA_results/TADA_${coding}_${tag}_${mut}
        noncoding_args=()
        else
        output_prefix=${TADA_dir}/${sample}/TADA_results/TADA_${coding}_${tag}_${mut}
        noncoding_args=( --noncoding_file "${noncoding_file}" )
        fi

        mkdir -p "$(dirname ${output_prefix})"

        echo "Running: sample=${sample} tag=${tag} MODE=${MODE} mut=${mut} VAR_SIGNAL=${VAR_SIGNAL}"

        Rscript src/2_TADA/run_TADA.r \
            --coding_file ${coding_file} \
            --prior_file ${prior_file} \
            --info_file ${info_file} \
            --sample_file ${sample_file} \
            --output_file "${output_prefix}.txt" \
            --mut_col ${mut} \
            "${r_flags[@]}" \
            "${noncoding_args[@]}"
        
        # fdr_threshold=0.05
        # filtered_output=${output_prefix}_filtered_f${fdr_threshold}.txt
        # awk -F "\t" -v fdr_th=${fdr_threshold} '
        # NR==1 {for(i=1;i<=NF;i++) if($i=="fdr") fdr_col=i; print; next}
        # fdr_col && $fdr_col+0 < fdr_th {print}
        # ' "${output_prefix}.txt" > "${filtered_output}"
        # echo "Filtered results saved to ${filtered_output}"
    done
done