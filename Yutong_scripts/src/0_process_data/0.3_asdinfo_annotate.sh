#!/bin/bash
#SBATCH --job-name=annotate_asdinfo
#SBATCH -o logs/annotate_asdinfo.%j.out
#SBATCH -e logs/annotate_asdinfo.%j.err
#SBATCH --mail-type=ALL
#SBATCH --mail-user=yz2676@cornell.edu

#SBATCH -N 1
#SBATCH --ntasks-per-node=1
#SBATCH --mem=8G
#SBATCH --time=48:00:00
#SBATCH --nodelist=cbsuhy02


# Final step of data preprocessing:
# 1) Exclude siblings with ASD risk
# 2) Exclude multiplex individuals for DNV analysis
# 3) Exclude LOF & damaging missense variants:
#    a) LOF: - frameshift, splicing, stopgain, stoploss, startloss
#    b) Damaging missense: - missense consequence with REVEL >= 0.5
# 4) Annotate necessary metadata (sex, etc.)

# == File paths ==
vep_anno_dir=/fs/cbsuhy01/storage/yz2676/projects/1-ASD/NetFlow_ASD/Yutong_scripts/data/annotated/
dnv_vep=$vep_anno_dir/dnm/UCSF_20211001_dnms_WGS_hg38.with_vep.vcf.gz
inh_vep=$vep_anno_dir/inh/All_inherited_WGS_hg38.with_vep.vcf.gz

pheno=/fs/cbsuhy01/storage/NetFlow_ASD/Week8_429/SSC+SPARK+KOR-STAR_parsed.txt
asd_family=/fs/cbsuhy01/storage/NetFlow_ASD/Week8_429/roles_index.csv
concern_individual=/fs/cbsuhy01/storage/NetFlow_ASD/Week8_429/individuals_registration.csv

# == Parameters ==
FAMILY_TAG=FAM
SAMPLE_TAG=SAMPLE
LOF="(frameshift_variant|splice_acceptor_variant|splice_donor_variant|stop_gained|stop_lost|start_lost)"
DAMAGING_MISSENSE="missense_variant"

# == Build exclusion lists ==
asd_dir=data/asdinfo
mkdir -p "$asd_dir"

# 1) Multiplex families
fam_list="$asd_dir/multiplex_families.txt"
csvgrep -c family_type -m Multiplex "$asd_family" \
| csvcut -c family_sf_id \
| tail -n +2 | sort | uniq \
> "$fam_list"


# 2) Individuals with ASD risk (based on concern reason)
#   - Siblings with potential developmental concerns
#   - Siblings who are too young for a reliable ASD diagnosis
samp_list="$asd_dir/samples_exclude.txt"

# Concern siblings (non-empty control_sibling_flag)
csvgrep -c control_sibling_flag -m '1' "$concern_individual" \
| csvcut -c subject_sp_id \
| tail -n +2 | sort | uniq \
> "$asd_dir/concern_siblings.txt"

# Young siblings from Unknown families:
# Father!="0" & Mother!="0" & Pheno=="control" & Fam in unknown_fam
csvgrep -c family_type -m Unknown "$asd_family" \
| csvcut -c family_sf_id \
| awk '$0!="SF0153305"' | sort | uniq \
> "$asd_dir/unknown_families.txt"

csvgrep -t -c 3 -i -r '^0$' "$pheno" \
| csvgrep -c 4 -i -r '^0$' \
| csvgrep -c 6 -r '^control$' \
> "$asd_dir/control_siblings.txt"

csvjoin -c Fam,family_sf_id \
  "$asd_dir/control_siblings.txt" "$asd_dir/unknown_families.txt" \
| csvcut -c SampleID \
| tail -n +2 | awk 'NF' | sort -u > "$asd_dir/young_sib.txt"

cat "$asd_dir/young_sib.txt" "$asd_dir/concern_siblings.txt" \
| sort -u > "$samp_list"


# == Filter the VCF/BCF ==
# Convert once to BCF for speed
tmp_vep="$asd_dir/UCSF_20211001_dnms_WGS_hg38.bcf"
bcftools view --threads 8 -Ob -o "$tmp_vep" "$dnv_vep"
bcftools index -f "$tmp_vep"

# Build a position list for filtering
cache="$asd_dir/cache_fam_sample.bed.gz"
bcftools query -f '%CHROM\t%POS\t%POS\t%INFO/FAM\t%INFO/SAMPLE\n' "$tmp_vep" \
| bgzip > "$cache"
tabix -s1 -b2 -e3 "$cache"

# Exclude if:
#  - INFO/FAM has a Multiplex family, OR
#  - INFO/SAMPLE is a risky sibling
filtered_asd="$asd_dir/UCSF_20211001_dnms_WGS_hg38.filtered.bcf"
bcftools view -T ^<(
  {
    # A) positions where FAM ∈ list
    zcat "$cache" \
    | awk -F'\t' 'NR==FNR{fam[$1]=1; next} { # load list into fam[]
        if($4 in fam) print $1"\t"$2"\t"$3
        }' "$fam_list" -

    # B) positions where SampleID ∈ list
    zcat "$cache" \
    | awk -F'\t' 'NR==FNR{ok[$1]=1; next}{
        if($5 in ok) print $1"\t"$2"\t"$3
        }' "$samp_list" -
  } | sort -u
) -Ob -o "$filtered_asd" "$asd_dir/UCSF_20211001_dnms_WGS_hg38.bcf"
bcftools index -f "$filtered_asd"

# Exclude if: (optional)
#  - LOF or damaging missense
noncoding_asd="$asd_dir/UCSF_20211001_dnms_WGS_hg38.filtered.noncoding.bcf"
# if any of the transcript consequence is not LOF or damaging missense, keep the variant
bcftools +split-vep \
    -c REVEL:Float \
    -i '(Consequence!~"'$LOF'") && (Consequence!~"missense_variant" || REVEL<=0.5)' \
    -Ou "$filtered_asd" \
| bcftools annotate -x 'INFO/Consequence,INFO/REVEL' \
    -Ob -o "$noncoding_asd"


# Sanity check: count variants
N_IN=$(bcftools view -H "$tmp_vep" | wc -l)
echo "Input variants: $N_IN"

N_OUT=$(bcftools view -H "$filtered_asd" | wc -l 2>/dev/null || echo 0)
echo "Filtered (Multiplex + concerning individuals) variants: $N_OUT"

N_FILTERED=$(bcftools view -H "$noncoding_asd" | wc -l 2>/dev/null || echo 0)
echo "Filtered (LOF + damaging missense) variants: $N_FILTERED"

# Annotate
final_asd="$asd_dir/UCSF_20211001_dnms_WGS_hg38.filtered.noncoding.with_asdinfo.bcf"
sex_map="$asd_dir/sex.map"
awk -F'\t' 'NR==1{
  for(i=1;i<=NF;i++){
    if($i=="SampleID") sid=i
    if(tolower($i)=="sex") sx=i
  }
  next
}
{print $sid "\t" $sx}' "$pheno" > "$sex_map"

{
bcftools view -h "$noncoding_asd" \
| awk '
      BEGIN { add=1 }
      /^##INFO=/ && $0 ~ /ID=SEX,/ { add=0 }
      /^#CHROM/ {                                     
        if (add)
          print "##INFO=<ID=SEX,Number=1,Type=String,Description=\"Sample sex from mapping\">"
        print; next
      }
      { print } 
  '
bcftools view -H "$noncoding_asd" \
  | awk -v FS="\t" -v OFS="\t" -v MAP="$sex_map" -v KEY="SAMPLE" '
      BEGIN{
      if (MAP=="") { print "ERROR: MAP is empty" > "/dev/stderr"; exit 1 }
        # Load SampleID -> sex
        while((getline < MAP) >= 0){
          sub(/\r$/,"");
          split($0,f,"\t"); if(f[1] && f[2]) { sex[f[1]]=f[2]; nmap++ }
        }
        close(MAP)
        print "INFO: loaded", nmap, "sample->sex pairs" > "/dev/stderr"
        changed=0; missing=0; noKey=0
      }
      function get(info,key, a){
        if (match(info, "(^|;)" key "=([^;]*)", a)) return a[2]
        return ""
      }
      {
        info = $8
        sid = get(info, KEY)
        if (sid=="") { noKey++; print; next }
        if(sid in sex){
          gsub("(^|;)SEX=[^;]*","",$8)
          $8 = info";SEX=" sex[sid]
        changed++
        } else {
          # SAMPLE not in map
          missing++
        }
        print
      }
      END{
        print "INFO: rows edited:", changed,
              " | rows with SAMPLE not in map:", missing,
              " | rows missing", KEY ":", noKey > "/dev/stderr"
      }'
} | bcftools view -Ob -o "$final_asd" -

bcftools index -f "$final_asd"

# Remove temporary files
rm "$asd_dir/concern_siblings.txt" "$asd_dir/unknown_families.txt" "$asd_dir/control_siblings.txt" "$asd_dir/young_sib.txt"
rm "$tmp_vep" "$tmp_vep".csi
