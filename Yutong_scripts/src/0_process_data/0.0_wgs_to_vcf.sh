#!/bin/bash
#SBATCH --job-name=wgs_to_vcf.parents
#SBATCH -o logs/wgs_to_vcf.parents.%j.out
#SBATCH -e logs/wgs_to_vcf.parents.%j.err
#SBATCH --mail-type=ALL
#SBATCH --mail-user=yz2676@cornell.edu

#SBATCH -N 1
#SBATCH --ntasks-per-node=1
#SBATCH --mem=64G
#SBATCH --time=48:00:00
#SBATCH --nodelist=cbsuhy02

# Transform WGS data (DNV, Inherited and Parents) to VCF format
# Data cleaning: keep only primary contigs

## Raw files
wgs_var_dir=/fs/cbsuhy01/storage/NetFlow_ASD/Week8_429/UCSF_raw

dnv_var=UCSF_20211001_dnms_WGS_hg38.tsv
inherited_var=("KOR-STAR" "SPARK_WGS_1" "SPARK_WGS_2" "SPARK_WGS_3" "SSC_WGS9209s")
parents_var=/fs/cbsuhy01/storage/NetFlow_ASD/Week8_429/Variants/UCSF_Parent.txt.gz

## Reference
ref=/fs/cbsuhy01/storage/yz2676/data/DetailedProCapNet/wdir/reference/hg38.22XY.fasta

## Outputs
data_dir=data/raw_processed
mkdir -p $data_dir

dnv_var_vcf=$data_dir/UCSF_20211001_dnms_WGS_hg38.vcf

### Process DNV first
echo "Processing DNV"
#### 1. Minimal header
{
  printf '##fileformat=VCFv4.2\n'
  printf '##INFO=<ID=TYPE,Number=1,Type=String,Description="Variant type">\n'
  printf '##INFO=<ID=LEN,Number=1,Type=Integer,Description="Variant length">\n'
  printf '##INFO=<ID=SAMPLE,Number=1,Type=String,Description="Sample_ID">\n'
  printf '##INFO=<ID=FAM,Number=1,Type=String,Description="Family_ID">\n'
  printf '##INFO=<ID=STATUS,Number=1,Type=String,Description="Status">\n'
  printf '##INFO=<ID=BATCH,Number=1,Type=String,Description="Batch_ID">\n'
  printf '#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\n'
} > "$data_dir/dnv_header.vcf"

#### 2. TSV -> VCF body
awk -F'\t' -v OFS='\t' '
NR==1 {
    for(i=1;i<=NF;i++){ k=$i; sub(/^#/,"",k); H[k]=i } next
}
{
    chrom = $(H["CHROM"]); sub(/^chr/,"",chrom)
    info =  "TYPE=" $(H["TYPE"]) \
            ";LEN=" $(H["LEN"]) \
            ";SAMPLE=" $(H["SAMPLE"]) \
            ";FAM=" $(H["FAM"]) \
            ";STATUS=" $(H["STATUS"]) \
            ";BATCH=" $(H["BATCH"])
  print chrom, $(H["POS"]), ".", $(H["REF"]), $(H["ALT"]), ".", "PASS", info
}' $wgs_var_dir/$dnv_var > $data_dir/dnv_body.vcf

#### 3. Combine header and body
cat $data_dir/dnv_header.vcf $data_dir/dnv_body.vcf > $dnv_var_vcf

awk -v OFS='' '
{
    id=$1
    sub(/^chr/,"",id)
    print "##contig=<ID="id",length="$2">"
}' "$ref.fai" > $data_dir/contigs.hdr

bcftools annotate -h $data_dir/contigs.hdr -Oz -o $dnv_var_vcf.unsorted.gz $dnv_var_vcf

bcftools sort -T $data_dir/tmp_sort \
  -Oz -o $dnv_var_vcf.gz \
  $dnv_var_vcf.unsorted.gz

#### 4. Index VCF
tabix -p vcf $dnv_var_vcf.gz

#### 5. Clean up
rm $data_dir/dnv_header.vcf $data_dir/dnv_body.vcf $dnv_var_vcf.unsorted.gz $dnv_var_vcf


### Same steps for inherited variants
kor_inherited_vcf=$data_dir/KOR-STAR_inherited_WGS_hg38.vcf
spark_inherited_vcf=$data_dir/SPARK_inherited_WGS_hg38.vcf
ssc_inherited_vcf=$data_dir/SSC_inherited_WGS_hg38.vcf
inh_var_vcf=$data_dir/All_inherited_WGS_hg38.vcf

# -----------------------------
# Build header
# -----------------------------
inh_header="$data_dir/inh_header.vcf"
{
  printf '##fileformat=VCFv4.2\n'
  printf '##FILTER=<ID=INDEL_TS90_99,Description="GATK VQSR INDEL tranche 90.00–99.00">\n'
  printf '##FILTER=<ID=INDEL_TS99_999,Description="GATK VQSR INDEL tranche 99.00–99.90">\n'
  printf '##FILTER=<ID=INDEL_TS999_100,Description="GATK VQSR INDEL tranche 99.90–100.00">\n'
  printf '##FILTER=<ID=INDEL_TS99_100,Description="GATK VQSR INDEL tranche 99.00–100.00">\n'
  printf '##INFO=<ID=SAMPLE,Number=1,Type=String,Description="Sample_ID">\n'
  printf '##INFO=<ID=FAM,Number=1,Type=String,Description="Family_ID">\n'
  printf '##INFO=<ID=FATHER,Number=1,Type=String,Description="Father_ID">\n'
  printf '##INFO=<ID=MOTHER,Number=1,Type=String,Description="Mother_ID">\n'
  printf '##INFO=<ID=CHILD_GT,Number=1,Type=String,Description="Child genotype">\n'
  printf '##INFO=<ID=FATHER_GT,Number=1,Type=String,Description="Father genotype">\n'
  printf '##INFO=<ID=MOTHER_GT,Number=1,Type=String,Description="Mother genotype">\n'
  printf '##INFO=<ID=NUM_ALLELES,Number=1,Type=Integer,Description="Number of alleles">\n'
  printf '##INFO=<ID=BATCH,Number=1,Type=String,Description="Batch_ID">\n'
  printf '#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\n'
} > "$inh_header"

# Build a regex set of allowed (primary) contigs for quick filtering
primary_set_regex='^([1-9]|1[0-9]|2[0-2]|X|Y)$'

# -----------------------------
# Helper function: tranche -> FILTER mapping
# -----------------------------
map_filter='
function map_filter(fr) {
    if (fr=="" || fr=="." || fr=="[]") return "PASS";
    gsub(/"/, "", fr); gsub(/^\[|\]$/, "", fr); gsub(/"/, "", fr);
    if      (fr=="VQSRTrancheINDEL90.00to99.00")  y="INDEL_TS90_99";
    else if (fr=="VQSRTrancheINDEL99.00to99.90")  y="INDEL_TS99_999";
    else if (fr=="VQSRTrancheINDEL99.90to100.00") y="INDEL_TS999_100";
    else if (fr=="VQSRTrancheINDEL99.00to100.00") y="INDEL_TS99_100";
    else y=fr;
    return y
}
'

# -----------------------------
# Function: TSV(.gz) -> primary-only VCF body (no header), with normalized CHROM names
# -----------------------------
tsv_to_vcf_body() {
  local in_tsv="$1" batch="$2" out_body="$3"

  gzip -cd "$in_tsv" \
  | awk -F'\t' -v OFS='\t' \
        -v BATCH="$batch" \
        -v PRIMARY="$primary_set_regex" '
    '"$map_filter"'
    NR==1 {
        for(i=1;i<=NF;i++){ k=$i; sub(/^#/,"",k); H[k]=i }
        # Required: locus, ref, alt, qual, filters, sample_id, fam_id, father_id, mother_id, child_GT, father_GT, mother_GT, num_alleles
        next
    }
    {
        split($1,a,":"); chrom=a[1]; sub(/^chr/,"",chrom); pos=a[2];

        # PRIMARY filter: keep only 1..22,X,Y
        c2=chrom;
        if (c2 !~ PRIMARY) next;
        filt = map_filter($(H["filters"]));
        info =  "SAMPLE=" $(H["sample_id"]) \
                    ";FAM=" $(H["fam_id"]) \
                    ";FATHER=" $(H["father_id"]) \
                    ";MOTHER=" $(H["mother_id"]) \
                    ";CHILD_GT=" $(H["child_GT"]) \
                    ";FATHER_GT=" $(H["father_GT"]) \
                    ";MOTHER_GT=" $(H["mother_GT"]) \
                    ";NUM_ALLELES=" $(H["num_alleles"]) \
                    ";BATCH=" BATCH
        print chrom, pos, ".", $(H["ref"]), $(H["alt"]), $(H["qual"]), filt, info
    }' > "$out_body"
}


# -----------------------------
# Function: assemble final VCF from body -> add contigs (22+XY) -> sort -> index
# -----------------------------
finalize_vcf() {
  local body="$1" out_vcf_gz="$2"

  {
    cat "$inh_header"
    cat "$body"
  } \
  | bcftools annotate -h "$data_dir/contigs.hdr" -Ou \
  | bcftools sort -Oz -o "$out_vcf_gz"

  tabix -f -p vcf "$out_vcf_gz"
}


# -----------------------------
# 1) Build primary-only, sorted VCFs per cohort
#    (SPARK shards become one SPARK VCF after concatenation)
# -----------------------------
tmp_bodies=()
for x in "${inherited_var[@]}"; do
  in_tsv="$wgs_var_dir/${x}.rare.inherited.tsv.gz"
  batch="$x"
  [[ $batch == "SSC_WGS9209s" ]] && batch="SSC_WGS_9209"

  body="$data_dir/${batch}.body.primary.vcf"
  tmp_bodies+=( "$body" )

  echo "[Convert] $x -> VCF body (primary only)"
  tsv_to_vcf_body "$in_tsv" "$batch" "$body"
done

# Combine the SPARK bodies into one body, then finalize once
if ls "$data_dir"/SPARK_WGS_*.body.primary.vcf >/dev/null 2>&1; then
  cat "$data_dir"/SPARK_WGS_*.body.primary.vcf > "$data_dir/SPARK_WGS.body.primary.vcf"
  rm -f "$data_dir"/SPARK_WGS_*.body.primary.vcf
fi

# Create per-cohort VCFs (sorted, indexed)
[[ -f "$data_dir/KOR-STAR.body.primary.vcf"     ]] && finalize_vcf "$data_dir/KOR-STAR.body.primary.vcf"     "$kor_inherited_vcf".gz
[[ -f "$data_dir/SPARK_WGS.body.primary.vcf"    ]] && finalize_vcf "$data_dir/SPARK_WGS.body.primary.vcf"    "$spark_inherited_vcf".gz
[[ -f "$data_dir/SSC_WGS_9209.body.primary.vcf" ]] && finalize_vcf "$data_dir/SSC_WGS_9209.body.primary.vcf" "$ssc_inherited_vcf".gz

# -----------------------------
# 2) Concatenate cohorts (already primary-only) -> all-inherited
# -----------------------------
inputs=()
[[ -f "$kor_inherited_vcf".gz  ]] && inputs+=( "$kor_inherited_vcf".gz )
[[ -f "$spark_inherited_vcf".gz ]] && inputs+=( "$spark_inherited_vcf".gz )
[[ -f "$ssc_inherited_vcf".gz  ]] && inputs+=( "$ssc_inherited_vcf".gz )

if ((${#inputs[@]})); then
  echo "[Concat] KOR + SPARK + SSC -> All-inherited (primary only)"
  bcftools concat -a -Oz -o "$inh_var_vcf".gz "${inputs[@]}"
  # Final tidy sort (ensures perfect contig/POS order across cohorts)
  bcftools sort -Oz -o "${inh_var_vcf%.vcf}.sorted.vcf.gz" "$inh_var_vcf".gz
  mv "${inh_var_vcf%.vcf}.sorted.vcf.gz" "$inh_var_vcf".gz
  tabix -f -p vcf "$inh_var_vcf".gz
fi

rm -f "$data_dir/"*.body.primary.vcf "$inh_header"

# Process Parents variants
echo "Processing Parents"
parents_var_vcf=$data_dir/UCSF_Parents_WGS_hg38.vcf
{
  printf '##fileformat=VCFv4.2\n'
  printf '##INFO=<ID=SAMPLE,Number=1,Type=String,Description="Sample_ID">\n'
  printf '#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\n'
} > "$data_dir/parents_header.vcf"

gzip -cd "$parents_var" \
| awk -F'\t' -v OFS='\t' -v PRIMARY="$primary_set_regex" '
NR==1{
  # map header -> index (robust to ordering)
  for(i=1;i<=NF;i++){k=$i; gsub(/^[#]+/,"",k); H[k]=i}
  next
}
{
  chrom=$(H["chr"]); sub(/^chr/,"",chrom)
  if (chrom !~ PRIMARY) next

  pos=$(H["pos"])
  ref=$(H["ref"])
  alt=$(H["alt"])
  sample=$(H["SampleID"])

  info="SAMPLE="sample;
  print chrom, pos, ".", ref, alt, ".", "PASS", info
}' > "$data_dir/parents_body.vcf"

cat "$data_dir/parents_header.vcf" "$data_dir/parents_body.vcf" \
| bcftools annotate -h "$data_dir/contigs.hdr" -Ou \
| bcftools sort -Oz -o "$parents_var_vcf".gz

tabix -f -p vcf "$parents_var_vcf".gz

rm -f "$data_dir/parents_header.vcf" "$data_dir/parents_body.vcf" "$data_dir/contigs.hdr"