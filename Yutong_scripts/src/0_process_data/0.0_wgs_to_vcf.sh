#!/bin/bash
# Transform WGS data (DNV, Inherited) to VCF format

## Raw files
wgs_var_dir=/fs/cbsuhy01/storage/NetFlow_ASD/Week8_429/UCSF_raw

dnv_var=UCSF_20211001_dnms_WGS_hg38.tsv
inherited_var=("KOR-STAR" "SPARK_WGS_1" "SPARK_WGS_2" "SPARK_WGS_3" "SSC_WGS9209s")
inherited_tasks=()
for x in "${inherited_var[@]}"; do
  inherited_tasks+=("$wgs_var_dir/${x}.rare.inherited.tsv.gz")
done

## Reference
ref=/fs/cbsuhy01/storage/yz2676/data/DetailedProCapNet/wdir/reference/hg38.22XY.fasta

## Outputs
data_dir=data/raw_processed
mkdir -p $data_dir

dnv_var_vcf=$data_dir/UCSF_20211001_dnms_WGS_hg38.vcf
inh_var_vcf=$data_dir/All_inherited_WGS_hg38.vcf

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

### 5. Clean up
rm $data_dir/dnv_header.vcf $data_dir/dnv_body.vcf $dnv_var_vcf.unsorted.gz $data_dir/contigs.hdr $dnv_var_vcf

