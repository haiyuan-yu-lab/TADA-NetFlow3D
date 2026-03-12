#!/bin/bash
#SBATCH --job-name=annotate_coding
#SBATCH -o logs/annotate_coding.%j.out
#SBATCH -e logs/annotate_coding.%j.err
#SBATCH --mail-type=ALL
#SBATCH --mail-user=yz2676@cornell.edu

#SBATCH -N 1
#SBATCH --ntasks-per-node=1
#SBATCH --cpus-per-task=8
#SBATCH --mem=24G
#SBATCH --time=64:00:00
#SBATCH --nodelist=cbsuhy02

source src/0_process_data/vep.env
source src/0_process_data/vep_common.sh
set -euo pipefail
log(){ echo "[$(date '+%F %T')] $*" >&2; }

# -----------------------------
# Parameters
# -----------------------------
TYPE="dnm"  # dnm | inh | parent

dnm_with_vep="$DNM_OUT_DIR/UCSF_20211001_dnms_WGS_hg38.with_vep.vcf.gz"
inherited_with_vep="$INH_OUT_DIR/All_inherited_WGS_hg38.with_vep.vcf.gz"
parent_with_vep="$PARENTS_OUT_DIR/UCSF_parents_WGS_hg38.with_vep.vcf.gz"

# -----------------------------
# Inputs per TYPE
# -----------------------------
while [[ $# -gt 0 ]]; do
  case "$1" in
    --type)               TYPE="$2"; shift 2;;
    *) echo "Unknown argument: $1" >&2; exit 1 ;;
    esac
done

case "$TYPE" in
  dnm)
    in_vcf="$dnm_with_vep"
    out_prefix="$DNM_OUT_DIR/UCSF_20211001_dnms_WGS_hg38.coding"
    ;;
  inh)
    in_vcf="$inherited_with_vep"
    out_prefix="$INH_OUT_DIR/All_inherited_WGS_hg38.coding"
    ;;
  parent)
    in_vcf="$parent_with_vep"
    out_prefix="$PARENTS_OUT_DIR/UCSF_parents_WGS_hg38.coding"
    ;;
  *) echo "ERROR: --type must be dnm or inh or parent" >&2; exit 1 ;;
esac

# -----------------------------
# Functions
# -----------------------------
workdir="$(mktemp -d -p "$TMP_DIR" annotate_coding.${TYPE}.XXXXXX)"
log "Working directory: $workdir"
# trap 'rm -rf "$workdir"' EXIT

cat > $workdir/addchr.txt <<'EOF'
1	chr1
2	chr2
3	chr3
4	chr4
5	chr5
6	chr6
7	chr7
8	chr8
9	chr9
10	chr10
11	chr11
12	chr12
13	chr13
14	chr14
15	chr15
16	chr16
17	chr17
18	chr18
19	chr19
20	chr20
21	chr21
22	chr22
X	chrX
Y	chrY
EOF

cat > $workdir/rmchr.txt <<'EOF'
chr1	1
chr2	2
chr3	3
chr4	4
chr5	5
chr6	6
chr7	7
chr8	8
chr9	9
chr10	10
chr11	11
chr12	12
chr13	13
chr14	14
chr15	15
chr16	16
chr17	17
chr18	18
chr19	19
chr20	20
chr21	21
chr22	22
chrX	X
chrY	Y
EOF

# -----------------------------
# Step 1) Extract hg38 coding (BIOTYPE=protein_coding) from VEP-annotated VCF
#         (and ensure ID = CHROM:POS:REF:ALT for join key)
# -----------------------------
coding_hg38="${out_prefix}.extracted.id.vcf.gz"
log "Extracting coding and setting ID: $coding_hg38"

tmp_extract="$workdir/coding.extracted.vcf.gz"
extract_coding "$in_vcf" "$tmp_extract"

bcftools annotate --set-id '%CHROM:%POS:%REF:%ALT' \
  -Oz -o "$coding_hg38" "$tmp_extract"
tabix -f -p vcf "$coding_hg38"

# -----------------------------
# Step 2) Build unique-sites VCF (hg38), then liftover to hg19
# -----------------------------
uniq_hg38="$workdir/coding.unique.hg38.vcf.gz"
uniq_hg38_addchr="$workdir/coding.unique.hg38.addchr.vcf.gz"
lifted_addchr_hg19="$workdir/coding.unique.hg38.addchr.hg19.vcf.gz"
lifted_hg19="$workdir/coding.unique.hg19.vcf.gz"
rejected_hg19="$workdir/coding.unique.unmapped.hg19.vcf.gz"

log "Creating unique-sites VCF (hg38) for liftover"
make_unique_sites_vcf "$coding_hg38" "$uniq_hg38"

bcftools annotate --rename-chrs "$workdir/addchr.txt" \
  -Oz -o "$uniq_hg38_addchr" "$uniq_hg38"
tabix -f -p vcf "$uniq_hg38_addchr"

log "Liftover unique-sites hg38 -> hg19 (Picard)"
picard -Xmx16g LiftoverVcf \
  I="$uniq_hg38_addchr" \
  O="$lifted_addchr_hg19" \
  CHAIN="$LIFTOVER_CHAIN" \
  REJECT="$rejected_hg19" \
  R="$HG19_FA" \
  RECOVER_SWAPPED_REF_ALT=true \
  WRITE_ORIGINAL_POSITION=true \
  WRITE_ORIGINAL_ALLELES=true
tabix -f -p vcf "$lifted_addchr_hg19"

# remove chr prefix back to match your GRCh37 VEP/cache convention
bcftools annotate --rename-chrs "$workdir/rmchr.txt" \
    -Oz -o "$lifted_hg19" "$lifted_addchr_hg19"
tabix -f -p vcf "$lifted_hg19"

log "Unmapped during liftover to hg19: $(bcftools view -H "$rejected_hg19" | wc -l)"

# -----------------------------
# Step 3) Run VEP (GRCh37) with MPC plugin on hg19 sites
# -----------------------------
hg19_vep="$workdir/hg19.vep.vcf.gz"
log "Running GRCh37 VEP+MPC on liftover hg19 unique sites"
run_vep_with_mpc_grch37 "$lifted_hg19" "$hg19_vep" 4 20000
tabix -f -p vcf "$hg19_vep"

# -----------------------------
# Step 4) Extract MPC per transcript, take max per variant ID
#         Output: mpc_byid = (ID \t MPCmax)
# -----------------------------
mpc_byid="$workdir/mpc_by_id.tsv.gz"
log "Computing MPC max per variant ID from hg19 VEP output"

# split-vep to pull MPC into INFO/MPC (comma-separated across transcripts) then query
bcftools +split-vep "$hg19_vep" -s all -c "MPC" -Ou \
| bcftools query -f '%ID\t%INFO/MPC\n' \
| awk -v OFS="\t" '
  function max_mpc(s,a,i,v,m,n){
    n=split(s,a,","); m=""
    for(i=1;i<=n;i++){
      if(a[i]=="." || a[i]=="NA" || a[i]=="") continue
      v=a[i]+0
      if(m=="" || v>m) m=v
    }
    return (m=="" ? "." : m)
  }
  { print $1, max_mpc($2) }
' \
| bgzip -c > "$mpc_byid"

log "mpc_byid lines: $(zcat "$mpc_byid" | wc -l)"
# -----------------------------
# Step 5) Join MPCmax back to hg38 coding VCF by ID and annotate INFO/MPC
#         Build annotation table: CHROM POS ID MPC
# -----------------------------
coding_hg38_vep="$workdir/hg38.coding.vep_addedinfo.vcf.gz"
log "Adding VEP INFO fields back to hg38 coding VCF for easier join"
run_vep_grch38_added "$coding_hg38" "$coding_hg38_vep" 8 80000
tabix -f -p vcf "$coding_hg38_vep"

mpc_annot_hg38="$workdir/hg38.mpc_annot.tsv.gz"
log "Joining MPC max back to hg38 coding variants by ID"

# left: CHROM POS REF ALT ID (sorted by ID)
# right: ID MPC (sorted by ID)
LC_ALL=C join -a1 -e '.' -t $'\t' -1 5 -2 1 \
  <(bcftools query -f '%CHROM\t%POS\t%REF\t%ALT\t%ID\n' "$coding_hg38_vep" \
    | LC_ALL=C sort -t$'\t' -k5,5) \
  <(zcat "$mpc_byid" | tr -d '\r' | LC_ALL=C sort -t$'\t' -k1,1) \
| awk -v OFS="\t" '{print $2,$3,$4,$5,$6}' \
| LC_ALL=C sort -t$'\t' -k1,1V -k2,2n \
| bgzip -c > "$mpc_annot_hg38"

tabix -f -s 1 -b 2 -e 2 "$mpc_annot_hg38"
log "hg38 MPC annotation rows: $(zcat "$mpc_annot_hg38" | wc -l)"

out_vcf="${out_prefix}.with_mpcmax.vcf.gz"
log "Annotating hg38 coding VCF with INFO/MPC -> $out_vcf"

# Add INFO header for MPCmax and annotate
mpc_hdr="$workdir/mpcmax.header.txt"
cat > "$mpc_hdr" <<'EOF'
##INFO=<ID=MPC,Number=1,Type=Float,Description="Max MPC across transcripts (computed from GRCh37 VEP MPC list; joined back to GRCh38 by variant ID)">
EOF

bcftools annotate \
  -a "$mpc_annot_hg38" \
  -c CHROM,POS,REF,ALT,INFO/MPC \
  -h "$mpc_hdr" \
  -Oz -o "$out_vcf" \
  "$coding_hg38_vep"
tabix -f -p vcf "$out_vcf"

# -----------------------------
# Sanity checks
# -----------------------------
log "Coding variants (hg38 extracted): $(bcftools view -H "$coding_hg38_vep" | wc -l)"
log "With MPC annotated (MPC != .): $(bcftools view -H -i 'INFO/MPC!="."' "$out_vcf" | wc -l)"
log "Done: $out_vcf"

# clean
rm -rf "$workdir"