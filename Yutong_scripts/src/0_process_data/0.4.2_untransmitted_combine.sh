#!/bin/bash
#SBATCH --job-name=untrans_finish
#SBATCH -o logs/untrans_finish.%j.out
#SBATCH -e logs/untrans_finish.%j.err

#SBATCH -N 1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=8
#SBATCH --mem=128G
#SBATCH --time=24:00:00
#SBATCH --nodelist=cbsuhy02

set -euo pipefail
log(){ echo "[$(date '+%F %T')] $*" >&2; }

# ---- CONFIG (edit these) ----
ASD_DIR="/fs/cbsuhy01/storage/yz2676/projects/1-ASD/NetFlow_ASD/Yutong_scripts/data/asdinfo"
OUT_PREFIX="$ASD_DIR/UCSF_parents_WGS_hg38.includesimplex.codingonly"

FINAL_PARENT_BCF="$OUT_PREFIX.filtered.vep.with_asdinfo.bcf"
OUT_BODY_DIR="$ASD_DIR/untrans_bodies"
NCHUNKS=25

UNTRANSMITTED_BCF="$OUT_PREFIX.untransmitted.bcf"
ALL_BODY="$OUT_PREFIX.untransmitted.body.all.vcf"
THREADS="${SLURM_CPUS_PER_TASK:-8}"
# -----------------------------

[[ -s "$FINAL_PARENT_BCF" ]] || { log "Missing FINAL_PARENT_BCF: $FINAL_PARENT_BCF"; exit 1; }

# Make sure all chunk outputs exist
missing=0
for i in $(seq -w 0 $((NCHUNKS-1))); do
  ii=$((10#$i))
  f=$(printf "%s/%s.untransmitted.body.%04d.vcf" "$OUT_BODY_DIR" "$(basename "$OUT_PREFIX")" "$ii")
  if [[ ! -s "$f" ]]; then
    log "Missing body: $f"
    missing=$((missing+1))
  fi
done
(( missing==0 )) || { log "ERROR: $missing chunks missing. Not finishing."; exit 1; }

log "Concatenating bodies -> $ALL_BODY"
: > "$ALL_BODY"
for i in $(seq -w 0 $((NCHUNKS-1))); do
  ii=$((10#$i))
  f=$(printf "%s/%s.untransmitted.body.%04d.vcf" "$OUT_BODY_DIR" "$(basename "$OUT_PREFIX")" "$ii")
  cat "$f" >> "$ALL_BODY"
done

# Header helper
add_parent_untrans_headers(){
  awk '
    BEGIN { add_off=1; add_status=1 }
    /^##INFO=/ && $0 ~ /ID=OFFSPRING,/   { add_off=0 }
    /^##INFO=/ && $0 ~ /ID=STATUS,/ { add_status=0 }
    /^#CHROM/ {
      if (add_off)
        print "##INFO=<ID=OFFSPRING,Number=1,Type=String,Description=\"Child SampleID this parent variant is compared against\">"
      if (add_status)
        print "##INFO=<ID=STATUS,Number=1,Type=String,Description=\"Child phenotype from metadata\">"
      print; next
    }
    { print }'
}

log "Building final BCF -> $UNTRANSMITTED_BCF"
tmpdir="$ASD_DIR/tmp"
mkdir -p "$tmpdir"
{
  bcftools view -h "$FINAL_PARENT_BCF" | add_parent_untrans_headers
  cat "$ALL_BODY"
} | bcftools view --threads "$THREADS" -Ob - \
  | bcftools sort -T "$tmpdir" -m 8G -Ob -o "$UNTRANSMITTED_BCF" -

bcftools index -f "$UNTRANSMITTED_BCF"
log "DONE: $UNTRANSMITTED_BCF"

# remove intermediate files
rm "$ALL_BODY"

# Deduplicate
log "Creating deduplicated per-variant-per-offspring untransmitted VCF -> $OUT_PREFIX.untransmitted.dedup.bcf"
DEDUP_BCF="$OUT_PREFIX.untransmitted.dedup.bcf"

bcftools view "$UNTRANSMITTED_BCF" \
| awk -v FS="\t" -v OFS="\t" '
    function get(info,key,a){ return match(info, "(^|;)"key"=([^;]*)", a) ? a[2] : "" }
    /^#/ { print; next }   # header lines: print and skip parsing
    {
    info = $8
    off=get(info, "OFFSPRING")
    var=$1"\t"$2"\t"$4"\t"$5   # CHROM POS REF ALT
    key=var"\t"off
    if (!(key in SEEN)){
        SEEN[key]=1
        print
    }
    }' \
| bcftools view --threads "$THREADS" -Ob -o "$DEDUP_BCF" -

bcftools index -f "$DEDUP_BCF"
log "DONE. Deduplicated untransmitted VCF: $DEDUP_BCF"
