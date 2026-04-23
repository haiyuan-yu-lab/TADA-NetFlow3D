#!/bin/bash
#SBATCH --job-name=untrans_add_batch
#SBATCH -o logs/untrans_add_batch.%j.out
#SBATCH -e logs/untrans_add_batch.%j.err

#SBATCH -N 1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=4
#SBATCH --mem=32G
#SBATCH --time=12:00:00
#SBATCH --nodelist=cbsuhy02

set -euo pipefail

log(){ echo "[$(date '+%F %T')] $*" >&2; }

# ---- DEFAULT CONFIG ----
ASD_DIR="/fs/cbsuhy01/storage/yz2676/projects/1-ASD/NetFlow_ASD/Yutong_scripts/data/asdinfo"
OUT_PREFIX=""
IN_BCF=""
OUT_BCF=""
THREADS="${SLURM_CPUS_PER_TASK:-4}"
# ------------------------

while [[ $# -gt 0 ]]; do
  case "$1" in
    --asd-dir)    ASD_DIR="$2"; shift 2 ;;
    --out-prefix) OUT_PREFIX="$2"; shift 2 ;;
    --in-bcf)     IN_BCF="$2"; shift 2 ;;
    --out-bcf)    OUT_BCF="$2"; shift 2 ;;
    --threads)    THREADS="$2"; shift 2 ;;
    *) echo "Unknown arg: $1" >&2; exit 1 ;;
  esac
done

if [[ -z "$OUT_PREFIX" ]]; then
  OUT_PREFIX="$ASD_DIR/UCSF_parents_WGS_hg38.includesimplex"
fi

if [[ -z "$IN_BCF" ]]; then
  IN_BCF="$OUT_PREFIX.untransmitted.dedup.bcf"
fi

if [[ -z "$OUT_BCF" ]]; then
  OUT_BCF="$OUT_PREFIX.untransmitted.dedup.with_batch.bcf"
fi

[[ -s "$IN_BCF" ]] || { log "Missing input BCF: $IN_BCF"; exit 1; }

out_dir="$(dirname "$OUT_BCF")"
mkdir -p "$out_dir"

add_batch_header(){
  awk '
    BEGIN { add_batch=1 }
    /^##INFO=/ && $0 ~ /ID=BATCH,/ { add_batch=0 }
    /^#CHROM/ {
      if (add_batch)
        print "##INFO=<ID=BATCH,Number=1,Type=String,Description=\"Batch inferred from SAMPLE prefix\">"
      print
      next
    }
    { print }'
}

log "Annotating INFO/BATCH from INFO/SAMPLE"
log "Input:  $IN_BCF"
log "Output: $OUT_BCF"

{
  bcftools view -h "$IN_BCF" | add_batch_header
  bcftools view -H "$IN_BCF" \
  | awk -v FS="\t" -v OFS="\t" '
      function get(info, key, a){
        return match(info, "(^|;)" key "=([^;]*)", a) ? a[2] : ""
      }
      function infer_batch(sample){
        if (sample ~ /^SP/) return "SPARK"
        if (sample ~ /^SSC/) return "SSC_WGS_9209"
        if (sample ~ /^SS/) return "SSC_WGS_9209"
        if (sample ~ /^KOR-/) return "KOR_STAR"
        if (sample ~ /^STAR-/) return "KOR_STAR"
        return ""
      }
      {
        sample = get($8, "SAMPLE")
        if (sample == ""){
          print "ERROR: missing SAMPLE in INFO for record " $1 ":" $2 ":" $4 ":" $5 > "/dev/stderr"
          exit 1
        }

        batch = infer_batch(sample)
        if (batch == ""){
          print "ERROR: unable to infer BATCH from SAMPLE=" sample " at " $1 ":" $2 ":" $4 ":" $5 > "/dev/stderr"
          exit 1
        }

        if (get($8, "BATCH") != ""){
          print "ERROR: input record already contains BATCH for SAMPLE=" sample " at " $1 ":" $2 ":" $4 ":" $5 > "/dev/stderr"
          exit 1
        }

        $8 = $8 ";BATCH=" batch
        print
      }'
} | bcftools view --threads "$THREADS" -Ob -o "$OUT_BCF" -

bcftools index -f "$OUT_BCF"
log "DONE. Annotated untransmitted deduplicated BCF: $OUT_BCF"
