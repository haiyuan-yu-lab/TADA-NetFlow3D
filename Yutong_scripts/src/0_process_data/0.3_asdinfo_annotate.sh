#!/bin/bash
#SBATCH --job-name=annotate_asdinfo
#SBATCH -o logs/annotate_asdinfo.%j.out
#SBATCH -e logs/annotate_asdinfo.%j.err


# Final step of data preprocessing:
# 1) Exclude siblings with ASD risk
# 2) Exclude multiplex individuals for DNV analysis
# 3) Exclude LOF & damaging missense variants:
#    a) LOF: - frameshift, splicing, stopgain, stoploss, startloss
#    b) Damaging missense: - missense consequence with REVEL >= 0.5
# 4) Annotate necessary metadata (sex, etc.)

set -euo pipefail
# -----------------------------
# Parameters
# -----------------------------
TYPE="dnm"  # dnm | inh
VEP_ANNO_DIR="/fs/cbsuhy01/storage/yz2676/projects/1-ASD/NetFlow_ASD/Yutong_scripts/data/annotated"
DNV_VEP_REL="dnm/UCSF_20211001_dnms_WGS_hg38.with_vep.vcf.gz"
INH_VEP_REL="inh/All_inherited_WGS_hg38.with_vep.vcf.gz"

PHENO="/fs/cbsuhy01/storage/NetFlow_ASD/Week8_429/SSC+SPARK+KOR-STAR_parsed.txt"
ASD_FAMILY="/fs/cbsuhy01/storage/NetFlow_ASD/Week8_429/roles_index.csv"
CONCERN_INDIV="/fs/cbsuhy01/storage/NetFlow_ASD/Week8_429/individuals_registration.csv"

FAMILY_TAG="FAM"
SAMPLE_TAG="SAMPLE"

# LOF & missense thresholds
LOF_REGEX="(frameshift_variant|splice_acceptor_variant|splice_donor_variant|stop_gained|stop_lost|start_lost)"
MISSENSE_TERM="missense_variant"
REVEL_THRESH="0.5"

# Work/output dir
ASD_DIR="/fs/cbsuhy01/storage/yz2676/projects/1-ASD/NetFlow_ASD/Yutong_scripts/data/asdinfo"

# Exclusion lists
FAM_LIST="$ASD_DIR/multiplex_families.txt"
SAMP_LIST="$ASD_DIR/samples_exclude.txt"
SEX_MAP="$ASD_DIR/sex.map"

THREADS="${THREADS:-8}"

# Toggles
EXCLUDE_MULTIPLEX="yes"     # yes|no
EXCLUDE_CONCERNING="yes"    # yes|no
EXCLUDE_CODING="yes"        # yes|no  (drop LOF + missense with REVEL>thres)

# -----------------------------
# Parse args
# -----------------------------
while [[ $# -gt 0 ]]; do
  case "$1" in
    --type)               TYPE="$2"; shift 2;;
    --asd-dir)            ASD_DIR="$2"; shift 2;;
    --vep-anno-dir)       VEP_ANNO_DIR="$2"; shift 2;;
    --dnv-vep)            DNV_VEP_REL="$2"; shift 2;;
    --inh-vep)            INH_VEP_REL="$2"; shift 2;;
    --pheno)              PHENO="$2"; shift 2;;
    --asd-family)         ASD_FAMILY="$2"; shift 2;;
    --concern-ind)        CONCERN_INDIV="$2"; shift 2;;
    --family-tag)         FAMILY_TAG="$2"; shift 2;;
    --sample-tag)         SAMPLE_TAG="$2"; shift 2;;
    --fam-list)           FAM_LIST="$2"; shift 2;;
    --samp-list)          SAMP_LIST="$2"; shift 2;;
    --sex-map)            SEX_MAP="$2"; shift 2;;
    --threads)            THREADS="$2"; shift 2;;
    --exclude-multiplex)  EXCLUDE_MULTIPLEX="$2"; shift 2;;
    --exclude-concerning) EXCLUDE_CONCERNING="$2"; shift 2;;
    --exclude-coding)     EXCLUDE_CODING="$2"; shift 2;;
    --revel)              REVEL_THRESH="$2"; shift 2;;
    --help|-h)            SHOW_HELP=1; shift;;
    --)                   shift; break;;
    *) echo "Unknown arg: $1" >&2; exit 1;;
  esac
done

mkdir -p "$ASD_DIR"
# -----------------------------
# Helpers
# -----------------------------
log(){ echo "[$(date '+%F %T')] $*" >&2; }
require_file(){ [[ -s "$1" ]] || { log "ERROR: Missing or empty file: $1"; exit 1; }; }
optional(){ [[ -s "$1" ]]; }

add_or_keep_sex_header(){
  awk '
    BEGIN { add=1 }
    /^##INFO=/ && $0 ~ /ID=SEX,/ { add=0 }
    /^#CHROM/ {
      if (add) print "##INFO=<ID=SEX,Number=1,Type=String,Description=\"Sample sex from mapping\">"
      print; next
    }
    { print }'
}

annotate_sex_stream(){
  local key="$1" mapfile="$2"
  awk -v FS="\t" -v OFS="\t" -v MAP="$mapfile" -v KEY="$key" '
    BEGIN{
      if (MAP=="") { print "ERROR: MAP is empty" > "/dev/stderr"; exit 1 }
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
      info=$8
      sid=get(info, KEY)
      if (sid==""){ noKey++; print; next }
      if (sid in sex){ $8=info";SEX="sex[sid]; changed++ } else { missing++ }
      print
    }
    END{
      print "INFO: rows edited:", changed,
            " | rows with " KEY " not in map:", missing,
            " | rows missing " KEY ":", noKey > "/dev/stderr"
    }'
}

build_multiplex_list(){
  # writes $FAM_LIST
  if [[ -s "$ASD_FAMILY" ]]; then
    csvgrep -c family_type -m Multiplex "$ASD_FAMILY" \
    | csvcut -c family_sf_id \
    | tail -n +2 | sort -u > "$FAM_LIST" || true
    log "Multiplex families -> $FAM_LIST ($(wc -l < "$FAM_LIST" 2>/dev/null || echo 0))"
  else
    log "WARN: ASD_FAMILY not found; cannot build multiplex list"
    : > "$FAM_LIST"
  fi
}

build_concerning_list(){
  # writes $SAMP_LIST (concern siblings + young siblings from Unknown families)
  local concern_sib="$ASD_DIR/concern_siblings.txt"
  local unknown_fam="$ASD_DIR/unknown_families.txt"
  local control_sib="$ASD_DIR/control_siblings.txt"
  local young_sib="$ASD_DIR/young_sib.txt"

  # concern_sib: control_sibling_flag==1
  if [[ -s "$CONCERN_INDIV" ]]; then
    csvgrep -c control_sibling_flag -m '1' "$CONCERN_INDIV" \
    | csvcut -c subject_sp_id \
    | tail -n +2 | sort -u > "$concern_sib" || true
  else
    : > "$concern_sib"
  fi

  # unknown_fam: family_type==Unknown (minus SF0153305)
  if [[ -s "$ASD_FAMILY" ]]; then
    csvgrep -c family_type -m Unknown "$ASD_FAMILY" \
    | csvcut -c family_sf_id \
    | awk '$0!="SF0153305"' | sort -u > "$unknown_fam" || true
  else
    : > "$unknown_fam"
  fi

  # control_sib: Father!="0" & Mother!="0" & Pheno=="control"  (columns 3,4,6)
  if [[ -s "$PHENO" ]]; then
    csvgrep -t -c 3 -i -r '^0$' "$PHENO" \
    | csvgrep -c 4 -i -r '^0$' \
    | csvgrep -c 6 -r '^control$' \
    > "$control_sib" || true
  else
    : > "$control_sib"
  fi

  # young_sib: join control_sib with unknown_fam on Fam, then take SampleID
  if [[ -s "$control_sib" && -s "$unknown_fam" ]]; then
    csvjoin -c Fam,family_sf_id "$control_sib" "$unknown_fam" \
    | csvcut -c SampleID \
    | tail -n +2 | awk 'NF' | sort -u > "$young_sib" || true
  else
    : > "$young_sib"
  fi

  cat "$young_sib" "$concern_sib" 2>/dev/null | sort -u > "$SAMP_LIST" || :
  log "Concerning samples -> $SAMP_LIST ($(wc -l < "$SAMP_LIST" 2>/dev/null || echo 0))"
}

# -----------------------------
# Inputs per TYPE
# -----------------------------
case "$TYPE" in
  dnm) IN_VCF="$VEP_ANNO_DIR/$DNV_VEP_REL"; OUT_PREFIX="$ASD_DIR/UCSF_20211001_dnms_WGS_hg38";;
  inh) IN_VCF="$VEP_ANNO_DIR/$INH_VEP_REL"; OUT_PREFIX="$ASD_DIR/All_inherited_WGS_hg38";;
  *)   echo "ERROR: --type must be dnm or inh" >&2; exit 1;;
esac

log "TYPE=$TYPE | IN_VCF=$IN_VCF | OUT_PREFIX=$OUT_PREFIX"
log "TOGGLES: exclude_multiplex=$EXCLUDE_MULTIPLEX, exclude_concerning=$EXCLUDE_CONCERNING, exclude_coding=$EXCLUDE_CODING"
log "THREADS=$THREADS"
log "OUTPUT dir: $ASD_DIR"

require_file "$IN_VCF"

# Build lists only if that exclusion is ON and no custom list provided (or default path used)
if [[ "$EXCLUDE_MULTIPLEX" == "yes" ]]; then
  if [[ ! -s "$FAM_LIST" || "$FAM_LIST" == "$ASD_DIR/multiplex_families.txt" ]]; then
    build_multiplex_list
  else
    log "Using provided multiplex family list: $FAM_LIST"
  fi
fi
if [[ "$EXCLUDE_CONCERNING" == "yes" ]]; then
  if [[ ! -s "$SAMP_LIST" || "$SAMP_LIST" == "$ASD_DIR/samples_exclude.txt" ]]; then
    build_concerning_list
  else
    log "Using provided concerning sample list: $SAMP_LIST"
  fi
fi

# -----------------------------
# Step 0: Convert to BCF + index
# -----------------------------
TMP_BCF="$OUT_PREFIX.tmp.bcf"
log "Converting to BCF: $TMP_BCF"
bcftools view --threads "$THREADS" -Ob -o "$TMP_BCF" "$IN_VCF"
bcftools index -f "$TMP_BCF"

# -----------------------------
# Step 1: Build position cache (CHROM,POS,POS,FAM,SAMPLE)
# -----------------------------
CACHE="$OUT_PREFIX.fam_sample.cache.bed.gz"
log "Building fam/sample cache: $CACHE"
bcftools query -f "%CHROM\t%POS\t%POS\t%INFO/${FAMILY_TAG}\t%INFO/${SAMPLE_TAG}\n" "$TMP_BCF" | bgzip > "$CACHE"
tabix -s1 -b2 -e3 "$CACHE"

# -----------------------------
# Step 2: Exclude multiplex families and concerning siblings
# -----------------------------
FILTERED_BCF="$OUT_PREFIX.filtered.bcf"
if [[ "$EXCLUDE_MULTIPLEX" == "yes" || "$EXCLUDE_CONCERNING" == "yes" ]]; then
  log "Applying selected exclusions -> $FILTERED_BCF"
  bcftools view -T ^<(
    {
      if [[ "$EXCLUDE_MULTIPLEX" == "yes" && -s "$FAM_LIST" ]]; then
        zcat "$CACHE" | awk -F'\t' 'NR==FNR{fam[$1]=1; next} { if($4 in fam) print $1"\t"$2"\t"$3 }' "$FAM_LIST" -
      fi
      if [[ "$EXCLUDE_CONCERNING" == "yes" && -s "$SAMP_LIST" ]]; then
        zcat "$CACHE" | awk -F'\t' 'NR==FNR{ok[$1]=1; next} { if($5 in ok) print $1"\t"$2"\t"$3 }' "$SAMP_LIST" -
      fi
    } | sort -u
  ) -Ob -o "$FILTERED_BCF" "$TMP_BCF"
else
  log "No FAM/SAMPLE exclusion selected; pass-through"
  bcftools view -Ob -o "$FILTERED_BCF" "$TMP_BCF"
fi
bcftools index -f "$FILTERED_BCF"

# =============================
# Step 3: Optional coding exclusion (LOF + high-REVEL missense)
# =============================
NONCODING_BCF="$OUT_PREFIX.filtered.noncoding.bcf"
if [[ "$EXCLUDE_CODING" == "yes" ]]; then
  log "Excluding LOF & missense(REVEL>$REVEL_THRESH) -> $NONCODING_BCF"
  bcftools +split-vep \
    -c REVEL:Float \
    -i "(Consequence!~\"$LOF_REGEX\") && (Consequence!~\"$MISSENSE_TERM\" || REVEL<=$REVEL_THRESH)" \
    -Ou "$FILTERED_BCF" \
  | bcftools annotate -x 'INFO/Consequence,INFO/REVEL' \
    -Ob -o "$NONCODING_BCF"
else
  log "Skipping coding exclusion; pass-through -> $NONCODING_BCF"
  bcftools view -Ob -o "$NONCODING_BCF" "$FILTERED_BCF"
fi
bcftools index -f "$NONCODING_BCF"

# =============================
# Step 4: Counts
# =============================
N_IN=$(bcftools view -H "$TMP_BCF" | wc -l || echo 0)
N_AFTER_FS=$(bcftools view -H "$FILTERED_BCF" | wc -l || echo 0)
N_AFTER_CODE=$(bcftools view -H "$NONCODING_BCF" | wc -l || echo 0)
log "Counts | Input=$N_IN | After FAM/SAMPLE=$N_AFTER_FS | After coding exclusion=$N_AFTER_CODE"

# =============================
# Step 5: SEX annotation
# =============================
FINAL_BCF="$OUT_PREFIX.filtered.noncoding.with_asdinfo.bcf"
log "Annotating SEX using $SEX_MAP -> $FINAL_BCF"
require_file "$SEX_MAP"

{
  bcftools view -h "$NONCODING_BCF" | add_or_keep_sex_header
  bcftools view -H "$NONCODING_BCF" | annotate_sex_stream "$SAMPLE_TAG" "$SEX_MAP"
} | bcftools view --threads "$THREADS" -Ob -o "$FINAL_BCF" -
bcftools index -f "$FINAL_BCF"

log "DONE. Final: $FINAL_BCF"

# =============================
# Cleanup
# =============================
log "Cleaning up temporary files"
rm -f "$TMP_BCF" "$TMP_BCF".csi
rm -f "$FILTERED_BCF" "$FILTERED_BCF".csi
rm -f "$NONCODING_BCF" "$NONCODING_BCF".csi
rm -f "$CACHE" "$CACHE".tbi
log "All done."
