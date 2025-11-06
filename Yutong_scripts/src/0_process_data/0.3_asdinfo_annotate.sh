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

# Usage:
# 1) For DNVs:
# sbatch src/0_process_data/0.3_asdinfo_annotate.sh --type dnm 
#    --exclude-multiplex yes --exclude-concerning yes --exclude-coding yes
# 2) For inherited variants:
# sbatch src/0_process_data/0.3_asdinfo_annotate.sh --type inh --exclude-multiplex no --keep-only-multiplex yes
# 3) DNV with coding variants kept:
# sbatch src/0_process_data/0.3_asdinfo_annotate.sh --type dnm --exclude-coding no --out-prefix UCSF_20211001_dnms_WGS_hg38.codingkept

set -eo pipefail
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
LOF_REGEX="frameshift_variant\|splice_acceptor_variant\|splice_donor_variant\|stop_gained\|stop_lost\|start_lost\|inframe_deletion\|inframe_insertion\|protein_altering_variant"
MISSENSE_TERM="missense_variant"
REVEL_THRESH=0.5

# Work/output dir
ASD_DIR="/fs/cbsuhy01/storage/yz2676/projects/1-ASD/NetFlow_ASD/Yutong_scripts/data/asdinfo"

THREADS="${THREADS:-8}"

# Toggles
EXCLUDE_MULTIPLEX="yes"     # yes|no
EXCLUDE_CONCERNING="yes"    # yes|no
EXCLUDE_CODING="yes"        # yes|no  (drop LOF + missense with REVEL>thres)
KEEP_ONLY_MULTIPLEX="no"    # yes|no

# -----------------------------
# Parse args
# -----------------------------
while [[ $# -gt 0 ]]; do
  case "$1" in
    --type)               TYPE="$2"; shift 2;;
    --asd-dir)            ASD_DIR="$2"; shift 2;;
    --out-prefix)         USER_PREFIX="$2"; shift 2;;
    --vep-anno-dir)       VEP_ANNO_DIR="$2"; shift 2;;
    --dnv-vep)            DNV_VEP_REL="$2"; shift 2;;
    --inh-vep)            INH_VEP_REL="$2"; shift 2;;
    --pheno)              PHENO="$2"; shift 2;;
    --asd-family)         ASD_FAMILY="$2"; shift 2;;
    --concern-ind)        CONCERN_INDIV="$2"; shift 2;;
    --family-tag)         FAMILY_TAG="$2"; shift 2;;
    --sample-tag)         SAMPLE_TAG="$2"; shift 2;;
    --threads)            THREADS="$2"; shift 2;;
    --exclude-multiplex)  EXCLUDE_MULTIPLEX="$2"; shift 2;;
    --exclude-concerning) EXCLUDE_CONCERNING="$2"; shift 2;;
    --exclude-coding)     EXCLUDE_CODING="$2"; shift 2;;
    --keep-only-multiplex) KEEP_ONLY_MULTIPLEX="$2"; shift 2;;
    --revel)              REVEL_THRESH="$2"; shift 2;;
    --help|-h)            SHOW_HELP=1; shift;;
    --)                   shift; break;;
    *) echo "Unknown arg: $1" >&2; exit 1;;
  esac
done

mkdir -p "$ASD_DIR"

# Exclusion lists
FAM_LIST="$ASD_DIR/multiplex_families.txt"
SAMP_LIST="$ASD_DIR/samples_exclude.txt"
SEX_MAP="$ASD_DIR/sex.map"
PHENO_MAP="$ASD_DIR/pheno.map"

# -----------------------------
# Helpers
# -----------------------------
log(){ echo "[$(date '+%F %T')] $*" >&2; }
require_file(){ [[ -s "$1" ]] || { log "ERROR: Missing or empty file: $1"; exit 1; }; }
optional(){ [[ -s "$1" ]]; }

add_or_keep_sex_pheno_headers(){
  awk '
    BEGIN { add_sex=1; add_status=1 }
    /^##INFO=/ && $0 ~ /ID=SEX,/   { add_sex=0 }
    /^##INFO=/ && $0 ~ /ID=STATUS,/ { add_status=0 }
    /^#CHROM/ {
      if (add_sex)
        print "##INFO=<ID=SEX,Number=1,Type=String,Description=\"Sample sex from mapping\">"
      if (add_status)
        print "##INFO=<ID=STATUS,Number=1,Type=String,Description=\"Phenotype from metadata\">"
      print; next
    }
    { print }'
}

annotate_sex_pheno_stream(){
  local key="$1" sexmap="$2" phenomap="$3"
  awk -v FS="\t" -v OFS="\t" \
      -v KEY="$key" -v S="$sexmap" -v P="$phenomap" '
    BEGIN{
      if (S=="") { print "ERROR: SEX MAP is empty" > "/dev/stderr"; exit 1 }
      while((getline < S) > 0){
        sub(/\r$/,"");
        split($0,f,"\t"); if(f[1] && f[2]) { sex[f[1]]=f[2]; nsex++ }
      }
      close(S)
      if (P!=""){
        while((getline < P) > 0){ sub(/\r$/,""); split($0,f,"\t"); if(f[1] && f[2]) { pheno[f[1]]=f[2]; nph++ } }
        close(P)
      }
      print "INFO: loaded", nsex, "sample->sex pairs" > "/dev/stderr"
      print "INFO: loaded", nph+0, "sample->pheno pairs" > "/dev/stderr"
      sex_added=0; pheno_added=0; sex_missing=0; pheno_missing=0; noKey=0
    }
    function get(info,key, a){
      if (match(info, "(^|;)" key "=([^;]*)", a)) return a[2]
      return ""
    }
    {
      info=$8
      sid=get(info, KEY)
      if (sid==""){ noKey++; print; next }
      if (sid in sex){ info=info";SEX="sex[sid]; sex_added++ } else { sex_missing++ }
      if (P!="" && (sid in pheno)){ info=info";STATUS="pheno[sid]; pheno_added++ } else { pheno_missing++ }
      $8=info
      print
    }
    END{
      print "INFO: rows missing " KEY ":", noKey+0 \
          , " | SEX added:", sex_added+0, "missing:", sex_missing+0 \
          , " | PHENO added:", pheno_added+0, "missing:", pheno_missing+0 > "/dev/stderr"
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

build_sex_map(){
  log "Building SEX map -> $SEX_MAP"
  if [[ -s "$PHENO" ]]; then
    awk -F'\t' 'NR==1{
    for(i=1;i<=NF;i++){
        gsub(/\r$/,"",$i)
        if($i=="SampleID") sid=i
        if($i=="sex") sx=i
    }
    next
    }
    {print $sid "\t" $sx}' "$PHENO" > "$SEX_MAP"
  else
    : > "$SEX_MAP"
  fi
}

build_pheno_map(){
  log "Building PHENO map -> $PHENO_MAP"
  if [[ -s "$PHENO" ]]; then
    awk -F'\t' 'NR==1{
    for(i=1;i<=NF;i++){
        gsub(/\r$/,"",$i)
        if($i=="SampleID") sid=i
        if($i=="Pheno") ph=i
    }
    next
    }
    {print $sid "\t" $ph}' "$PHENO" > "$PHENO_MAP"
  else
    : > "$PHENO_MAP"
  fi
}

# -----------------------------
# Inputs per TYPE
# -----------------------------
case "$TYPE" in
  dnm) IN_VCF="$VEP_ANNO_DIR/$DNV_VEP_REL";;
  inh) IN_VCF="$VEP_ANNO_DIR/$INH_VEP_REL";;
  *)   echo "ERROR: --type must be dnm or inh" >&2; exit 1;;
esac
require_file "$IN_VCF"

if [[ -n "$USER_PREFIX" ]]; then
    OUT_PREFIX="$ASD_DIR/$(basename "$USER_PREFIX")"
else
    case "$TYPE" in
    dnm) OUT_PREFIX="$ASD_DIR/UCSF_20211001_dnms_WGS_hg38";;
    inh) OUT_PREFIX="$ASD_DIR/All_inherited_WGS_hg38";;
    esac
fi

log "[INFO] TYPE=$TYPE | IN_VCF=$IN_VCF"
log "[INFO] OUT_PREFIX=$OUT_PREFIX"
log "[INFO] TOGGLES: exclude_multiplex=$EXCLUDE_MULTIPLEX, exclude_concerning=$EXCLUDE_CONCERNING, exclude_coding=$EXCLUDE_CODING, keep_only_multiplex=$KEEP_ONLY_MULTIPLEX"
log "[INFO] THREADS=$THREADS"

# Build lists only if that exclusion is ON and no custom list provided (or default path used)
if [[ "$EXCLUDE_MULTIPLEX" == "yes" && "$KEEP_ONLY_MULTIPLEX" == "yes" ]]; then
  log "ERROR: exclude_multiplex and keep_only_multiplex are mutually exclusive"; exit 1
fi

if [[ "$EXCLUDE_MULTIPLEX" == "yes" ]]; then
  if [[ ! -s "$FAM_LIST" || "$FAM_LIST" == "$ASD_DIR/multiplex_families.txt" ]]; then
    build_multiplex_list
  else
    log "Using provided multiplex family list: $FAM_LIST"
  fi
fi
if [[ "$KEEP_ONLY_MULTIPLEX" == "yes" ]]; then
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
# Step 1: Exclude/keep multiplex families and concerning siblings
# -----------------------------
FILTERED_BCF="$OUT_PREFIX.filtered.bcf"
log "Applying selected exclusions -> $FILTERED_BCF"
tmp_kept_vcf="$OUT_PREFIX.kept.tmp.vcf.gz"

{ bcftools view -h "$IN_VCF"
  bcftools view -H "$IN_VCF" \
  | awk -v FS="\t" -v OFS="\t" \
         -v fam_list="$FAM_LIST" -v samp_list="$SAMP_LIST" \
         -v use_fam="$EXCLUDE_MULTIPLEX" -v use_samp="$EXCLUDE_CONCERNING" \
         -v keep_only="$KEEP_ONLY_MULTIPLEX" \
         -v fam_tag="$FAMILY_TAG" -v samp_tag="$SAMPLE_TAG" '
  BEGIN{
    need_fam = (tolower(use_fam)=="yes" || tolower(keep_only)=="yes")
    if (need_fam && system("[ -s \""fam_list"\" ]")==0){
      while ((getline l < fam_list)>0){ sub(/\r$/,"",l); if(l!="") F[l]=1 }
      close(fam_list)
    }
    if (use_samp=="yes" && system("[ -s \""samp_list"\" ]")==0){
      while ((getline l < samp_list)>0){ sub(/\r$/,"",l); if(l!="") S[l]=1 }
      close(samp_list)
    }
  }
  function get(info,key,a){ return match(info, "(^|;)"key"=([^;]*)", a) ? a[2] : "" }
  # does any comma-separated token exist in assoc array A?
  function any_in(val,A){ if(val=="") return 0; { if(val in A) return 1 } return 0 }
  {
    info = $8
    fam  = get(info, fam_tag)
    samp = get(info, samp_tag)
    drop = 0
    if (tolower(keep_only)=="yes"){ if (!any_in(fam, F)) drop=1}
    if (!drop && tolower(use_fam)=="yes"){ if (any_in(fam, F)) drop=1}
    if (!drop && tolower(use_samp)=="yes"){ if (any_in(samp, S)) drop=1 }
    if (!drop) print
  }'
} | bgzip -c > "$tmp_kept_vcf"

tabix -f -p vcf "$tmp_kept_vcf"
bcftools view -Ob -o "$FILTERED_BCF" "$tmp_kept_vcf"
bcftools index -f "$FILTERED_BCF"

# =============================
# Step 2: Optional coding exclusion (LOF + high-REVEL missense)
# =============================  
NONCODING_BCF="$OUT_PREFIX.filtered.vep.bcf"
if [[ "$EXCLUDE_CODING" == "yes" ]]; then
  log "Excluding LOF & missense(REVEL>=$REVEL_THRESH) -> $NONCODING_BCF"
  bcftools +split-vep \
    -c Consequence,REVEL:Float \
    -e "(Consequence~\"$LOF_REGEX\") || (Consequence~\"$MISSENSE_TERM\" && REVEL>=$REVEL_THRESH)" \
    -Ou "$FILTERED_BCF" \
    | bcftools annotate -x 'INFO/Consequence,INFO/REVEL' \
      -Ob -o "$NONCODING_BCF"
else
  log "Skipping coding exclusion; pass-through -> $NONCODING_BCF"
  bcftools view -Ob -o "$NONCODING_BCF" "$FILTERED_BCF"
fi
bcftools index -f "$NONCODING_BCF"

# =============================
# Step 3: Counts
# =============================
N_IN=$(bcftools view -H "$IN_VCF" | wc -l || echo 0)
N_AFTER_FS=$(bcftools view -H "$FILTERED_BCF" | wc -l || echo 0)
N_AFTER_CODE=$(bcftools view -H "$NONCODING_BCF" | wc -l || echo 0)
log "[INFO] Counts | Input=$N_IN | After FAM/SAMPLE=$N_AFTER_FS | After coding exclusion=$N_AFTER_CODE"

# =============================
# Step 4: SEX annotation
# =============================
FINAL_BCF="$OUT_PREFIX.filtered.vep.with_asdinfo.bcf"
if [[ "$TYPE" == "inh" ]]; then
  log "Annotating SEX + PHENO -> $FINAL_BCF"
  if [[ ! -s "$SEX_MAP" ]]; then build_sex_map; fi
  if [[ ! -s "$PHENO_MAP" ]]; then build_pheno_map; fi
  require_file "$SEX_MAP"
  require_file "$PHENO_MAP"
  {
    bcftools view -h "$NONCODING_BCF" | add_or_keep_sex_pheno_headers
    bcftools view -H "$NONCODING_BCF" | annotate_sex_pheno_stream "$SAMPLE_TAG" "$SEX_MAP" "$PHENO_MAP"
  } | bcftools view --threads "$THREADS" -Ob -o "$FINAL_BCF" -
else
  log "Annotating SEX only -> $FINAL_BCF"
  if [[ ! -s "$SEX_MAP" ]]; then build_sex_map; fi
  require_file "$SEX_MAP"
  {
    bcftools view -h "$NONCODING_BCF" | add_or_keep_sex_pheno_headers
    bcftools view -H "$NONCODING_BCF" | annotate_sex_pheno_stream "$SAMPLE_TAG" "$SEX_MAP" ""   # empty PHENO map
  } | bcftools view --threads "$THREADS" -Ob -o "$FINAL_BCF" -

fi
bcftools index -f "$FINAL_BCF"

log "DONE. Final: $FINAL_BCF"

# =============================
# Cleanup
# =============================
log "Cleaning up temporary files"
rm -f "$tmp_kept_vcf" "$tmp_kept_vcf.tbi"
rm -f "$FILTERED_BCF" "$FILTERED_BCF".csi
rm -f "$NONCODING_BCF" "$NONCODING_BCF".csi
rm -f "$CACHE" "$CACHE".tbi
log "All done."
