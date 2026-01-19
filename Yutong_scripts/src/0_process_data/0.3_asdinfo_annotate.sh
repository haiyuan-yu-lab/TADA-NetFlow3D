#!/bin/bash
#SBATCH --job-name=annotate_asdinfo
#SBATCH -o logs/annotate_asdinfo.%j.out
#SBATCH -e logs/annotate_asdinfo.%j.err
#SBATCH --mail-type=ALL
#SBATCH --mail-user=yz2676@cornell.edu

#SBATCH -N 1
#SBATCH --ntasks-per-node=1
#SBATCH --cpus-per-task=1
#SBATCH --mem=128G
#SBATCH --time=480:00:00
#SBATCH --nodelist=cbsuhy02


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
# 4) For parent variants with untransmitted annotation:
# sbatch src/0_process_data/0.3_asdinfo_annotate.sh --type parent --make-untransmitted yes --exclude-multiplex no --keep-only-multiplex yes --out-prefix UCSF_parents_WGS_hg38.onlymultiplex
# sbatch src/0_process_data/0.3_asdinfo_annotate.sh --type parent --make-untransmitted yes --exclude-multiplex no --inh-for-untrans All_inherited_WGS_hg38.includesimplex --out-prefix UCSF_parents_WGS_hg38.includesimplex


set -eo pipefail
# -----------------------------
# Parameters
# -----------------------------
TYPE="dnm"  # dnm | inh | parent
VEP_ANNO_DIR="/fs/cbsuhy01/storage/yz2676/projects/1-ASD/NetFlow_ASD/Yutong_scripts/data/annotated"
DNV_VEP_REL="dnm/UCSF_20211001_dnms_WGS_hg38.with_vep.vcf.gz"
INH_VEP_REL="inh/All_inherited_WGS_hg38.with_vep.vcf.gz"
PARENT_VEP_REL="parents/UCSF_parents_WGS_hg38.with_vep.vcf.gz"

INH_FOR_UNTRANS_REL="All_inherited_WGS_hg38"

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
MAKE_UNTRANSMITTED="no"    # yes|no

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
    --parent-vep)         PARENT_VEP_REL="$2"; shift 2;;
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
    --make-untransmitted)  MAKE_UNTRANSMITTED="$2"; shift 2;;
    --inh-for-untrans)    INH_FOR_UNTRANS_REL="$2"; shift 2;;
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

has_info_tag(){
  local vcf="$1" tag="$2"
  bcftools view -h "$vcf" | grep -qE "^##INFO=<ID=${tag},"
}

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

add_parent_untrans_headers(){
  awk '
    BEGIN { add_off=1 }
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

build_multiplex_parent_samples(){
  local out="$1"
  require_file "$PHENO"
  require_file "$FAM_LIST"

  awk -F'\t' -v OFS='\t' '
    NR==1{
      for(i=1;i<=NF;i++){
        gsub(/\r$/,"",$i)
        if($i=="Fam") fam=i
        if($i=="Father") f=i
        if($i=="Mother") m=i
      }
      next
    }
    { print $fam, $f, $m }
  ' "$PHENO" \
  | awk -F'\t' -v OFS='\t' -v FL="$FAM_LIST" '
    BEGIN{
      while((getline l < FL)>0){ sub(/\r$/,"",l); if(l!="") F[l]=1 }
      close(FL)
    }
    ($1 in F){
      if($2!="0" && $2!="") print $2
      if($3!="0" && $3!="") print $3
    }
  ' | sort -u > "$out"

  log "Multiplex parent samples -> $out ($(wc -l < "$out" 2>/dev/null || echo 0))"
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

build_trio_table(){
  # output: child \t father \t mother \t child_pheno
  local out="$1"
  require_file "$PHENO"
  awk -F'\t' -v OFS='\t' '
    NR==1{
      for(i=1;i<=NF;i++){
        gsub(/\r$/,"",$i)
        if($i=="SampleID") c=i
        if($i=="Father")   f=i
        if($i=="Mother")   m=i
        if($i=="Pheno")    p=i
      }
      next
    }
    # skip incomplete trios
    ($f!="0" && $m!="0" && $c!="" && $p!="") { print $c, $f, $m, $p }
  ' "$PHENO" > "$out"
  log "Trio table -> $out ($(wc -l < "$out" 2>/dev/null || echo 0))"
}

# Build a tabix-indexed TSV keyed by INFO/SAMPLE (chrom-like), so tabix can fetch by sample id.
# Output columns:
#   1) SAMPLE_ID   2) POS   3) END   4..= original VCF columns (CHROM POS ID REF ALT QUAL FILTER INFO ...)
# tabix uses col1 as "chrom", col2/3 as "pos/end".
build_sample_indexed_tsv(){
  local in_bcf="$1"
  local out_tsv_bgz="$2"
  local sample_tag="$3"

  if [[ -s "$out_tsv_bgz" && -s "${out_tsv_bgz}.tbi" ]]; then
    log "Found existing indexed TSV: $out_tsv_bgz"
    return 0
  fi

  log "Building sample-indexed TSV -> $out_tsv_bgz"

  local tmpdir="$ASD_DIR/tmp"
  mkdir -p "$tmpdir"
  export TMPDIR="$tmpdir"

  # Write to temp output then mv (avoid half-written final file)
  local out_tmp="${out_tsv_bgz}.tmp.$(date +%s).$$"
  trap 'rm -f "$out_tmp" "${out_tmp}.tbi"' RETURN

  # Stream BCF records once and attach (sample_id, line_no, line_no) in front.
  bcftools view -H "$in_bcf" \
  | awk -v FS="\t" -v OFS="\t" -v KEY="$sample_tag" '
      function get(info,key, a){
        return match(info, "(^|;)" key "=([^;]*)", a) ? a[2] : ""
      }
      {
        sid = get($8, KEY)
        if (sid=="") next
        n++
        # tabix needs integer coordinates; use line number as position
        print sid, n, n, $0
      }
    ' | LC_ALL=C sort -T "$tmpdir" -k1,1 -k2,2n \
      | bgzip -c > "$out_tmp"

  tabix -f -s 1 -b 2 -e 3 "$out_tmp"

  mv -f "$out_tmp" "$out_tsv_bgz"
  mv -f "${out_tmp}.tbi" "${out_tsv_bgz}.tbi"

  trap - RETURN
  log "Done indexing: $out_tsv_bgz"
}

# -----------------------------
# Inputs per TYPE
# -----------------------------
case "$TYPE" in
  dnm) IN_VCF="$VEP_ANNO_DIR/$DNV_VEP_REL";;
  inh) IN_VCF="$VEP_ANNO_DIR/$INH_VEP_REL";;
  parent) IN_VCF="$VEP_ANNO_DIR/$PARENT_VEP_REL";;
  *)   echo "ERROR: --type must be dnm or inh or parent" >&2; exit 1;;
esac
require_file "$IN_VCF"

if [[ -n "$USER_PREFIX" ]]; then
    OUT_PREFIX="$ASD_DIR/$(basename "$USER_PREFIX")"
else
    case "$TYPE" in
    dnm) OUT_PREFIX="$ASD_DIR/UCSF_20211001_dnms_WGS_hg38";;
    inh) OUT_PREFIX="$ASD_DIR/All_inherited_WGS_hg38";;
    parent) OUT_PREFIX="$ASD_DIR/UCSF_parents_WGS_hg38";;
    esac
fi

log "[INFO] TYPE=$TYPE | IN_VCF=$IN_VCF"
log "[INFO] OUT_PREFIX=$OUT_PREFIX"
log "[INFO] TOGGLES: exclude_multiplex=$EXCLUDE_MULTIPLEX, exclude_concerning=$EXCLUDE_CONCERNING, exclude_coding=$EXCLUDE_CODING, keep_only_multiplex=$KEEP_ONLY_MULTIPLEX, make_untransmitted=$MAKE_UNTRANSMITTED"
log "[INFO] THREADS=$THREADS"

KEEP_SAMPLE_LIST=""   # allowlist (optional)
DROP_SAMPLE_LIST=""   # denylist (optional)

# Build lists only if that exclusion is ON and no custom list provided (or default path used)
if [[ "$EXCLUDE_MULTIPLEX" == "yes" && "$KEEP_ONLY_MULTIPLEX" == "yes" ]]; then
  log "ERROR: exclude_multiplex and keep_only_multiplex are mutually exclusive"; exit 1
fi

if has_info_tag "$IN_VCF" "$FAMILY_TAG"; then
  HAVE_FAM=1
else
  HAVE_FAM=0
fi
log "INFO tag check: ${FAMILY_TAG} present? $HAVE_FAM"

if [[ "$EXCLUDE_MULTIPLEX" == "yes" || "$KEEP_ONLY_MULTIPLEX" == "yes" ]]; then
  if [[ ! -s "$FAM_LIST" || "$FAM_LIST" == "$ASD_DIR/multiplex_families.txt" ]]; then
    build_multiplex_list
  else
    log "Using provided multiplex family list: $FAM_LIST"
  fi
fi

if [[ "$KEEP_ONLY_MULTIPLEX" == "yes" ]]; then
  if [[ "$HAVE_FAM" -eq 1 ]]; then
    # fam-based filtering will happen inside awk using FAM_LIST
    :
  else
    # No INFO/FAM, so we must use a sample allowlist derived from PHENO
    if [[ "$TYPE" == "parent" ]]; then
      MULTIPLEX_PARENT_SAMPLES="$ASD_DIR/multiplex_parent_samples.txt"
      if [[ ! -s "$MULTIPLEX_PARENT_SAMPLES" ]]; then
        build_multiplex_parent_samples "$MULTIPLEX_PARENT_SAMPLES"
      fi
      KEEP_SAMPLE_LIST="$MULTIPLEX_PARENT_SAMPLES"
      log "Using multiplex parent sample allowlist: $KEEP_SAMPLE_LIST"
    else
      log "ERROR: --keep-only-multiplex yes but INFO/$FAMILY_TAG is missing and TYPE=$TYPE is not 'parent'."
      log "       Provide a sample allowlist or add INFO/$FAMILY_TAG to the VCF."
      exit 1
    fi
  fi
fi

if [[ "$EXCLUDE_CONCERNING" == "yes" ]]; then
  if [[ ! -s "$SAMP_LIST" || "$SAMP_LIST" == "$ASD_DIR/samples_exclude.txt" ]]; then
    build_concerning_list
  else
    log "Using provided concerning sample list: $SAMP_LIST"
  fi
  DROP_SAMPLE_LIST="$SAMP_LIST"
fi

# -----------------------------
# Step 1: Exclude/keep multiplex families and concerning siblings
# -----------------------------
FILTERED_BCF="$OUT_PREFIX.filtered.bcf"
log "Applying selected exclusions -> $FILTERED_BCF"

{ bcftools view -h "$IN_VCF"
  bcftools view -H "$IN_VCF" \
  | awk -v FS="\t" -v OFS="\t" \
        -v have_fam="$HAVE_FAM" \
        -v fam_list="$FAM_LIST" \
        -v keep_list="$KEEP_SAMPLE_LIST" -v drop_list="$DROP_SAMPLE_LIST" \
        -v use_fam="$EXCLUDE_MULTIPLEX" \
        -v keep_only="$KEEP_ONLY_MULTIPLEX" \
        -v fam_tag="$FAMILY_TAG" -v samp_tag="$SAMPLE_TAG" '
  BEGIN{
    # Load fam list if we are going to use fam filtering
    if (have_fam==1 && (tolower(use_fam)=="yes" || tolower(keep_only)=="yes") &&
        system("[ -s \""fam_list"\" ]")==0){
      while ((getline l < fam_list)>0){ sub(/\r$/,"",l); if(l!="") F[l]=1 }
      close(fam_list)
    }
    # Load keep/drop sample lists if provided
    if (keep_list!="" && system("[ -s \""keep_list"\" ]")==0){
      while ((getline l < keep_list)>0){ sub(/\r$/,"",l); if(l!="") KEEP[l]=1 }
      close(keep_list); use_keep=1
    } else use_keep=0

    if (drop_list!="" && system("[ -s \""drop_list"\" ]")==0){
      while ((getline l < drop_list)>0){ sub(/\r$/,"",l); if(l!="") DROP[l]=1 }
      close(drop_list); use_drop=1
    } else use_drop=0
  }
  function get(info,key,a){ return match(info, "(^|;)"key"=([^;]*)", a) ? a[2] : "" }
  # does any comma-separated token exist in assoc array A?
  function any_in(val,A){ if(val=="") return 0; { if(val in A) return 1 } return 0 }
  {
    info = $8
    samp=get(info, samp_tag)
    if (samp=="") next

    # sample filters always apply if lists exist
    if (use_keep && !(samp in KEEP)) next
    if (use_drop && (samp in DROP)) next

    # fam filters apply only if have_fam==1
    if (have_fam==1){
      fam=get(info, fam_tag)
      if (tolower(keep_only)=="yes" && !any_in(fam, F)) next
      if (tolower(use_fam)=="yes" && any_in(fam, F)) next
    }

    print
  }'
} | bcftools view --threads "$THREADS" -Ob -o "$FILTERED_BCF" -
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
# Step 4: SEX/PHENO annotation
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
# Step 5: Optional untransmitted parent variants
# =============================
if [[ "$MAKE_UNTRANSMITTED" == "yes" && "$TYPE" != "parent" ]]; then
  log "ERROR: --make-untransmitted yes requires --type parent"; exit 1
fi

if [[ "$MAKE_UNTRANSMITTED" == "yes" && "$TYPE" == "parent" ]]; then
  log "Emitting untransmitted parent variants -> $OUT_PREFIX.untransmitted.bcf"
  TRIO_TABLE="$ASD_DIR/trio_table.txt"
  build_trio_table "$TRIO_TABLE"
  
  if [[ "$KEEP_ONLY_MULTIPLEX" == "no" ]]; then
    log "Processing all trios for untransmitted variants, use parallelization."; 
    NCHUNKS=25
    CHUNK_DIR="$ASD_DIR/trio_chunks"
    mkdir -p "$CHUNK_DIR"

    split -n l/$NCHUNKS -d -a 4 "$TRIO_TABLE" "$CHUNK_DIR/trios."
    exit 1
  fi

  UNTRANSMITTED_BCF="$OUT_PREFIX.untransmitted.bcf"
  INH_FOR_UNTRANS_BCF="$ASD_DIR/$INH_FOR_UNTRANS_REL.filtered.vep.with_asdinfo.bcf"
  require_file "$INH_FOR_UNTRANS_BCF"
  log "Using inherited variants for transmitted set: $INH_FOR_UNTRANS_BCF"
  
  # Build indexed TSVs once
  tag="$(basename "$INH_FOR_UNTRANS_BCF" .bcf)"
  TX_TSV="$ASD_DIR/${tag}.by_child.tsv.gz"

  tagp="$(basename "$FINAL_BCF" .bcf)"
  PR_TSV="$ASD_DIR/${tagp}.by_parent.tsv.gz"

  build_sample_indexed_tsv "$INH_FOR_UNTRANS_BCF" "$TX_TSV" "$SAMPLE_TAG"
  build_sample_indexed_tsv "$FINAL_BCF" "$PR_TSV" "$SAMPLE_TAG"

  # Generate VCF body with python (no header)
  BODY_TXT="$OUT_PREFIX.untransmitted.body.vcf"
  python3 src/0_process_data/make_untransmitted.py \
    --trios "$TRIO_TABLE" \
    --transmitted-tsv "$TX_TSV" \
    --parent-tsv "$PR_TSV" \
    --out-body "$BODY_TXT"
  
  {
    bcftools view -h "$FINAL_BCF" | add_parent_untrans_headers
    cat "$BODY_TXT"
  } | bcftools view --threads "$THREADS" -Ob - \
  | bcftools sort -m 8G -Ob -o "$UNTRANSMITTED_BCF" -

  bcftools index -f "$UNTRANSMITTED_BCF"
  log "DONE. Untransmitted parent variants: $UNTRANSMITTED_BCF"
fi

# =============================
# Step 6: Optional create deduplicated per variant per offspring untransmitted VCF
# =============================
if [[ "$MAKE_UNTRANSMITTED" == "yes" && "$TYPE" == "parent" ]]; then
  log "Creating deduplicated per-variant-per-offspring untransmitted VCF -> $OUT_PREFIX.untransmitted.dedup.bcf"
  UNTRANSMITTED_BCF="$OUT_PREFIX.untransmitted.bcf"
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
fi

# =============================
# Cleanup
# =============================
log "Cleaning up temporary files"
rm -f "$tmp_kept_vcf" "$tmp_kept_vcf.tbi"
rm -f "$FILTERED_BCF" "$FILTERED_BCF".csi
rm -f "$NONCODING_BCF" "$NONCODING_BCF".csi
rm -f "$CACHE" "$CACHE".tbi
rm -f "$BODY_TXT"
log "All done."
