#!/usr/bin/env Rscript
# Fit selected coding and noncoding TADA evidence from prepared tables.
suppressPackageStartupMessages({
  library(GenomicRanges)
  library(stringr)
})
rm(list=ls())

args <- commandArgs(trailingOnly = TRUE)
if (length(args) == 1L && args[[1]] %in% c("--help", "-h")) {
  cat(paste(
    "Usage: Rscript src/tada/run_tada.R [options]",
    "Required: --info-file FILE --prior-file FILE --sample-file FILE --output-file FILE",
    "Inputs: --coding-file FILE; --noncoding-file FILE --mut-col COLUMN",
    "Select evidence explicitly (all categories default to disabled):",
    "  --dn-coding --in-coding          Enable all coding categories for that component",
    "  --dn-noncoding --in-noncoding    Enable the corresponding enhancer component",
    "  --dn-ptv --in-ptv --dn-misb --in-misb --dn-misa --in-misa",
    "  --dn-em --in-em                 Individual enhancer switches",
    "  --in-floor                     Floor inherited BFs (default: off)",
    "  --help, -h                     Show this help",
    "Only hyphenated long options are accepted. Paths are relative to your working directory.",
    sep = "\n"
  ), "\n")
  quit(status = 0)
}

file_arg <- grep("^--file=", commandArgs(trailingOnly = FALSE), value = TRUE)
script_dir <- dirname(normalizePath(sub("^--file=", "", file_arg[[1]])))
source(file.path(script_dir, "tada_functions.R"))

value_options <- c(
  "--coding-file" = "coding_file", "--noncoding-file" = "noncoding_file",
  "--info-file" = "info_file", "--prior-file" = "prior_file",
  "--sample-file" = "sample_file", "--output-file" = "output_file",
  "--mut-col" = "mut_col"
)
flag_options <- c(
  "--dn-coding" = "dn_coding", "--in-coding" = "in_coding",
  "--dn-noncoding" = "dn_noncoding", "--in-noncoding" = "in_noncoding",
  "--dn-ptv" = "dn_ptv", "--in-ptv" = "in_ptv",
  "--dn-misb" = "dn_misb", "--in-misb" = "in_misb",
  "--dn-misa" = "dn_misa", "--in-misa" = "in_misa",
  "--dn-em" = "dn_em", "--in-em" = "in_em", "--in-floor" = "in_floor"
)
for (variable in unname(value_options)) assign(variable, NULL)
for (variable in unname(flag_options)) assign(variable, FALSE)
i <- 1L
while (i <= length(args)) {
  option <- args[[i]]
  if (option %in% names(value_options)) {
    if (i == length(args) || startsWith(args[[i + 1L]], "--")) {
      stop("Missing value for ", option, call. = FALSE)
    }
    # Match the original get_next_arg(): the first occurrence supplies the value.
    if (is.null(get(value_options[[option]]))) {
      assign(value_options[[option]], args[[i + 1L]])
    }
    i <- i + 2L
  } else if (option %in% names(flag_options)) {
    assign(flag_options[[option]], TRUE)
    i <- i + 1L
  } else {
    stop("Unsupported argument: ", option, ". Use --help; underscore options are not accepted.", call. = FALSE)
  }
}

# group flags
if (dn_coding) dn_ptv <- dn_misb <- dn_misa <- TRUE
if (in_coding) in_ptv <- in_misb <- in_misa <- TRUE
if (dn_noncoding) dn_em <- TRUE
if (in_noncoding) in_em <- TRUE

### check required inputs
coding_enabled <- any(dn_ptv, dn_misb, dn_misa, in_ptv, in_misb, in_misa)
if (coding_enabled && is.null(coding_file))
  stop("Coding analysis enabled but --coding-file is missing.")

noncoding_enabled <- any(dn_em, in_em)
if (noncoding_enabled && is.null(noncoding_file))
  stop("Noncoding analysis enabled but --noncoding-file is missing.")
if (noncoding_enabled && is.null(mut_col))
  stop("Noncoding analysis enabled but --mut-col is missing.")
if (!coding_enabled && !noncoding_enabled)
  stop("No analysis enabled. Use coding and/or noncoding flags.")

if (coding_enabled) {
  message("Input coding file: ", coding_file)
}
if (noncoding_enabled) {
  message("Input noncoding file: ", noncoding_file)
}

# Required shared tables and destination; coding/noncoding requirements above depend on selection.
for (option in c("--info-file", "--prior-file", "--sample-file", "--output-file")) {
  if (is.null(get(value_options[[option]]))) stop("Missing required option: ", option, call. = FALSE)
}
for (input in c(info_file, prior_file, sample_file, coding_file, noncoding_file)) {
  if (!file.exists(input)) stop("Input file not found: ", input, call. = FALSE)
}

### load data
info <- read.table(info_file, header = TRUE, sep = "\t")
for (i in c(prior_file, sample_file)) {
  df <- read.table(i, header = TRUE, sep = "\t")
  for (j in 1:nrow(df)) {
    assign(df$variable[j], df$value[j], envir = .GlobalEnv)
  }
}

if (!is.null(coding_file)) {
  coding <- read.table(coding_file, header = TRUE, sep = "\t")
} else {
  message("No coding file provided; using info file as gene backbone.")
  coding <- info
}
if (!"Hugo_Symbol" %in% names(coding)) {
  if ("gene" %in% names(coding)) {
    coding$Hugo_Symbol <- coding$gene
  } else if ("gene_gencodeV33" %in% names(coding)) {
    coding$Hugo_Symbol <- coding$gene_gencodeV33
  } else {
    stop("Coding file is missing a gene symbol column: expected Hugo_Symbol, gene, or gene_gencodeV33.")
  }
}
if (!is.null(noncoding_file)) {
  noncoding <- read.table(noncoding_file, header = TRUE, sep = "\t")
}

# ### estimate prior
# pred = proportion_risk_gene(tada=coding, col="ptv", N=n_dn_case_coding, interval=c(2, 50))
# # proportion_risk_gene(tada=coding, col="misb", N=n_dn_case_coding, interval=c(2, 15))
# # proportion_risk_gene(tada=coding, col="misa", N=n_dn_case_coding, interval=c(2, 15))
# # proportion_risk_gene(tada=noncoding, col="em", N=n_dn_case_noncoding, interval=c(1, 15))

### coding analysis
for (var_type in c("ptv", "misb", "misa")) {
  var_name <- paste0("dn_", var_type)
  if (get(var_name)) {
    message(paste(var_name, "analysis is enabled."))
    coding[[paste0("bf.dn.", var_type)]] <- BF_DN(
      count_case = coding[[paste0("dn.", var_type)]],
      count_con = coding[[paste0("dn.", var_type, ".sib")]],
      n_case = n_dn_case_coding,
      n_con = n_dn_control_coding,
      mut = coding[[paste0("mut.", var_type)]],
      gamma.dn = coding[[paste0("prior.dn.", var_type)]],
      beta.dn = beta.dn
    )
  }
  var_name = paste0("in_", var_type)
  if (get(var_name)) {
    message(paste(var_name, "analysis is enabled."))
    coding[[paste0("bf.in.", var_type)]] <- BF_CC(
      count_case = coding[[paste0("in.t.", var_type)]],
      count_con = coding[[paste0("in.u.", var_type)]],
      n_case = n_in_case_coding,
      n_con = n_in_case_coding,
      mut = coding[[paste0("mut.", var_type)]],
      gamma.cc = coding[[paste0("prior.in.", var_type)]],
      floor = in_floor
    )
    ## combine de novo + inherited BF
    if (get(paste0("dn_", var_type)) && get(paste0("in_", var_type))) {
      # coding[[paste0("bf.", var_type)]] <- pmax(coding[[paste0("bf.dn.", var_type)]] * coding[[paste0("bf.in.", var_type)]], 1)
      coding[[paste0("bf.", var_type)]] <- coding[[paste0("bf.dn.", var_type)]] * coding[[paste0("bf.in.", var_type)]]
    }
  }
}

### noncoding analysis
if (!is.null(noncoding_file)) {
  noncoding[[mut_col]] <- as.numeric(na_if(noncoding[[mut_col]], "."))
  noncoding <- noncoding[!is.na(noncoding[[mut_col]]), ]
  ## de novo
  if(dn_em) {
    message("dn.em analysis is enabled.")
    noncoding$bf.dn.em <- BF_DN(
      count_case=noncoding$dn.em,
      count_con=noncoding$dn.em.sib,
      n_case=n_dn_case_noncoding,
      n_con=n_dn_control_noncoding,
      mut=noncoding[[mut_col]],
      gamma.dn=rep(gamma.dn.em, times = nrow(noncoding)),
      beta.dn=beta.dn
    )
    lBF <- do.call(rbind, lapply(1:nrow(noncoding), function(i) split_weight(log(noncoding$bf.dn.em[i]), noncoding$gene_id[i], info)))
    lBF <- na.omit(lBF)
    BF <- exp(by(lBF$lBF, lBF$gene, sum))
    df_output <- data.frame(gene_id = names(BF), bf.dn.em = as.vector(BF))
    coding <- merge(coding, df_output, by = "gene_id", all.x = TRUE)
    coding$bf.dn.em[is.na(coding$bf.dn.em)] <- 1
  }

  ## inherited
  if (in_em){
    message("in.em analysis is enabled.")
    noncoding$bf.in.em <- BF_CC(
      count_case=noncoding$in.t.em,
      count_con=noncoding$in.u.em,
      n_case=n_in_case_noncoding,
      n_con=n_in_case_noncoding,
      mut=noncoding[[mut_col]],
      gamma.cc=rep(gamma.in.em, times = nrow(noncoding)),
      floor=FALSE
    )
    lBF <- do.call(rbind, lapply(1:nrow(noncoding), function(i) split_weight(log(noncoding$bf.in.em[i]), noncoding$gene_id[i], info)))
    lBF <- na.omit(lBF)
    BF <- exp(by(lBF$lBF, lBF$gene, sum))
    df_output <- data.frame(gene_id = names(BF), bf.in.em = as.vector(BF))
    if (in_floor) {
      df_output$bf.in.em = pmax(df_output$bf.in.em, 1) # floor BF at 1
    }
    coding <- merge(coding, df_output, by = "gene_id", all.x = TRUE)
    coding$bf.in.em[is.na(coding$bf.in.em)] <- 1
  }

  ## combine de novo + inherited enhancer BF
  if (dn_em && in_em) {
    # coding$bf.em <- pmax(coding$bf.dn.em * coding$bf.in.em, 1)
    coding$bf.em <- coding$bf.dn.em * coding$bf.in.em
  }
}

## TADA identifications
bf_cols <- grep("^bf\\.(dn|in)\\.", names(coding), value = TRUE)
coding$bf <- apply(coding[, bf_cols, drop = FALSE], 1, function(x) prod(x, na.rm = TRUE))
coding$fdr <- Bayesian.FDR(coding$bf, pi0=pi0)
qval <- setNames(coding$fdr, coding$gene_id)
coding$pval <- q2p(qval)
bf_cols <- grep("^bf", names(coding), value=TRUE)
write.table(
  coding[, c("gene_id", "Hugo_Symbol", "LOEUF", "pLI", bf_cols, "fdr", "pval")],
  file = output_file,
  sep = "\t",
  row.names = FALSE,
  quote = FALSE
)
