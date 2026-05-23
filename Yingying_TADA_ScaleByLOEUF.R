library(GenomicRanges); library(stringr); 
rm(list=ls())
source('/Genomics/argo/users/yz8667/Autism_network_integration/scripts/TADA/Yingying_functions.R')

### user input
# Initialize variables
args <- commandArgs(trailingOnly = TRUE)
coding_file <- NULL
noncoding_file <- NULL
output_file <- NULL
prior_file <- NULL
info_file <- NULL
sample_file <- NULL
dn_ptv <- TRUE
in_ptv <- TRUE
dn_misb <- TRUE
in_misb <- TRUE
dn_misa <- TRUE
in_misa <- TRUE
dn_em <- FALSE
in_em <- FALSE
in_floor <- FALSE

# Process each flag and its value
for (i in seq_along(args)) {
  switch(args[i],
         "--coding_file" = { coding_file <- get_next_arg("--coding_file", args) },
         "--noncoding_file" = { noncoding_file <- get_next_arg("--noncoding_file", args) },
         "--output_file" = { output_file <- get_next_arg("--output_file", args) },
         "--prior_file" = { prior_file <- get_next_arg("--prior_file", args) },
         "--info_file" = { info_file <- get_next_arg("--info_file", args) },
         "--sample_file" = { sample_file <- get_next_arg("--sample_file", args) },
         # "--dn_ptv" = { dn_ptv <- TRUE },
         # "--in_ptv" = { in_ptv <- TRUE },
         # "--dn_misb" = { dn_misb <- TRUE },
         # "--in_misb" = { in_misb <- TRUE },
         # "--dn_misa" = { dn_misa <- TRUE },
         # "--in_misa" = { in_misa <- TRUE },
         "--dn_em" = { dn_em <- TRUE },
         "--in_em" = { in_em <- TRUE },
         "--in_floor" = { in_floor <- TRUE }
  )
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
}
if (!is.null(noncoding_file)) {
  noncoding <- read.table(noncoding_file, header = TRUE, sep = "\t")
}

### estimate prior
# proportion_risk_gene(tada=coding, col="ptv", N=n_dn_case_coding, interval=c(2, 50))
# proportion_risk_gene(tada=coding, col="misb", N=n_dn_case_coding, interval=c(2, 15))
# proportion_risk_gene(tada=coding, col="misa", N=n_dn_case_coding, interval=c(2, 15))
# proportion_risk_gene(tada=noncoding, col="em", N=n_dn_case_noncoding, interval=c(1, 15))

### coding analysis
if (!is.null(coding_file)) {
  for (var_type in c("ptv", "misb", "misa")) {
    var_name <- paste0("dn_", var_type)
    if (get(var_name)) {
      print(paste(var_name, "analysis is enabled."))
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
      print(paste(var_name, "analysis is enabled."))
      coding[[paste0("bf.in.", var_type)]] <- BF_CC(
        count_case = coding[[paste0("in.t.", var_type)]],
        count_con = coding[[paste0("in.u.", var_type)]],
        n_case = n_in_case_coding,
        n_con = n_in_case_coding,
        mut = coding[[paste0("mut.", var_type)]],
        gamma.cc = coding[[paste0("prior.in.", var_type)]],
        floor=in_floor
      )
    }
  }
} else {
  coding <- info
}


### noncoding analysis
if (!is.null(noncoding_file)) {
  ## de novo
  if(dn_em) {
    print("dn.em analysis is enabled.")
    noncoding$bf.dn.em <- BF_DN(count_case=noncoding$dn.em, count_con=noncoding$dn.em.sib, n_case=n_dn_case_noncoding, n_con=n_dn_control_noncoding, mut=noncoding$mut.em, gamma.dn=rep(gamma.dn.em, times = nrow(noncoding)), beta.dn=beta.dn)
    lBF <- do.call(rbind, lapply(1:nrow(noncoding), function(i) split_weight(log(noncoding$bf.dn.em[i]), noncoding$gene_id[i], info)))
    lBF <- na.omit(lBF)
    BF <- exp(by(lBF$lBF, lBF$gene, sum))
    df_output <- data.frame(gene_id = names(BF), bf.dn.em = as.vector(BF))
    coding <- merge(coding, df_output, by = "gene_id", all.x = TRUE)
    coding$bf.dn.em[is.na(coding$bf.dn.em)] <- 1
  }
  
  ## inherited
  if (in_em){
    print("in.em analysis is enabled.")
    noncoding$bf.in.em <- BF_CC(
      count_case=noncoding$in.t.em, 
      count_con=noncoding$in.u.em, 
      n_case=n_in_case_noncoding, 
      n_con=n_in_case_noncoding, 
      mut=noncoding$mut.em, 
      gamma.cc=rep(gamma.in.em, times = nrow(noncoding)),
      floor=FALSE)
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
}

## TADA identifications
bf_cols <- grep("^bf", names(coding))
coding$bf <- apply(coding[, bf_cols], 1, function(x) prod(x, na.rm = TRUE))
coding$fdr <- Bayesian.FDR(coding$bf, pi0=pi0)
qval <- setNames(coding$fdr, coding$gene_id)
coding$pval <- q2p(qval)
bf_cols <- grep("^bf", names(coding), value=TRUE)
write.table(coding[, c("gene_id", "gene", "LOEUF", "pLI", bf_cols, "fdr", "pval")], file = output_file, sep = "\t", row.names = FALSE)

