# Copyright 2022 
# Original code and conception: Xin He
# Modifications and extensions: Jiebiao Wang, Lambertus Klei, and Jack Fu
# Distributed under terms of the MIT license.
### Change to working directory where the folder was cloned
setwd('/fs/cbsuhy01/storage/NetFlow_ASD/Week8_429/Jinjin')
source('./functions.R')
source('./prior_PTV_LOEUF.R')

run_noncode<-function(
    dn_n_control = 2938, 
    dn_n_case = 4467,
    dn_name = 'all', # one from c('none','all','female','male'),
    
    trans_n_control = 3103, 
    trans_n_case = 4794,
    trans_name = 'none', # one from c('none', 'all','female','male'), NULL means don't consider
    
    untrans_n_control = 3103, 
    untrans_n_case = 4794,
    untrans_name = 'none', # one from c('none', 'all','female','male'), NULL then don't consider
    
    mutability_name = '20230320anno', # one from c('20230320anno', 'samocha14' ,'MR', 'MC')
    
    withFu = FALSE # whether combine fu's bf or not
){
  
  mydata <- read.delim('./TADA_input_scaled_withsib.txt', header=TRUE)
  mydata = na.omit(mydata)
  
  # match the gene id
  mydata$gene_id = sapply(mydata$gene_id, function(x){
    strsplit(x,split='[.]')[[1]][1]
  })
  
  # intersect the data to get louef score
  data_loufe <- read.xlsx('./fu_Table.xlsx', sheet = 5)
  geneid = intersect(mydata$gene_id, data_loufe$gene_id)

  # get the data input for bayes factor calculation

  data = data.frame(gene = data_loufe$gene[match(geneid, data_loufe$gene_id)],
                    geneid = geneid,
                    mut = mydata[[paste0("Noncoding_Mutability_",mutability_name)]][match(geneid, mydata$gene_id)])
                    
  if(dn_name!='none'){
    gamma = getprior(mut_name = mutability_name, dn_name = dn_name, dn_n = dn_n_case, info = 'dn', plot = FALSE)
    data[["dn_case"]] = mydata[[paste0("Noncoding_DNV_",dn_name,"_case")]][match(geneid, mydata$gene_id)]
    data[["dn_control"]] = mydata[[paste0("Noncoding_DNV_",dn_name,"_control")]][match(geneid, mydata$gene_id)]
    data[['gamma']] = gamma
  }
  if(trans_name!='none'){
    gammatrans = getprior(mut_name = mutability_name, dn_name = trans_name, dn_n = trans_n_case, info = 'trans', plot = FALSE)
    data[['trans_case']] = mydata[[paste0("Noncoding_Transmitted_",trans_name,"_case")]][match(geneid, mydata$gene_id)]
    data[['trans_control']] = mydata[[paste0("Noncoding_Transmitted_",trans_name,"_control")]][match(geneid, mydata$gene_id)]
    data[['trans_gamma']] = gammatrans
  }
  if(untrans_name!='none'){
    gammauntrans = getprior(mut_name = mutability_name, dn_name = untrans_name, dn_n = untrans_n_case, info = 'untrans', plot = FALSE)
    data[['untrans_case']] = mydata[[paste0("Noncoding_Untransmitted_",untrans_name,"_case")]][match(geneid, mydata$gene_id)]
    data[['untrans_control']] = mydata[[paste0("Noncoding_Untransmitted_",untrans_name,"_control")]][match(geneid, mydata$gene_id)]
    data[['untrans_gamma']] = gammauntrans
  }
  
  #=====compute bayes factor
  # for deno
  if(dn_name!='none'){
    BF_dn = BF_DN_SNV(count_case = data$dn_case, count_con = data$dn_control, 
                    n_case = dn_n_case, n_con = dn_n_control, 
                    mut = data$mut, gamma.dn = data$gamma, 
                    beta.dn = 0.2)
  }else{
    BF_dn = rep(1,length(data$gene))
  }
  
  if(trans_name!='none'){
    BF_trans = BF_DN_SNV(count_case = data$trans_case, count_con = data$trans_control, 
                         n_case = trans_n_case, n_con = trans_n_control, 
                         mut = data$mut, gamma.dn = data$gamma, 
                         beta.dn = 0.2)
  }else{
    BF_trans = rep(1,length(data$gene))
  }
  
  if(untrans_name!='none'){
    BF_untrans = BF_DN_SNV(count_case = data$untrans_case, count_con = data$untrans_control, 
                           n_case = untrans_n_case, n_con = untrans_n_control, 
                           mut = data$mut, gamma.dn = data$gamma, 
                           beta.dn = 0.2)
  }else{
    BF_untrans = rep(1,length(data$gene))
  }
  
  BF_asd = cbind(BF_dn, BF_trans, BF_untrans)
  
  if(withFu){
    fudata = openxlsx::read.xlsx('./fu_Table.xlsx', sheet=11)
    BF_asd = cbind(BF_asd,
                   bf_code_dd_ptv  = pmax(10**fudata$l10BF_DD_PTV[match(data$gene,fudata$gene)], 1), 
                   bf_code_dd_misa  = pmax(10**fudata$l10BF_DD_misA[match(data$gene,fudata$gene)], 1), 
                   bf_code_dd_mibs  = pmax(10**fudata$l10BF_DD_misB[match(data$gene,fudata$gene)], 1))
  }  
  
  BF_asd = na.omit(BF_asd)
  
  # multiply all the bayes factors 
  BF = apply(BF_asd, 1, prod)
  
  # get q-value
  qval <- Bayesian.FDR(BF, pi0 = 1 - 0.06)
  names(qval) = data$gene
  
  # TADA-ASD significant genes at various thresholds
  tada_genes = sapply(c(0.001, 0.01, 0.05, 0.1), function(x) sum(qval<x))
  print(paste0('get ', tada_genes, ' tada genes at q-val cutoff ', c(0.001, 0.01, 0.05, 0.1), '.'))
  
  # back-transform FDR to p-values
  ec_threshold <- .05/length(qval)
  pval <- q2p(qval)
  print(paste0('get ', sum(pval <= ec_threshold), ' tada genes at p-val cutoff ', ec_threshold, '.'))
  pval[is.infinite(pval)] = 1
  
  results = data.frame(gene = data$gene,
                       geneid = data$geneid,
                       FDR = qval,
                       pval = pval,
                       log10_BF = log10(BF))
  
  
  filename = paste0('./results_mut-',mutability_name)
  if(dn_name!='none'){
    filename = paste0(filename,'_dn-',dn_name)
    results$gamma = data$gamma
  }
  if(trans_name!='none'){
    filename = paste0(filename, '_trans-', trans_name)
    results$gamma_trans = data$gamma_trans
  }
  if(untrans_name!='none'){
    filename = paste0(filename, '_untrans-', untrans_name)
    results$gamma_untrans = data$gamma_untrans
  }
  if(withFu)filename = paste0(filename, '_withFu')
  write.csv(results, file=paste0('./', filename, '.csv'), row.names=FALSE)
}

#========command line
# Rscript ./run_noncode.R $dn-name $trans-name $untrans-anme, $mut_name, $withfu
args=(commandArgs(trailingOnly = TRUE))
run_noncode(dn_name =  as.character(args[1]), 
            trans_name = as.character(args[2]), 
            untrans_name = as.character(args[3]), 
            mutability_name = as.character(args[4]), 
            withFu = as.logical(args[5]))


