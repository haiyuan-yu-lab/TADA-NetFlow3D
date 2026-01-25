suppressPackageStartupMessages({
    library(stringr)
    library(dplyr)
});

get_next_arg <- function(flag, args) {
  flag_index <- match(flag, args)
  if (!is.na(flag_index) && flag_index < length(args)) {
    return(args[flag_index + 1])
  }
  return(NULL)
}

proportion_risk_gene <- function(tada, col, cons_interval=NULL, N, thresh=1, interval=c(2, 15), gamma_adj=1, gforce=NULL){
  # thresh: used to define multi-hit genes
  ind <- which(!is.na(tada[[paste("mut", col, sep = ".")]]))
  p1 <- 0.05
  ghat <- get_gamma(tada, col, ind, N, frac=p1)
  # ind <- which(tada[[paste("dn", col, sep = ".")]] > 0)
  # quant <- quantile(tada$LOEUFu, cons_interval, na.rm=TRUE)
  # ind <- which(tada$LOEUFu>=quant[1] & tada$LOEUFu<=quant[2])
  # message(length(ind), "\t", sum(tada$dn.ptv[ind]), "\t", sum(tada$dn.ptv[ind]>thresh))
  
  #if(length(gforce)==0){
  #  ghat <- uniroot(eq, interval=interval, tada=tada, col=col, ind=ind, N=N, thresh=thresh, gamma_adj=gamma_adj)$root
  #}else{
    #ghat <- gforce
  #}
  #p1 <- equation29(ghat, tada=tada, ind=ind, N=N)/length(ind) # proportion of disease genes
  # return(c(p1, as.vector(quant), mean(quant), length(ind), cons_interval, ghat))
  return(c(p1, ghat))
}

get_gamma <- function(tada, col, ind, N, frac){
  # N: Total number of probands
  C <- sum(tada[[paste("dn", col, sep = ".")]][ind]) # Total number of DN PTVs (in proband)
  mu <- tada[[paste("mut", col, sep = ".")]][ind] # Mutation rate of each gene
  m <- length(ind) # Number of genes
  return((C/(2*N*mean(mu))-m)/(frac*m)+1)
}

equation29 <- function(g, tada, col, ind, N){
  # N: Total number of probands
  C <- sum(tada[[paste("dn", col, sep = ".")]][ind]) # Total number of DN PTVs (in proband)
  mu <- tada[[paste("mut", col, sep = ".")]][ind] # Mutation rate of each gene
  m <- length(ind) # Number of genes
  return((C/(2*N*mean(mu))-m)/(g-1)) # k
}

equation33 <- function(g, tada, col, ind, N, thresh=1, gamma_adj=1){
  # N: Total number of probands
  M <- sum(tada[[paste("dn", col, sep = ".")]][ind]>thresh) # Number of multi-hit genes
  m <- length(ind) # Number of genes
  mu <- tada[[paste("mut", col, sep = ".")]][ind] # Mutation rate of each gene
  # Expected number of genes with more than *thresh* count under the null model (RR==1)
  B <- mean(1-ppois(thresh, 2*N*mu))
  # Expected number of genes with more than *thresh( count under the alt model (RR>1 ~ gamma(A,B)))
  A <- mean(sapply(ind, function(x) integrate(f, lower=0, upper=1000, g=g, mui=tada[[paste("mut", col, sep = ".")]][x], N=N, adj=gamma_adj, thresh=thresh)[[1]]))
  return((M-m*B)/(A-B)) # k
}

f <- function(gi, g, mui, N, adj, thresh){
  (1-ppois(thresh, 2*N*gi*mui))*dgamma(gi, g/adj, 1/adj)  # adj: 1 / beta; g: gamma bar
}

eq <- function(g, tada, col, ind, N, thresh, gamma_adj){
  return(equation29(g=g, tada=tada, col=col, ind=ind, N=N) - equation33(g=g, tada=tada, col=col, ind=ind, N=N, thresh=thresh, gamma_adj=gamma_adj))
}

################################################################################################
### Helper function to bayes.factor.cc
### Calculates the evidence of case-control/inherited counts under null model
## Input: 
# x.cc is a vector of case/transmitted and control/untransmitted variant counts in this gene
# n.cc is the number of case and control samples (or number of probands when evaluating inherited variants)
# Prior distribution of q|H0: Gamma(rho0, nu0)
## Output:
# cn: marginal likelihood of control/untransmitted variant counts under null model
# ca: marginal likelihood of case/transmitted variant counts under null model
################################################################################################
evidence.null.cc <- function(x.cc, n.cc, rho0, nu0) {
  marglik0.ctrl.log <- log(dnbinom(x.cc$cn, rho0, nu0/(nu0+n.cc$cn)))
  marglik0.case.log <- log(dnbinom(x.cc$ca, rho0+x.cc$cn, (nu0+n.cc$cn)/(nu0+n.cc$cn+n.cc$ca)))
  marglik0.log <- marglik0.ctrl.log + marglik0.case.log
  
  return (list(cn=exp(marglik0.ctrl.log), ca=exp(marglik0.case.log)))#, total=exp(marglik0.log)))
}

################################################################################################
### Helper function to bayes.factor.cc
### Calculates the evidence of case-control/inherited counts under alternate model
## Input: 
# x.cc is a vector of case/transmitted and control/untransmitted variant counts in this gene
# n.cc is the number of case and control samples (or number of probands when evaluating inherited variants)
# gamma.cc is a vector of the prior mean on risk
## Output:
# cn: marginal likelihood of control/untransmitted variant counts under alternate model
# ca: marginal likelihood of case/transmitted variant counts under alternate model
################################################################################################
evidence.alt.cc <- function(x.cc, n.cc, gamma.cc, rho1, nu1){
  marglik1.ctrl <- dnbinom(x.cc$cn, rho1, nu1/(nu1+n.cc$cn))
  marglik1.case <- dnbinom(x.cc$ca, rho1+x.cc$cn, (nu1+n.cc$cn)/(nu1+n.cc$ca*gamma.cc+n.cc$cn))
  marglik1 <- marglik1.ctrl * marglik1.case
  return (list(cn=marglik1.ctrl, ca=marglik1.case))#, total=marglik1))
}

################################################################################################
### Bayes factor calculation for case-control/inherited contribution to a gene
## Input: 
# x.cc is a vector of case/transmitted and control/untransmitted variant counts in this gene
# n.cc is the number of case and control samples (or number of probands when evaluating inherited variants)
# gamma.cc is a vector of the prior mean on risk
# Prior distribution of q|H1: Gamma(rho1, nu1)
# Prior distribution of q|H0: Gamma(rho0, nu0)
## Output:
# BF: Bayes factor evidence for case-control/inherited contribution to a gene
################################################################################################
bayes.factor.cc <- function(x.cc, n.cc, gamma.cc, rho1, nu1, rho0, nu0){
  marglik0.cc <- evidence.null.cc(x.cc, n.cc, rho0, nu0)
  marglik1.cc <- evidence.alt.cc(x.cc, n.cc, gamma.cc, rho1, nu1)
  BF.cn <- marglik1.cc$cn / marglik0.cc$cn
  ## Assuming null and alt of the control model is identical
  if(is.na(BF.cn)){BF.cn<-1} 
  BF.ca <- marglik1.cc$ca / marglik0.cc$ca
  BF <- BF.cn * BF.ca
  return(BF=BF)
}

################################################################################################
### log Bayes factor calculation for de novo variant contribution to a gene/TRE
## Input:
# x.dn: the de novo variant count 
# n.dn: the sample size (number of trios)
# mu: the mutation rate (of this type of mutational events in this gene/TRE)
# Prior distribution of risk: gamma ~ Gamma(gamma.dn*beta.dn, beta.dn)
## Output:
# lBF: log Bayes factor evidence for de novo variants of this type in this gene/TRE
################################################################################################
log.bayes.factor.dn <- function(x.dn, n.dn, mu,gamma.dn, beta.dn){
  marg.lik0 <- dpois(x.dn, 2*n.dn*mu, log=TRUE)
  marg.lik1 <- dnbinom(x.dn, gamma.dn*beta.dn, beta.dn/(beta.dn+2*n.dn*mu), log=TRUE)
  lBF <- marg.lik1-marg.lik0
  return (lBF=lBF)
}

################################################################################################
### Wrapper function to obtain de novo Bayes Factors for SNV/indel/DEM, over a list of genes/TREs for a given variant type
### Makes calls to log.bayes.factor.dn for each gene/TRE
## Input: 
# count_case: a vector of de novo variant counts over each gene/TRE in cases
# count_con: a vector of de novo variant counts over each gene/TRE in controls
# n_case: the number of cases
# n_con: the number of controls
# mut: a vector of mutation rates for each gene/TRE
# Prior distribution of RR: gamma ~ Gamma(gamma.dn*beta.dn, beta.dn)
## Output:
# BF: a vector of Bayes factors of de novo evidence from SNV/indel/DEM for each gene/TRE for a given variant type
################################################################################################
BF_DN <- function(count_case, count_con, n_case, n_con, mut, gamma.dn, beta.dn){
  ### If there are siblings, evaluate BF evidence for siblings to offset evidence for probands
  if(n_con>0){
    BF_con <- exp(sapply(1:length(mut), function(i) log.bayes.factor.dn(x.dn=count_con[i], n.dn=n_con, mu=mut[i], beta.dn=beta.dn, gamma.dn=gamma.dn[i])))
  }else{
    BF_con <- rep(1, length(mut))
  }
#   BF_con[BF_con<1] <- 1
  
  ### BF evidence calculation for probands
  BF_case <- exp(sapply(1:length(mut), function(i) log.bayes.factor.dn(x.dn=count_case[i], n.dn=n_case, mu=mut[i], beta.dn=beta.dn, gamma.dn=gamma.dn[i])))
#   BF_case[BF_case<1] <- 1
  
  ### Combine proband and sibling BF evidence
  ### Enforce BF back to 1 if no observations in probands or siblings
  BF <- pmax(BF_case/BF_con, 1)
  BF[which(count_case==0 & count_con==0)] <- 1
  BF[is.na(BF)] <- 1
  return(BF)
}

################################################################################################
### Wrapper function to obtain case-control/inherited Bayes Factors for SNV/indel/DEM for a given variant type
### Makes calls to bayes.factor.cc for each gene/TRE
## Input: 
# count_case: a vector of case/transmitted variant counts over each gene/TRE
# count_con: a vector of control/untransmitted variant counts over each gene/TRE
# n_case: the number of cases
# n_con: the number of controls (case-control analysis) or number of cases (inherited analysis)
# mut: a vector of mutation rates for each gene/TRE
# gamma.cc: a vector/number of the prior mean on risk
# nu: a hyperparameter that controls nu0 and nu1 in determining prior distributions, passed to bayes.factor.cc
## Output:
# BF: a vector of Bayes factors of case-control/inherited evidence from SNV/indel/DEM for each gene/TRE for a given variant type
################################################################################################
BF_CC <- function(count_case, count_con, n_case, n_con, mut, gamma.cc, nu=5000){
  rho.in = nu * sum(count_con, na.rm=T) / (2*length(mut)*n_con)
  rho.in = as.numeric(rho.in)*mut/mean(mut, na.rm=TRUE)
  BF = sapply((1:length(mut)), function(i){
    bayes.factor.cc(x.cc = data.frame(ca=count_case[i], cn=count_con[i]), n.cc=data.frame(ca=n_case, cn=n_con), gamma.cc=gamma.cc[i], rho1=rho.in[i], nu1=nu, rho0=rho.in[i], nu0=nu)
    })

  ### Enforce BF back to 1 if no mutation rate or no observed variant counts
  BF[BF<1] <- 1
  BF[mut==0] <- 1
  BF[which((count_case+count_con)==0)] <- 1
  return(BF)
}

################################################################################################
### Helper function to distribute BF evidence inversely to LOEUF score for CNV evidence
## Input: 
# x: a vector of log Bayes factors evaluated for CNV events, with names labeled by a "," concatenated string of constituent genes
# info: matrix that contains gene level information (Supplementary Table 8: gene names, constraint)
## Output:
# evidence: a matrix containing the log Bayes (lBF) factor contributed to each gene(gene)from evaluated CNVs
################################################################################################
split_weight <- function(lBF, target, info){
  genes <- unlist(str_split(target, ","))
  oes <- info$LOEUF[match(genes, info$gene_id)]
  weights <- sum(oes, na.rm = TRUE)/oes
#   message("summed weights: ", sum(weights, na.rm = TRUE))
  weights <- weights/sum(weights, na.rm = TRUE)
  return(evidence=data.frame(gene=genes, weight=weights, lBF=lBF*weights))
}

################################################################################################
### Function to transform total BF evidence to FDR values
### Bayesian FDR control (PMID:19822692, Section2.3)
## Input: 
# BF: total Bayes factor matrix (each row is a gene, multiple columns allowed)
# pi0: estimated proportion of non-risk genes
## Output:
# FDR: a vector of FDR values corresponding to the rows of the input matrix
################################################################################################
Bayesian.FDR <- function(BF, pi0) {
  # order the BF in decreasing order, need to retain order to get results back in proper order 
  i.order=order(BF, decreasing = T)
  BF=BF[i.order]
  # convert BFs to PPA (posterior probability of alternative model)
  pi <- 1-pi0
  q <- pi*BF/(1-pi+pi*BF) # PPA
  q0 <- 1 - q # posterior probability of null model
  
  # the FDR at each PPA cutoff
  FDR=cumsum(q0)/(1:length(BF))
  
  # reorder to the original order
  FDR[i.order]=FDR
  
  return (FDR=FDR)
}

################################################################################################
### Function to backtransform FDR to be on p-value scale
## Input: 
# x: a vector of FDRs
## Output:
# p: a vector of p-values
################################################################################################
q2p <- function(x){
    ec <- x
    ec <- sort(ec)
    ec <- ec*(1:length(ec))/(length(ec))/max(ec)    
    return(p=ec[names(x)])
}




