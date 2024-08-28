rm(list=ls())
setwd('/fs/cbsuhy01/storage/NetFlow_ASD/Week8_429/Jinjin')
source('./functions.R')

########################################################################
### our data
########################################################################
getprior<-function(mut_name, dn_name, dn_n, info, plot = TRUE){
  if(!file.exists(paste0('./prior_',info,'-',dn_name,'_mut-', mut_name,'.rds'))){
    data_loufe <- read.xlsx('./fu_Table.xlsx', sheet = 5)
    str(data_loufe)
    loufe = data.frame(genename = data_loufe$gene, gene_id = data_loufe$gene_id,
                       LOEUF = data_loufe$LOEUF,
                       dncount = data_loufe$dn.ptv,
                       mutation = data_loufe$mut.ptv)
    loufe = na.omit(loufe)
    
    mydata <- read.delim('./TADA_input_scaled_withsib.txt', header=TRUE)
    str(mydata)
    # match the gene id
    mydata$gene_id = sapply(mydata$gene_id, function(x){
      strsplit(x,split='[.]')[[1]][1]
    })
    # intersect the data
    geneid = intersect(mydata$gene_id, loufe$gene_id)
    
    if(info=='dn')type = 'DNV'
    if(info=='trans')type = 'Transmitted'
    if(info=='untrans')type = 'Untransmitted'
    
    # get the combined data
    data = data.frame(gene = loufe$genename[match(geneid, loufe$gene_id)],
                      geneid = geneid,
                      LOEUF = loufe$LOEUF[match(geneid, loufe$gene_id)],
                      mut.ptv = mydata[[paste0("Noncoding_Mutability_", mut_name)]][match(geneid, mydata$gene_id)],
                      dn.ptv = mydata[[paste0("Noncoding_",type,"_",dn_name,"_case")]][match(geneid, mydata$gene_id)])
    
    data = na.omit(data)
    str(data)
    
    data$LOEUFu <- data$LOEUF/2
    
    nbins=100
    spread=20
    sliding_p_non <- NULL
    for(i in 0:(nbins-spread)){
      res <- tryCatch(proportion_risk_gene(tada=data, cons_interval=c(i/nbins, (i+spread)/nbins),
                                           N=dn_n, thresh=1, interval=c(2, 30)),
                      error=function(e){})
      if(length(res)>0){
        sliding_p_non <- rbind(sliding_p_non, res)
      }
    }
    mod_p_non <- scam::scam(y ~ s(x, bs = "mpd", k=4), data = data.frame(y=sliding_p_non[,1], x=qnorm(sliding_p_non[,4])))
    pred_p_non <- predict(mod_p_non, newdata=data.frame(x=qnorm(data$LOEUFu[!is.na(data$LOEUFu)])))
    
    mod_g_non <- scam::scam(y ~ s(x, bs = "mpi", k=4), data = data.frame(y=sliding_p_non[,8], x=qnorm(sliding_p_non[,4])))
    pred_g_non <- predict(mod_g_non, newdata=data.frame(x=qnorm(data$LOEUFu[!is.na(data$LOEUFu)])))
    
    res = list(allinfo=sliding_p_non, risk_proportion = pred_p_non, gamma = pred_g_non)
    saveRDS(res, paste0('./prior_',info,'-',dn_name,'_mut-',mut_name,'.rds'))
    
    if(plot){
      pred_p_non = res[['risk_proportion']]
      pred_g_non = res[['gamma']]
      sliding_p_non = res[['allinfo']]
      
      par(mfrow=c(1,2), mar=rep(4,4))
      plot(y=pred_g_non, x=data$LOEUF[!is.na(data$LOEUF)], ylim=c(0, 35), ylab="Gamma", xlab="LOEUF")
      points(y=sliding_p_non[,8], x=sliding_p_non[,4]*2, col=2)
      plot(y=pred_p_non, x=data$LOEUF[!is.na(data$LOEUF)], ylim=c(0, 1), ylab="Proportion Risk Genes", xlab="LOEUF")
      points(y=sliding_p_non[,1], x=sliding_p_non[,4]*2, col=2)
      
    }
    
  }else{
    res = readRDS(paste0('./prior_',info,'-',dn_name,'_mut-',mut_name,'.rds'))
  }

  return(res$gamma)
}

