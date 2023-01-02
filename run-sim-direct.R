options(stringsAsFactors = F)
library(dplyr)
library(tidyr)
library(ggplot2)
library(Rcapture)

setwd("~/Documents/results/cr/")


# loads gff2 and snpstogenes
load("1000g/genes_to_snps.RData")


# runs the C-R, formats the results, does some error handling
runCR = function(obsgenes, detp){
  # Processing for C-R:
  listall = unique(unlist(obsgenes))
  capt.hist = sapply(obsgenes, function(x) as.numeric(listall %in% x))
  
  if(length(listall)<2){
    # mark NAs if only one gene reported
    mod = data.frame(abundance=NA, stderr=NA, deviance=NA, df=NA, AIC=NA, BIC=NA, infoFit=NA)
  } else {
    mod = data.frame(closedp.0(capt.hist)$results)
    # mark NAs in no-overlap cases:
    if(mod$abundance[1]>1e7){
      mod$abundance = NA
      mod$stderr = NA
    }
  }
  
  mod$detp = detp
  return(mod)
}

NITER = 300

for(NCAUSALGENES in c(30, 100, 300)){
  crres = tibble()
  for(detp in c(0.1, 0.3, 0.5)){
    for(i in 1:NITER){
      print(i)
      causal = sample_n(gff2[,c("gene", "snps")], NCAUSALGENES)
      causalgenes = unique(causal$gene)
      obsgenes = list()
      
      for(run in c(1, 2)){
        foundgenes = sample(causal$gene, round(detp*NCAUSALGENES))
        foundgenes = unique(foundgenes)
        obsgenes[[length(obsgenes)+1]] = foundgenes
      }
      
      # C-R:
      mod = runCR(obsgenes, detp)
      crres = bind_rows(crres, mod)
    }
  }
  write.table(crres, paste0("crres_direct_",NCAUSALGENES,"g.csv"),
              sep="\t", quote=F, col.names=T, row.names=F)
}
