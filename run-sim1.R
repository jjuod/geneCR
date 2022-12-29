library(dplyr)
library(tidyr)
library(ggplot2)
library(Rcapture)

setwd("~/Documents/results/cr/")

# loads gff2 and snpstogenes
load("1000g/genes_to_snps.RData")
load("1000g/snps_to_genes.RData")

NCAUSALGENES = 50

# Load the genotypes
gt3000 = data.table::fread("1000g/simulated_3000R.raw", h=T)
gt3000 = t(as.matrix(gt3000))
dim(gt3000)

gt2000 = data.table::fread("1000g/simulated_2000R.raw", h=T)
gt2000 = t(as.matrix(gt2000))
dim(gt2000)


runGWAS = function(ys, ids, bfile){
  ys = data.frame("FID"=ids, "IID"=ids, "Y"=ys)
  write.table(ys, "simphenos/tmpy.pheno", sep="\t", col.names = T, row.names = F, quote=F)
  
  # run plink assoc:
  system(paste0("plink --bfile 1000g/",bfile," --pheno simphenos/tmpy.pheno --assoc --out simphenos/tmp_assoc"),
         ignore.stdout=T)
  
  # analyze the results:
  res = read.table("simphenos/tmp_assoc.qassoc", h=T)
  return(which(res$P<5e-8))
}

# ------------------------------------------
# two samples, M0 only

summaries = data.frame(h2s=NULL, numsignhits=NULL, run=NULL, i=NULL)
crres = tibble()

NITER = 500

for(sigma in c(30, 40)){
  print(paste("------ sigma:", sigma, " -------------"))
  for(i in 1:NITER){
    print(i)
    causal = sample_n(gff2[,c("gene", "snps")], NCAUSALGENES)
    causalgenes = unique(causal$gene)
    obsgenes = list()
    
    for(run in c(1, 2)){
      causal$geneeff = rnorm(nrow(causal), 0, 0.5)
      causalsnps = unnest(causal, snps)
      causalsnps$snpeff = rnorm(nrow(causalsnps), causalsnps$geneeff, 0.5)
      
      # simulate phenos
      xbeta = NULL
      ids = NULL
      bfile = NULL
      if(run==1){
        xbeta = as.numeric(causalsnps$snpeff %*% gt3000[causalsnps$snps+6,1:2000])
        ids = 1:2000
        bfile = "simulated_3000"
      } else {
        xbeta = as.numeric(causalsnps$snpeff %*% gt2000[causalsnps$snps+6,1:2000])
        ids = 1:2000
        bfile = "simulated_2000"
      }
      ys = xbeta + rnorm(length(xbeta), 0, 40)
      
      ressnps = runGWAS(ys, ids, bfile)
      foundgenes = snpstogenes$gene[ressnps]
      foundgenes = unique(foundgenes[!is.na(foundgenes)])
      obsgenes[[length(obsgenes)+1]] = foundgenes
      
      # store heritability etc:
      newrow = data.frame(h2s=var(xbeta)/var(ys), numsignhits=length(ressnps),
                          numsigngenes=length(foundgenes), run=run, i=i)
      summaries = bind_rows(summaries, newrow)
    }
    # print(summaries)
    
    # Processing for C-R:
    listall = unique(unlist(obsgenes))
    capt.hist = sapply(obsgenes, function(x) as.numeric(listall %in% x))
    # print(capt.hist)
    
    # Run C-R + handle exceptions:  
    if(length(listall)==1){
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
    
    crres = bind_rows(crres, mod)
  }
  
  write.table(crres, paste0("crres_2000x2_s", sigma, ".csv"), sep="\t", quote=F, col.names=T, row.names=F)
  write.table(summaries, paste0("studysummaries_2000x2_s", sigma, ".csv"), sep="\t", quote=F, col.names=T, row.names=F)
}

crres = read.table("crres_2000x2_s40.csv", h=T)
summaries = read.table("studysummaries_2000x2_s40.csv", h=T)

summaries
crres

mean(crres$abundance, na.rm=T)
median(crres$abundance, na.rm=T)
table(is.na(crres$abundance))

mean(summaries$h2s)
median(summaries$numsignhits)


