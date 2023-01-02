options(stringsAsFactors = F)
library(dplyr)
library(tidyr)
library(ggplot2)
library(Rcapture)

setwd("~/Documents/results/cr/")

# loads gff2 and snpstogenes
load("1000g/genes_to_snps.RData")
load("1000g/snps_to_genes.RData")

NCAUSALGENES = 200

# Load the genotypes
gt3000 = data.table::fread("1000g/simulated_3000R.raw", h=T)
gt3000 = t(as.matrix(gt3000))
dim(gt3000)

gt2000 = data.table::fread("1000g/simulated_2000R.raw", h=T)
gt2000 = t(as.matrix(gt2000))
dim(gt2000)

# lists of independent snps
fullbim = data.table::fread("1000g/simulated_1000.bim", h=F)

runGWAS = function(ys, ids, bfile){
  ys = data.frame("FID"=ids, "IID"=ids, "Y"=ys)
  write.table(ys, "simphenos/tmpy.pheno", sep="\t", col.names = T, row.names = F, quote=F)
  
  # run plink assoc:
  system(paste0("plink --bfile 1000g/",bfile," --pheno simphenos/tmpy.pheno --assoc --out simphenos/tmp_assoc"),
         ignore.stdout=T)
  
  # run plink clumping
  system(paste0("plink --bfile 1000g/",bfile,
                " --clump simphenos/tmp_assoc.qassoc --clump-p1 5e-08 ",
                "--clump-p2 0.00001 --clump-r2 0.10 --out simphenos/tmp_clumped"),
         ignore.stdout=T)

  # analyze the results:
  res = read.table("simphenos/tmp_assoc.qassoc", h=T)
  return(which(res$P<5e-8))
}

# ------------------------------------------
# two samples, M0 only

NITER = 300


for(sigma in seq(5, 35, by=5)){
  print(paste("------ sigma:", sigma, " -------------"))
  
  summaries = data.frame(h2s=NULL, numsignhits=NULL, run=NULL, i=NULL, sigma=NULL, indep=NULL)
  crres = tibble()
  
  for(i in 1:NITER){
    print(i)
    causal = sample_n(gff2[,c("gene", "snps")], NCAUSALGENES)
    causalgenes = unique(causal$gene)
    obsgenes = list()
    obsgenesI = list()
    
    for(run in c(1, 2)){
      causal$geneeff = rnorm(nrow(causal), 0, 5.0)
      # causal$geneeff = runif(nrow(causal), -10, 10)
      # can simulate variable snp effects, or just select one here:
      causalsnps = unnest(causal, snps) %>%
        group_by(gene) %>%
        sample_n(1)
      causalsnps$snpeff = causalsnps$geneeff

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
      ys = xbeta + rnorm(length(xbeta), 0, sigma)
      
      ressnps = runGWAS(ys, ids, bfile)

      foundgenes = snpstogenes$gene[ressnps]
      foundgenes = unique(foundgenes[!is.na(foundgenes)])
      obsgenes[[length(obsgenes)+1]] = foundgenes

      if(run==1){
        noverl = NA 
      } else {
        noverl = length(intersect(foundgenes, obsgenes[[1]]))
      }

      # store heritability etc:
      newrow = data.frame(h2s=var(xbeta)/var(ys), numsignhits=length(ressnps),
                          numsigngenes=length(foundgenes), noverl,
                          ntps=sum(foundgenes %in% causalgenes),
                          run=run, i=i, sigma=sigma, indep="all")
      summaries = bind_rows(summaries, newrow)
      
      # same but only w/ independent snps
      res = read.table("simphenos/tmp_clumped.clumped", h=T)
      ressnps = match(res$SNP, fullbim$V2)
      # cat("found cl ")
      # print(ressnps)

      foundgenes = snpstogenes$gene[ressnps]
      foundgenes = unique(foundgenes[!is.na(foundgenes)])
      obsgenesI[[length(obsgenesI)+1]] = foundgenes
      # print(table(foundgenes %in% causalgenes))
      # print(foundgenes[!foundgenes %in% causalgenes])

      if(run==1){
        noverl = NA
      } else {
        noverl = length(intersect(foundgenes, obsgenesI[[1]]))
      }

      # store heritability etc:
      newrow = data.frame(h2s=var(xbeta)/var(ys), numsignhits=length(ressnps),
                          numsigngenes=length(foundgenes), noverl,
                          ntps=sum(foundgenes %in% causalgenes),
                          run=run, i=i, sigma=sigma, indep="clumped")
      summaries = bind_rows(summaries, newrow)
    }

    
    # Processing for C-R:
    listall = unique(unlist(obsgenes))
    capt.hist = sapply(obsgenes, function(x) as.numeric(listall %in% x))
    
    
    # Run C-R + handle exceptions:  
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
    mod$indep = "all"
    mod$sigma = sigma
    crres = bind_rows(crres, mod)
    
    # same but only with independent genes:
    listall = unique(unlist(obsgenesI))
    capt.hist = sapply(obsgenesI, function(x) as.numeric(listall %in% x))
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

    mod$indep = "clumped"
    mod$sigma = sigma
    crres = bind_rows(crres, mod)

  }
  
  write.table(crres, paste0("crres_2000x2_200g_s", sigma, ".csv"), sep="\t", quote=F, col.names=T, row.names=F)
  write.table(summaries, paste0("studysummaries_2000x2_200g_s", sigma, ".csv"), sep="\t", quote=F, col.names=T, row.names=F)
}


