library(dplyr)
library(tidyr)
library(ggplot2)
library(Rcapture)

setwd("~/Documents/results/cr/")

# loads gff2 and snpstogenes
load("1000g/genes_to_snps.RData")
load("1000g/snps_to_genes.RData")

# lists of independent snps
fullbim = data.table::fread("1000g/simulated_1000.bim", h=F)

NCAUSALGENES = 200

# Load the genotypes
gt9000 = data.table::fread("1000g/simulated_9000R.raw", h=T)
gt9000 = t(as.matrix(gt9000))
dim(gt9000)

gt3000 = data.table::fread("1000g/simulated_3000R.raw", h=T)
gt3000 = t(as.matrix(gt3000))
dim(gt3000)

gt2000 = data.table::fread("1000g/simulated_2000R.raw", h=T)
gt2000 = t(as.matrix(gt2000))
dim(gt2000)

gt1000 = data.table::fread("1000g/simulated_1000R.raw", h=T)
gt1000 = t(as.matrix(gt1000))
dim(gt1000)


runGWAS = function(ys, ids, bfile){
  ys = data.frame("FID"=ids, "IID"=ids, "Y"=ys)
  write.table(ys, "simphenos/tmpy.pheno", sep="\t", col.names = T, row.names = F, quote=F)
  
  # run plink assoc:
  system(paste0("plink --bfile 1000g/",bfile," --pheno simphenos/tmpy.pheno --assoc --out simphenos/tmp_assoc"),
         ignore.stdout=T)
  
  # run plink clumping
  system(paste0("echo 'SNP\nNA' > simphenos/tmp_clumped.clumped; ",
                "plink --bfile 1000g/",bfile,
                " --clump simphenos/tmp_assoc.qassoc --clump-p1 5e-08 ",
                "--clump-p2 0.00001 --clump-r2 0.10 --out simphenos/tmp_clumped"),
         ignore.stdout=T)
  
  # analyze the results:
  # res = read.table("simphenos/tmp_assoc.qassoc", h=T)
  res = read.table("simphenos/tmp_clumped.clumped", h=T)$SNP
  # return(which(res$P<5e-8))
  return(res)
}

# runs the C-R, formats the results, does some error handling
runCR = function(obsgenes){
  # Processing for C-R:
  listall = unique(unlist(obsgenes))
  capt.hist = sapply(obsgenes, function(x) as.numeric(listall %in% x))
  
  if(length(listall)<2 | dim(capt.hist)[2]<2){
    # mark NAs if only one gene reported
    # or if only one study had detections
    mod = data.frame(abundance=NA, stderr=NA, deviance=NA, df=NA,
                     AIC=NA, BIC=NA, infoFit=NA)
  } else {
    mod = closedp(capt.hist)
    mod = data.frame(mod$results)
    # mark NAs in no-overlap cases:
    notconverged = which(mod$abundance>1e7 | mod$abundance<0)
    mod$abundance[notconverged] = NA
    mod$stderr[notconverged] = NA
  }
  mod$model = rownames(mod)
  return(mod)
}

# ------------------------------------------
# three samples, with list heterogeneity, M0 and others tested

summaries = data.frame(h2s=NULL, numsignhits=NULL, run=NULL, i=NULL)
crres = tibble()

NITER = 300

for(i in 1:NITER){
  print(i)
  causal = sample_n(gff2[,c("gene", "snps")], NCAUSALGENES)
  causalgenes = unique(causal$gene)
  obsgenes = list()
  
  for(run in c(1000, 3000, 9000)){
    # causal$geneeff = runif(nrow(causal), -10, 10)
    # causalsnps = unnest(causal, snps)
    # causalsnps$snpeff = rnorm(nrow(causalsnps), causalsnps$geneeff, 0.5)
    
    causal$geneeff = rnorm(nrow(causal), 0, 5.0)
    causalsnps = unnest(causal, snps) %>%
      group_by(gene) %>%
      sample_n(1)
    causalsnps$snpeff = causalsnps$geneeff
    
    # simulate phenos
    xbeta = NULL
    ids = 1:run
    bfile = paste0("simulated_", run)
    
    if(run==1000){
      xbeta = as.numeric(causalsnps$snpeff %*% gt1000[causalsnps$snps+6,])
    } else if(run==2000){
      xbeta = as.numeric(causalsnps$snpeff %*% gt2000[causalsnps$snps+6,])
    } else if(run==3000){
      xbeta = as.numeric(causalsnps$snpeff %*% gt3000[causalsnps$snps+6,])
    } else {
      xbeta = as.numeric(causalsnps$snpeff %*% gt9000[causalsnps$snps+6,])
    }
    ys = xbeta + rnorm(length(xbeta), 0, 50)
    
    res = runGWAS(ys, ids, bfile)
    ressnps = match(res, fullbim$V2)
    
    foundgenes = snpstogenes$gene[ressnps]
    foundgenes = unique(foundgenes[!is.na(foundgenes)])
    obsgenes[[length(obsgenes)+1]] = foundgenes

    # store heritability etc:
    newrow = data.frame(h2s=var(xbeta)/var(ys), numsignhits=length(ressnps),
                        numsigngenes=length(foundgenes),
                        ntps=sum(foundgenes %in% causalgenes),
                        run=run, i=i)
    summaries = bind_rows(summaries, newrow)
  }

  # C-R:
  mod = runCR(obsgenes)
  crres = bind_rows(crres, mod)
}

write.table(crres, "crres_x3b_g200_s50.csv", sep="\t", quote=F, col.names=T, row.names=F)
write.table(summaries, "studysummaries_x3b_g200_s50.csv", sep="\t", quote=F, col.names=T, row.names=F)
