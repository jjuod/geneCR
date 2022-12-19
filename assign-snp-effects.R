library(dplyr)
library(tidyr)
library(ggplot2)
library(Rcapture)

setwd("~/Documents/results/cr/")

bim = read.table("1000g/simulated_3000.bim")

nrow(bim)

gff = read.table("Homo_sapiens.GRCh37.genes", sep="\t", quote="")
colnames(gff) = c("CHROM", "SOURCE", "TYPE", "START", "END", "V6", "strand", "V8", "ANN")
gff = filter(gff, CHROM==1)
nrow(gff)

gff$gene = gsub("ID=gene:(.*?);.*", "\\1", gff$ANN)

# evaluate the snp density
gff$nsnp = gff$nsnp10kb = gff$nsnp100kb = NA
for(i in seq_along(gff$START)){
  gff$nsnp[i] = filter(bim, between(V4, gff[i,"START"], gff[i, "END"])) %>% nrow
  gff$nsnp10kb[i] = filter(bim, between(V4, gff[i,"START"]-1e4, gff[i, "END"]+1e4)) %>% nrow
  gff$nsnp100kb[i] = filter(bim, between(V4, gff[i,"START"]-1e5, gff[i, "END"]+1e5)) %>% nrow
}

# <4 % genes have no snps with this cutoff
mean(gff$nsnp10kb==0)

allgenes = unique(gff$gene)
length(allgenes)  # 2677
nrow(gff)

# collect the snps for each gene
gff$snps = NA
for(i in seq_along(gff$START)){
  nearbysnps = which(between(bim$V4, gff[i,"START"]-1e4, gff[i, "END"]+1e4))
  gff$snps[i] = list(nearbysnps)
}
# gff = unnest(gff, snps)
# nrow(gff)


snpstogenes = unnest(gff[,c("snps", "gene", "START", "END")], snps)
snpstogenes$snppos = bim$V4[snpstogenes$snps]
# somewhat hacky assignment of snps to closest genes, but it doesn't really matter:
snpstogenes = mutate(snpstogenes, dist = abs((START+END)/2-snppos)) %>%
  group_by(snps) %>% 
  slice_min(dist, n=1, with_ties=F) %>%
  arrange(snps)
nrow(snpstogenes)
snpstogenes = ungroup(snpstogenes) %>%
  complete(snps=1:nrow(bim))
nrow(snpstogenes)

# on average each gene's snps map to 2 genes however,
# which would inflate the polygenicity estimate
unnest(gff[,c("snps", "gene", "START", "END")], snps) %>%
  left_join(snpstogenes, by="snps", suffix=c(".orig", ".snpmap")) %>%
  group_by(gene.orig) %>%
  summarize(ntarg=n_distinct(gene.snpmap)) %>%
  summarize(mean(ntarg), median(ntarg))

# so reverse the mapping, i.e. snp is assigned to the closest gene
gff2 = filter(snpstogenes, !is.na(gene)) %>%
  group_by(gene) %>%
  summarize(snps=list(snps))

save(gff2, file="1000g/genes_to_snps.RData")
save(snpstogenes, file="1000g/snps_to_genes.RData")

NCAUSALGENES = 50

# Load the genotypes
gt3000 = data.table::fread("1000g/simulated_3000R.raw", h=T)
gt3000 = t(as.matrix(gt3000))
dim(gt3000)
# gt3000[1:10, 1:10]

gt2000 = data.table::fread("1000g/simulated_2000R.raw", h=T)
gt2000 = t(as.matrix(gt2000))
dim(gt2000)
# gt2000[1:10, 1:10]

gt1000 = data.table::fread("1000g/simulated_1000R.raw", h=T)
gt1000 = t(as.matrix(gt1000))
dim(gt1000)
# gt1000[1:10, 1:10]


# ------------------------------------------
# SIMULATIONS

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

NITER = 50

for(i in 1:NITER){
  print(i)
  causal = sample_n(gff2[,c("gene", "snps")], NCAUSALGENES)
  causalgenes = unique(causal$gene)
  obsgenes = list()
  
  for(run in c(1000, 2000, 3000)){
    causal$geneeff = rnorm(nrow(causal), 0, 0.5)
    causalsnps = unnest(causal, snps)
    causalsnps$snpeff = rnorm(nrow(causalsnps), causalsnps$geneeff, 0.5)
    
    # simulate phenos
    xbeta = NULL
    ids = 1:run
    bfile = paste0("simulated_", run)
    
    if(run==1000){
      xbeta = as.numeric(causalsnps$snpeff %*% gt1000[causalsnps$snps+6,])
    } else if(run==2000){
      xbeta = as.numeric(causalsnps$snpeff %*% gt2000[causalsnps$snps+6,])
    } else {
      xbeta = as.numeric(causalsnps$snpeff %*% gt3000[causalsnps$snps+6,])
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
    mod = data.frame(closedp(capt.hist)$results)
    # mark NAs in no-overlap cases:
    notconverged = which(mod$abundance>1e7 | mod$abundance<0)
    mod$abundance[notconverged] = NA
    mod$stderr[notconverged] = NA
  }
  mod$i = i
  mod$model = rownames(mod)
  
  crres = bind_rows(crres, mod)
}

# (the unnamed one is s=30)
write.table(crres, "crres_x3_s40.csv", sep="\t", quote=F, col.names=T, row.names=F)
write.table(summaries, "studysummaries_x3_s40.csv", sep="\t", quote=F, col.names=T, row.names=F)

crres = read.table("crres_x3_s40.csv", h=T)
summaries = read.table("studysummaries_x3_s40.csv", h=T)

summaries
crres

mean(crres$abundance, na.rm=T)
median(crres$abundance, na.rm=T)
table(is.na(crres$abundance))

mean(summaries$h2s)
median(summaries$numsignhits)

group_by(crres, model) %>%
  summarize(me = mean(abundance, na.rm=T), md= median(abundance, na.rm=T))


# ----------------
# analyze the runs

crres = read.table("crres_2000x2.csv", h=T)
summaries = read.table("studysummaries_2000x2.csv", h=T)

median(crres$abundance, na.rm=T)
mean(crres$abundance, na.rm=T)
