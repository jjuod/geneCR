library(dplyr)
library(tidyr)
library(ggplot2)

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

NCAUSALGENES = 50

# Load the genotypes
gt3000 = data.table::fread("1000g/simulated_3000R.raw", h=T)
gt3000 = t(as.matrix(gt3000))
dim(gt3000)
# gt3000[1:10, 1:10]

gt1000 = data.table::fread("1000g/simulated_1000R.raw", h=T)
gt1000 = t(as.matrix(gt1000))
dim(gt1000)
# gt1000[1:10, 1:10]



# ------------------------------------------
# simulation rounds

# TODO store these:
h2s = c()
causalgenes = list()
numsignhits = c()
obsgenes = list()


causal = sample_n(gff[,c("gene", "snps")], NCAUSALGENES)
causal$geneeff = rnorm(nrow(causal), 0, 0.5)
causal = unnest(causal, snps)
causal$snpeff = rnorm(nrow(causal), causal$geneeff, 0.5)

# simulate phenos
xbeta = as.numeric(causal$snpeff %*% gt3000[causal$snps+6,])
ys = xbeta + rnorm(length(ys), 0, 50)

# heritability:
var(xbeta)/var(ys)

ys = data.frame("FID"=1:length(ys), "IID"=1:length(ys), "Y"=ys)
write.table(ys, "simphenos/tmpy.pheno", sep="\t", col.names = T, row.names = F, quote=F)

# run plink assoc:
system("plink --bfile 1000g/simulated_3000 --pheno simphenos/tmpy.pheno --assoc --out simphenos/tmp_assoc")

# analyze the results:
res = read.table("simphenos/tmp_assoc.qassoc", h=T)
nrow(causal)
table(res$P<5e-8)

head(causal)
which(res$P<5e-8)

unique(causal$gene)
foundgenes = snpstogenes$gene[which(res$P<5e-8)]
unique(foundgenes[!is.na(foundgenes)])

# TODO try out & set up Rcapture analysis part
list1 = unique(causal$gene)
list2 = unique(foundgenes[!is.na(foundgenes)])

listall = union(list1, list2)
sapply(list(list1, list2), function(x) listall %in% x)
