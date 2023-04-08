library(dplyr)
library(tidyr)
library(ggplot2)
library(Rcapture)

setwd("~/Documents/results/cr/")

# check mapping of SNPs to genes and create their conversion files

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

