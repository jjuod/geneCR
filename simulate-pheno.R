library(bigsnpr)

setwd("~/Documents/results/cr/")

bed = snp_readBed("1000g/simulated_cm.bed")
bed = snp_attach(bed)
y = snp_simuPheno(bed$genotypes, h2=0.3, M=2000)

newfam = bed$fam[,1:2]
newfam$Y = y$pheno

write.table(newfam, file="simphenos/cont.pheno", quote=F, sep="\t", col.names=F, row.names=F)
