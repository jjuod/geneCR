library(dplyr)
library(ggplot2)
library(bigsnpr)

setwd("~/Documents/results/cr/")

bed = snp_readBed("1000g/1000G_sm.bed")
bed = snp_attach(bed)
y = snp_simuPheno(bed$genotypes, h2=0.8, M=1000)

newfam = bed$fam[,1:2]
newfam$Y = y$pheno
write.table(newfam, file="simphenos/cont.pheno", quote=F, sep="\t", col.names=F, row.names=F)


library(sim1000G)
# simulate genotypes. Split in 5 blocks to make life easier for the package
for(i in 1:5){
  blocklen = 50e6
  vcf = readVCF("1000g/1000G_EUR_chr1.vcf.gz", region_start=blocklen*(i-1), region_end=blocklen*i)
  startSimulation(vcf, totalNumberOfIndividuals=3000)
  ids = generateUnrelatedIndividuals(2000)
  genotypes = retrieveGenotypes(ids)
  
  fam = data.frame(gtindex=ids, ids, pid=0, mid=0, sex=1, pheno=0)
  writePED(vcf, fam, filename=paste0("1000g/simulated",i))
}
# can check MAFs by
#plot(rowMeans(vcf$gt1+vcf$gt2)/2, colMeans(genotypes)/2)

