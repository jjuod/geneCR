library(dplyr)
library(ggplot2)
library(sim1000G)

setwd("~/Documents/results/cr/")

nind = commandArgs(T)
nind = as.numeric(nind[1])

# simulate genotypes. Split in blocks to make life easier for the package
# (14th block is empty, centromere)
for(i in c(1:13,15:25)){
  blocklen = 10e6
  gc()
  vcf = readVCF("1000g/1000G_EUR_chr1.vcf.gz", region_start=blocklen*(i-1), region_end=blocklen*i,
                maxNumberOfVariants=1e5)
  startSimulation(vcf, totalNumberOfIndividuals=nind+1000)
  ids = generateUnrelatedIndividuals(nind)
  genotypes = retrieveGenotypes(ids)
  
  fam = data.frame(gtindex=ids, ids, pid=0, mid=0, sex=1, pheno=0)
  writePED(vcf, fam, filename=paste0("1000g/tmp/simulated",i))
}
# can check MAFs by
#plot(rowMeans(vcf$gt1+vcf$gt2)/2, colMeans(genotypes)/2)