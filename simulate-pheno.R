library(bigsnpr)

setwd("~/Documents/results/cr/")

bed = snp_readBed("1000g/simulated.bed")
bed = snp_attach(bed)
y = snp_simuPheno(bed$genotypes, h2=0.3, M=2000)

newfam = bed$fam[,1:2]
newfam$Y = y$pheno

write.table(newfam, file="simphenos/cont.pheno", quote=F, sep="\t", col.names=F, row.names=F)

K = 0.10

newfam$Y = as.numeric(newfam$Y>quantile(newfam$Y, 1-K))+1
write.table(newfam, file="simphenos/bin.pheno", quote=F, sep="\t", col.names=F, row.names=F)


# meta-analyse blindly
dflog = read.table("simphenos/res_bin.assoc.logistic",h=T)
dflin = read.table("simphenos/res_cont.qassoc", h=T)

dflin$W = 1/dflin$SE^2
dflog$W = 1/(dflog$BETA/dflog$STAT)^2  # z=beta/SE => SE=beta/z

# Estimates the same as:
# metafor::rma(yi=c(dflin$BETA[1], dflog$BETA[1]), sei=c(dflin$SE[1], dflog$BETA[1]/dflog$STAT[1]), method="FE")
metab = (dflin$BETA*dflin$W + dflog$BETA*dflog$W)/(dflin$W+dflog$W)
metase = sqrt(1/(dflin$W + dflog$W))
metaz = metab/metase

meta = data.frame(SNP=dflin$SNP, A1="A1", A2="A2", N=4000, Z=metaz)

write.table(meta, file="simphenos/res_meta.ldready", quote=F, row.names=F)
