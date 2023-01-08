library(dplyr)
library(tidyr)
library(ggplot2)
library(Rcapture)
setwd("~/Documents/results/cr/")
setwd("/mnt/GU/results/cr")


# --------------
# sim-direct results

res = tibble()

for(ngenes in c(30, 100, 300)){
  res1 = read.table(paste0("crres_direct_", ngenes, "g.csv"), h=T)
  res1$ncausalgenes = ngenes
  res = bind_rows(res, res1)
}


nafracs = res %>%
  group_by(detp, ncausalgenes) %>%
  summarize(fracna=round(mean(is.na(abundance))*100))
res %>%
  ggplot(aes(x=factor(detp))) +
  geom_hline(data=data.frame(ncausalgenes=c(30,100,300)),
             aes(yintercept=ncausalgenes), lty="dashed") +
  geom_boxplot(aes(y=abundance), outlier.size=0.5) +
  geom_text(data=nafracs, aes(y=6, label=fracna), col="#5E8C84") + 
  labs(tag="% non-\nusable") +
  facet_wrap(~ncausalgenes, scales="free_x", labeller=label_both) +
  scale_y_log10() + xlab("detection prob.") + ylab("estimate") +
  scale_fill_brewer(palette="Set3", direction=-1, name=NULL) + theme_bw() +
  theme(plot.tag.position = c(0.9,0.3), plot.tag=element_text(size=11, colour="#5E8C84"),
        legend.text=element_text(size=11))
ggsave("plots/fig1a-simdirect.png", width=4, height=3)



# --------------
# sim1 results

res = tibble()

for(sigma in seq(10, 40, by=5)){
  res1 = read.table(paste0("crres_2000x2_s", sigma, ".csv"), h=T)
  res = bind_rows(res, res1)
}
res$ncausalgenes = 50

for(sigma in seq(5, 35, by=5)){
  res1 = read.table(paste0("crres_2000x2_200g_s", sigma, ".csv"), h=T)
  res1$ncausalgenes = 200
  res = bind_rows(res, res1)
}


# study summaries
stsums = tibble()

for(sigma in seq(10, 40, by=5)){
  stsum1 = read.table(paste0("studysummaries_2000x2_s", sigma, ".csv"), h=T)
  stsums = bind_rows(stsums, stsum1)
}
stsums$ncausalgenes = 50

for(sigma in seq(5, 35, by=5)){
  stsum1 = read.table(paste0("studysummaries_2000x2_200g_s", sigma, ".csv"), h=T)
  stsum1$ncausalgenes = 200
  stsums = bind_rows(stsums, stsum1)
}

# summaries
group_by(stsums, sigma, indep, ncausalgenes) %>%
  summarize(median(h2s), nhitsnps=median(numsignhits), numgenes=median(numsigngenes),
            noverl=median(noverl, na.rm=T), ntp=median(ntps)) %>%
  print.data.frame
# h2 range 26-99 %

# other info:
# group_by(stsums, sigma, indep) %>%
#   summarize(median(h2s), nhitsnps=median(numsignhits), numgenes=median(numsigngenes),
#             noverl=median(noverl, na.rm=T), ntp=median(ntps)) %>%
#   print.data.frame

nafracs = filter(res, indep=="clumped") %>%
  group_by(sigma, ncausalgenes, indep) %>%
  summarize(fracna=round(mean(is.na(abundance))*100))
res %>%
  ggplot(aes(x=factor(sigma), fill=indep)) +
  geom_hline(data=data.frame(ncausalgenes=c(50,200)),
             aes(yintercept=ncausalgenes), lty="dashed") +
  geom_boxplot(aes(y=abundance), outlier.size=0.5) +
  geom_text(data=nafracs, aes(y=8.5, label=fracna), col="#5E8C84") + 
  labs(tag="% non-\nusable") +
  facet_wrap(~ncausalgenes, scales="free_x", labeller=label_both) +
  scale_y_log10() + xlab("residual SD") + ylab("estimate") +
  scale_fill_brewer(palette="Set3", direction=-1, name=NULL) + theme_bw() +
  theme(plot.tag.position = c(0.87,0.2), plot.tag=element_text(size=11, colour="#5E8C84"),
        legend.text=element_text(size=11))
ggsave("plots/fig1-simres.png", width=7, height=3)


# ----------------
# sim 2

res1 = read.table(paste0("crres_x3b_g200_s50.csv"), h=T, sep="\t")
res1$ncausalgenes = 200
res2 = read.table(paste0("crres_x3b_s50.csv"), h=T, sep="\t")
res2$ncausalgenes = 50

# drop behavioral-component models
res = bind_rows(res1, res2) %>% filter(!model %in% c("Mb", "Mbh"))
ggplot(res, aes(y=abundance, x=model, group=model)) +
  geom_hline(data=data.frame(ncausalgenes=c(50,200)),
             aes(yintercept=ncausalgenes), lty="dashed") +
  geom_boxplot() + facet_wrap(~ncausalgenes) + 
  scale_y_log10() + ylab("estimate") +
  coord_cartesian(ylim=c(1e-1,1e5)) + 
  theme_bw() + theme(axis.text.x = element_text(angle=60))
ggsave("plots/fig3-sim2res.png", width=7, height=3)

stsums1 = read.table("studysummaries_x3b_g200_s50.csv", h=T)
stsums1$ncausalgenes = 200
stsums2 = read.table("studysummaries_x3b_s50.csv", h=T)
stsums2$ncausalgenes = 50

bind_rows(stsums1, stsums2) %>%
  group_by(ncausalgenes) %>%
  summarize(hits=median(numsignhits), genes=median(numsigngenes), tps=median(ntps))

##### -------------- tests -------------


load("/mnt/GU/results/cr/tmp_obsgenes_x3b_g200_s50.RData")
load("/mnt/GU/results/cr/tmp_causalgenes_x3b_g200_s50.RData")
obsgenes1 = obsgenes_all[[50]][[3]]
table(obsgenes1 %in% causalgenes_all[[1]])
gff = gff2
gff$causal = gff$gene %in% causalgenes_all[[50]]
gff$detected = gff$gene %in% obsgenes1
gff$pos = 1:nrow(gff)
filter(gff, causal | detected)

genenums = t(sapply(obsgenes_all, function(x) sapply(x, length) ))

abvec2 = sapply(obsgenes_all, function(x) runCR(x)$abundance[2])
genenums = data.frame(genenums, abvec)
genenums$abvecT = abvec2
group_by(genenums, X1>10) %>%
  summarize(median(abvec), median(abvecT))

genes3k = lapply(obsgenes_all, "[[", 3)
genes3k = mapply(list, genes3k[seq(1, 50, 2)], genes3k[seq(2,50,2)], SIMPLIFY = F)
abvec = sapply(genes3k, function(x) runCR(x)$abundance[1])

match(obsgenes_all[[1]][[3]], gff2$gene)

#### ----------- more tests --------------

gff22 = gff2
gff22$snps = lapply(gff2$snps, "[", 1)
gff22$snps = lapply(1:nrow(gff2), function(x) (5*x):(5*x+4))
snpstogenes2 = unnest(gff22, snps) %>%
  right_join(data.frame(snps=1:88060), by="snps") %>%
  arrange(snps)

NCAUSALGENES = 200

NITER = 70

DETP = 0.2
sigma = 20

crres = tibble()

for(i in 1:NITER){
  print(i)
  causal = sample_n(gff22[,c("gene", "snps")], NCAUSALGENES)
  causalgenes = unique(causal$gene)
  obsgenes = list()
  
  for(run in c(1, 2)){
    causal$geneeff = rnorm(nrow(causal), 0, 0.7)
    causalsnps = unnest(causal, snps)
    causalsnps = sample_frac(causalsnps, DETP)
    # causalsnps = sample_frac(causal, DETP) %>% unnest(snps)
    # causalsnps$snpeff = rnorm(nrow(causalsnps), causalsnps$geneeff, 0.3)
    
    # ressnps = causalsnps$snps[abs(causalsnps$snpeff)>1.3]
    # ressnps = causalsnps$snps[rbinom(nrow(causalsnps), 1, DETP)==1]
    ressnps = causalsnps$snps
    foundgenes = snpstogenes2$gene[ressnps]
    
    foundgenes = unique(foundgenes[!is.na(foundgenes)])
    print(paste("snps", length(ressnps), " -- genes", length(foundgenes)))
    obsgenes[[length(obsgenes)+1]] = foundgenes
  }

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

median(crres$abundance, na.rm=T)
summary(crres$abundance)
table(is.na(crres$abundance))

