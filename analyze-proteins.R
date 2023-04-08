library(dplyr)
library(tidyr)
library(ggplot2)
library(Rcapture)
library(ggvenn)
library(knitr)
setwd("~/Documents/gitrep/geneCR/")

## --------- GENOMICS --------------

# 1. only Swedish meta
genes = list(ga = c("WNT4","HIVEP3","FAF1","TET3",
                    "LSM3","ADCY5","EEFSEC","MRPS22",
                    "ZBTB38","KCNAB1","LEKR","KDR",
                    "HAND2","EBF1","HLA-DQA1","GDAP1",
                    "FBXO32","COL27A1","TFAP4","MYOCD",
                    "TCEA2","AGTR2","RAP2C"),
             ptd = c("WNT4", "EEFSEC", "KCNAB1", "HAND2",
                     "EBF1", "HLA-DQA1", "LRP5"))
listall = unique(unlist(genes))
capt.hist = sapply(genes, function(x) as.numeric(listall %in% x))
mod = closedp(capt.hist)
mod

# 2. only Finnish meta
genes.fin = list(ga = c("WNT4","WNT3A","TET3", "ADCY5","EEFSEC",
                    "ZBTB38","KCNAB1","HAND2","KCNN2",
                    "EBF1","RHAG","COBL",
                    "GNAQ","COL27A1","AGTR2"),
             ptd = c("EEFSEC", "GC", "EBF1", "LINC02824"))
listall = unique(unlist(genes.fin))
capt.hist = sapply(genes.fin, function(x) as.numeric(listall %in% x))
mod.fin = closedp(capt.hist)
mod.fin

# 3. both metas
genes.both = c(genes, genes.fin)
listall = unique(unlist(genes.both))
capt.hist = sapply(genes.both, function(x) as.numeric(listall %in% x))
mod.both = closedp(capt.hist)
mod.both

# Make Table 1:
tidy_cr = function(mod){
  res = data.frame(mod$results)
  res$mod = rownames(res)
  rownames(res) = NULL
  res = res[,c("mod", "abundance", "stderr","AIC","infoFit")]
}
table1 = bind_rows("swed"=tidy_cr(mod),
          "fin"=tidy_cr(mod.fin),
          "both"=tidy_cr(mod.both), .id="study")
table1 %>%
  format(digits=3)
# drop poorly fitted models
table1 %>%
  filter(mod %in% c("M0", "Mt", "Mh Poisson2", "Mth Poisson2")) %>%
  select(!one_of("infoFit")) %>%
  mutate(stderr=round(stderr, 1)) %>%
  format(digits=3) %>%
  kable(format="latex")

# numbers for text
sapply(genes, length)
sapply(genes.fin, length)
table(unique(unlist(genes)) %in% unlist(genes.fin))
table(unique(unlist(genes.fin)) %in% unlist(genes))


## --------- PROTEOMICS --------------

# read in, convert to q-vals, cutoff

tfiles = list.files("./prot_lists/", pattern="t*csv")
allts = list()
for(tfile in tfiles){
  print(tfile)
  t = read.table(paste0("./prot_lists/", tfile), h=T, sep="\t", quote="")
  t$fname = substr(tfile, 1, nchar(tfile)-4)
  allts[[length(allts)+1]] = t
}

# t1
allts[[1]]$p.adj = p.adjust(allts[[1]]$p.value, method="fdr", n=460)
min(allts[[1]]$p.adj)

# t10
allts[[2]]$p.adj = p.adjust(allts[[2]]$p.value, method="fdr", n=780)
min(allts[[2]]$p.adj)

# t11a
allts[[3]]$p.adj = p.adjust(allts[[3]]$p.value, method="fdr", n=501)
min(allts[[3]]$p.adj)

# t11b
allts[[4]]$p.adj = p.adjust(allts[[4]]$p.value, method="fdr", n=501)
min(allts[[4]]$p.adj)

# t12
numpvals = allts[[5]]$p.value
numpvals = ifelse(numpvals=="< 0.00001", 0.000005, numpvals)
numpvals = as.numeric(numpvals)
allts[[5]]$p.adj = p.adjust(numpvals, method="fdr", n=624)
min(allts[[5]]$p.adj)

# t13
allts[[6]]$p.adj = p.adjust(allts[[6]]$p.value, method="fdr", n=777)
min(allts[[6]]$p.adj)

# t14
allts[[7]]$p.adj = p.adjust(allts[[7]]$p.value, method="fdr", n=325)
min(allts[[7]]$p.adj)

# t15
# (pre-adjusted by local FDR)
allts[[8]]$p.adj = allts[[8]]$p.value
max(allts[[8]]$p.adj)

# t3a
# not clear how many to adjust for, but seems like all would be
# lost anyway, even assuming ~200 spots:
allts[[9]]$p.adj = p.adjust(allts[[9]]$p.value, method="fdr", 200)
min(allts[[9]]$p.adj)

# t3b
# same
allts[[10]]$p.adj = p.adjust(allts[[10]]$p.value, method="fdr", 200)
min(allts[[10]]$p.adj)

# t5
# already at FDR<5%
allts[[11]]
allts[[11]]$p.adj = allts[[11]]$p.value

# t6
# already at FDR<5%
allts[[12]]
allts[[12]]$p.adj = 0.05

# t7
allts[[13]]$p.adj = p.adjust(allts[[13]]$p.value, method="fdr", 244)
min(allts[[13]]$p.adj)

# t9
# (note that this includes hits from several comparisons)
# the q-values approx match an FDR adjustment for 1850
allts[[14]]$p.adj = allts[[14]]$q.value
min(allts[[14]]$p.adj, na.rm=T)


# hit counts for the table
allts.pass = lapply(allts,
              function(x) select(x, -one_of("p.value", "logP", "q.value", 
                                            "Fold.change..group.Y.group.D."))) %>%
  bind_rows() %>%
  filter(p.adj<=0.05)

allts.pass %>% group_by(fname) %>%
  summarize(n())

# drop one gene with seemingly wrong ID
allts.pass = filter(allts.pass, !is.na(uniprot.id))

# make sure all synonymous IDs are converted to the same primary ID
allts.pass$uniprot.id = gsub("-[0-9]", "", allts.pass$uniprot.id)

# C-R

allts.wide.df = mutate(allts.pass[,c("fname", "uniprot.id")], det=1) %>%
  distinct() %>%
  spread(value="det", key="fname")

nrow(allts.wide.df) # 912 total
length(unique(filter(allts.pass, fname!="t5", fname!="t9-wids")$uniprot.id))
length(filter(allts.pass, fname!="t5", fname!="t9-wids")$uniprot.id)
# 311/385 unique/total

allts.wide = !is.na(as.matrix(allts.wide.df[,2:ncol(allts.wide.df)]))
allts.wide.df = cbind("id"=allts.wide.df$uniprot.id, data.frame(allts.wide))

mod = closedp(allts.wide)
mod

mod.no59 = closedp(allts.wide[,-c(7,10)])
mod.no59

closedpCI.0(allts.wide[,-c(7,10)])
closedpCI.t(allts.wide[,-c(7,10)], m="Mt")

# t5, t9 is newborn, should be removed
ggvenn(allts.wide.df, c("t9.wids", "t15.wids", "t14.wids"),
       show_percentage=F)
ggvenn(allts.wide.df, c("t9.wids", "t5", "t14.wids", "t15.wids"),
       show_percentage=F)


# Plot the overlaps for the main figure
vennlabels = tribble(~x,   ~y,   ~hjust, ~vjust,
          -1.0, -1.3, 1,      1,
          -0.9,  1.2, 0.5,    0,
          1.0,  1.2, 0.5,    0,
          0.7, -1.6, 0,      1)
vennlabels$text = c("Govia et al.", "Lee et al. (2021)", "Lee et al. (2020)", "Romero et al.")

rename(allts.wide.df, `Govia et al.`=t15.wids, `Lee et al. (2021)`=t7,  
       `Romero et al.`=t14.wids, `Lee et al. (2020)`=t12) %>%
  ggvenn(vennlabels$text, show_percentage=F, set_name_size=0) +
  geom_text(data = vennlabels,
            aes(x = x, y = y, label = text, hjust = hjust, vjust = vjust),
            size = 3.6)
ggsave("~/Documents/results/cr/plots/fig4-venn.png", width=3, height=3)
