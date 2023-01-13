library(dplyr)
library(tidyr)
library(ggplot2)
library(Rcapture)
setwd("~/Documents/gitrep/geneCR/")

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

# t4
allts[[11]]$p.adj = p.adjust(allts[[11]]$p.value, method="fdr", 690)
min(allts[[11]]$p.adj)

# t5
# already at FDR<5%
allts[[12]]
allts[[12]]$p.adj = allts[[12]]$p.value

# t6
# already at FDR<5%
allts[[13]]
allts[[13]]$p.adj = 0.05

# t7
allts[[14]]$p.adj = p.adjust(allts[[14]]$p.value, method="fdr", 244)
min(allts[[14]]$p.adj)

# t8
# actually unclear how the p-values were generated,
# but none would pass, assuming unadjusted
allts[[15]]$p.adj = p.adjust(allts[[15]]$p.value, method="fdr", 1393)
min(allts[[15]]$p.adj)

# t9
# (note that this includes hits from several comparisons)
# the q-values approx match an FDR adjustment for 1850
allts[[16]]$p.adj = allts[[16]]$q.value
min(allts[[16]]$p.adj)

# hit counts for the table
lapply(allts, function(x) filter(x, p.adj<=0.05) %>% summarize(n=n(), fname=max(fname))) %>%
  bind_rows %>% arrange(fname)
