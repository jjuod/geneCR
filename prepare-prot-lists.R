library(dplyr)
library(tidyr)
library(ggplot2)

# Rough script to show what methods were needed for parsing the mini-review
# protein lists into the ones provided in prot_lists/.
# Basically, ensuring a matchable uniprot id for each protein.

setwd("./prot_lists/")

# read all ids
uniprotids = read.table("~/Documents/results/cr/review/uniprot-id-list.tsv",
                        h=T, sep="\t", quote = "", comment.char = "")
uniprotids = separate_rows(uniprotids, Gene.Names, sep=" ")
nrow(uniprotids)
colnames(uniprotids)[1] = "uniprot.id"
colnames(uniprotids)[3] = "protein"

uniprotids$Entry.Name = gsub("_HUMAN", "", uniprotids$Entry.Name)

# find prots w/ synonyms
altaccs = readLines("~/Documents/results/cr/review/uniprot-accessions.txt")
altaccs = strsplit(altaccs, ";", fixed=T)
altaccs = altaccs[sapply(altaccs, length)>1]

altaccs = data.frame(pri=sapply(altaccs, "[[", 1), 
           sec=sapply(altaccs, function(x) paste(x[2:length(x)], collapse=",")))
altaccs = separate_rows(altaccs, sec, sep=",")
altaccs$sec = trimws(altaccs$sec)

# get a unique id list for converting names
uids.uniq = distinct(uniprotids[,c("uniprot.id", "protein")])
# pick an arbitrary id, if for some reason there are multiple
# protein accessions w/ identical names in uniprot:
uids.uniq = group_by(uids.uniq, protein) %>% top_n(wt=desc(uniprot.id), n=1)

# same for converting gene names
gids.uniq = distinct(uniprotids[,c("uniprot.id", "protein", "Gene.Names")])
gids.uniq = group_by(gids.uniq, Gene.Names) %>% top_n(wt=desc(uniprot.id), n=1)

# same for converting "entry ids" ~= gene names
gids.uniq = distinct(uniprotids[,c("uniprot.id", "protein", "Entry.Name")])


# read data
allcsvs = list.files(pattern="*.csv")

dfs = list()
for(f in allcsvs){
  print(f)
  df = read.table(f, h=T, sep="\t", quote="")
  print(head(df))
  dfs[[length(dfs)+1]] = df
}

# convert ids to primary ones
tmpdf = dfs[[18]]

head(tmpdf)

table(tmpdf$uniprot.id %in% uniprotids$uniprot.id)
anti_join(tmpdf, uniprotids, by=c("uniprot.id"))

table(tmpdf$protein %in% uniprotids$protein)
anti_join(tmpdf, uniprotids, by="protein")

table(tmpdf$gene %in% uniprotids$Gene.Names)

table(tmpdf$gene %in% uniprotids$Entry.Name)


# export w/ added ids if needed
# tmpdf = left_join(tmpdf, uids.uniq, by="protein")[,c("uniprot.id", "protein", "p.value")]
# tmpdf = left_join(tmpdf, gids.uniq, by=c("gene"="Gene.Names"),
#                   suffix=c(".orig", ""))[,c("uniprot.id", "protein", "gene", "protein.orig", "p.value")]
tmpdf$uniprot.id = sapply(strsplit(tmpdf$Protein.IDs, ";"), "[[", 1)

tmpdf = tmpdf[,c("uniprot.id", "protein.names", "Protein.IDs", "Gene.names.SINGLE", "logP", "q.value")]

tmpdf = left_join(tmpdf, gids.uniq, by=c("gene"="Entry.Name"),
                  suffix=c(".orig", ""))[,c("uniprot.id", "protein", "gene", "protein.orig", "p.value")]
tmpdf = left_join(tmpdf, gids.uniq, by=c("gene"="Entry.Name"),
                  suffix=c(".orig", ""))[,c("uniprot.id", "protein", "gene", "p.value")]
tmpdf = left_join(tmpdf, gids.uniq, by=c("gene"="Gene.Names"),
                  suffix=c(".orig", ""))[,c("uniprot.id", "protein", "gene")]
tmpdf = left_join(dfs[[16]], dfs[[17]], by="Spot.label")[,c("uniprot.id", "protein", "gene", "p.value")]

write.table(tmpdf, "~/Documents/gitrep/geneCR/prot_lists/t9-wids.csv", row.names=F, col.names=T, quote=F, sep="\t")


# remove isoforms
tmpdf$uniprot.id = gsub("-[0-9]", "", tmpdf$uniprot.id)

anti_join(tmpdf, uniprotids, by=c("uniprot.id"))

tmpdf = left_join(tmpdf, altaccs, by=c("uniprot.id"="sec"))
tmpdf$pri = ifelse(is.na(tmpdf$pri), tmpdf$uniprot.id, tmpdf$pri)

table(tmpdf$pri %in% uniprotids$uniprot.id)  # nice, sorted

anti_join(tmpdf, uniprotids, by=c("pri"="uniprot.id"))


