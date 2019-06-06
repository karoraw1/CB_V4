
library(ape)
library(picante)
library(MicEco)
library(doSNOW)

#tree_90 = read.tree(file = "../otu_data/clustered_sequences/cluster_members.100.phylip32.mod.nh")
tree_90 = read.tree(file = "../otu_data/clustered_sequences/cluster_members.90.phylip32.mod.nh")
#tsv.data <- read.delim("../otu_data/clustered_sequences/abundances.100.tsv", row.names=1)
tsv.data <- read.delim("../otu_data/clustered_sequences/abundances.c90.tsv", row.names=1)
prunedphy <- prune.sample(tsv.data, tree_90)
phydist <- cophenetic(prunedphy)
#library(feather)
#phydf <- read_feather("../otu_data/clustered_sequences/cluster_members.100.distmat.feather")
#phydist = as.matrix(phydf)
#rownames(phydist) <- colnames(phydist)


start_time <- Sys.time()
MNTD_w <- comdistnt.par(tsv.data[1:10,], phydist, abundance.weighted = TRUE)
end_time <- Sys.time()
end_time - start_time

start_time <- Sys.time()
mntd_scores <- ses.mntd.par(tsv.data[1:10,], phydist, null.model = "taxa.labels", abundance.weighted = TRUE, runs = 10, iterations=50, cores=8)
end_time <- Sys.time()
end_time - start_time


start_time <- Sys.time()
MNTD_w <- comdistnt.par(tsv.data, phydist, abundance.weighted = TRUE, cores=8)
end_time <- Sys.time()
end_time - start_time


library(GUniFrac)
unifracs <- GUniFrac(tsv.data[1:10,], prunedphy, alpha = c(0, 0.5, 1))



# these both work fine. PCN stuff needs to be looked into more. 
# comdistnt is very slow though, it might take 100 hours optimistically....
#start_time <- Sys.time()
#MNTD_uw <- comdistnt(tsv.data[1:9,], phydist, abundance.weighted = TRUE)
#end_time <- Sys.time()
#end_time - start_time

#4.823542/sum(1:4)
#7.755663/sum(1:6)