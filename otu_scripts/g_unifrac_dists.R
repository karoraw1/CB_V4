
library(ape)
library(picante)
tree_90 = read.tree(file = "/Volumes/KeithSSD/CB_V4/otu_data/clustered_sequences/cluster_members.100.phylip32.mod.nh")
tsv.data <- read.delim("/Volumes/KeithSSD/CB_V4/otu_data/clustered_sequences/abundances.c90.tsv", row.names=1)
prunedphy <- prune.sample(tsv.data, tree_90)
phydist <- cophenetic(prunedphy)
MNTD_w <- comdistnt(tsv.data[1:10,], phydist, abundance.weighted = TRUE)
 <- unifrac(, prunedphy)

library(GUniFrac)
unifracs <- GUniFrac(tsv.data[1:10,], prunedphy, alpha = c(0, 0.5, 1))



# these both work fine. PCN stuff needs to be looked into more. 
# comdistnt is very slow though, it might take 100 hours optimistically....
#start_time <- Sys.time()
#MNTD_uw <- comdistnt(tsv.data[1:9,], phydist, abundance.weighted = TRUE)
#end_time <- Sys.time()
#nd_time - start_time

#4.823542/sum(1:4)
#7.755663/sum(1:6)