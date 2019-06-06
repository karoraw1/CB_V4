
library(ape)
library(picante)
library(MicEco)
library(doSNOW)


#tree_90 = read.tree(file = "../otu_data/clustered_sequences/cluster_members.100.phylip32.mod.nh")
#tsv.data <- read.delim("../otu_data/clustered_sequences/abundances.100.tsv", row.names=1)
#prunedphy <- prune.sample(tsv.data, tree_90)
#phydist <- cophenetic(prunedphy)

tsv.data <- read.delim("../otu_data/clustered_sequences/pplacer_abundances.tsv", row.names=1)
phydf <- read.delim("../otu_data/clustered_sequences/fixed_pplacer_distmat.tsv", row.names=1)
strata_df <- read.delim("../otu_data/clustered_sequences/strata.tsv", row.names=1)
phydist = as.matrix(phydf)
strat_vec <- unname(unlist(strata_df[,'CollectionAgency']))

start_time <- Sys.time()
MNTD_w <- comdistnt.par(tsv.data, phydist, abundance.weighted = TRUE, cores=24)
save(MNTD_w, file="../otu_data/clustered_sequences/mean_nearest_taxon_distance.RData")
end_time <- Sys.time()
end_time - start_time

start_time <- Sys.time()
mntd_scores <- ses.comdistnt2(tsv.data, phydist, strata = strat_vec, abundance.weighted = TRUE, runs = 333, iterations=333, cores=24)
save(mntd_scores, file="../otu_data/clustered_sequences/standard_effects_mntd.RData")
end_time <- Sys.time()
end_time - start_time



