library(vegan)
library(ape)
library(picante)
library(MicEco)
library(doSNOW)


# Read in data
tsv.data <- read.delim("../otu_data/dispersal_selection_data/final_rarefied_table.tsv", row.names=1)

# Remove zeros columns
tsv.data =  tsv.data[, colSums(tsv.data != 0) > 0]

# read in phylo dist mat
phydf <- read.delim("../otu_data/dispersal_selection_data/not_full_tree_distances.tsv", row.names=1)

# subset it to observed abundances and convert to matrix
phydf = phydf[ colnames(tsv.data), colnames(tsv.data) ]
phydist = as.matrix(phydf)

# create simulated community
sim.tsv <- permatswap(tsv.data, "quasiswap", times=10)$perm[[5]]

# create distance matrices
exp.comm.dist = as.matrix(vegdist(t(tsv.data), method='bray'))
sim.comm.dist = as.matrix(vegdist(t(sim.tsv), method="bray"))

# create real correlogram 
start_time <- Sys.time()
exp.correlog = mantel.correlog(exp.comm.dist, D.geo=phydist, mult="BH", nperm=999)
write(exp.correlog, stdout())
save(exp.correlog, file="../otu_data/dispersal_selection_data/exp_correlog.RData")
end_time <- Sys.time()

# create simulated correlogram
start_time <- Sys.time()
sim.correlog = mantel.correlog(sim.comm.dist, D.geo=phydist, mult="BH", nperm=999)
write(sim.correlog, stdout())
save(sim.correlog, file="../otu_data/dispersal_selection_data/sim_correlog.RData")
end_time <- Sys.time()


### Calculate NTR and NTI
## Test Data
start_time <- Sys.time()
MNTD_w <- comdistnt.par(tsv.data, phydist, abundance.weighted = TRUE, cores=1)
write(MNTD_w, stdout())
save(MNTD_w, file="../otu_data/clustered_sequences/test_mean_nearest_taxon_distance.RData")
mntd_scores <- ses.comdistnt2(tsv.data, phydist, method = "quasiswap", strata = NULL, abundance.weighted = TRUE, runs = 5, cores=1)
write(mntd_scores$comdistnt.obs.z, stdout())
save(mntd_scores, file="../otu_data/clustered_sequences/test_standard_effects_mntd.RData")
end_time <- Sys.time()
end_time - start_time

## Real Data
start_time <- Sys.time()
tsv.data <- read.delim("../otu_data/clustered_sequences/pplacer_abundances.tsv", row.names=1)
phydf <- read.delim("../otu_data/clustered_sequences/fixed_pplacer_distmat.tsv", row.names=1)
phydist = as.matrix(phydf)
MNTD_w <- comdistnt.par(tsv.data, phydist, abundance.weighted = TRUE, cores=24)
save(MNTD_w, file="../otu_data/clustered_sequences/mean_nearest_taxon_distance.RData")
end_time <- Sys.time()
end_time - start_time
mntd_scores <- ses.comdistnt2(tsv.data, phydist, method = "quasiswap", strata = NULL, abundance.weighted = TRUE, runs = 500, cores=24)
save(mntd_scores, file="../otu_data/clustered_sequences/standard_effects_mntd.RData")
end_time <- Sys.time()
end_time - start_time


