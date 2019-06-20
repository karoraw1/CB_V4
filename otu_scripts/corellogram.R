library(vegan)

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

# create distance objects and equivalent matrices
exp.eco.dist = vegdist(t(tsv.data), method='bray')
exp.comm.dist = as.matrix(exp.eco.dist)
sim.eco.dist = vegdist(t(sim.tsv), method="bray")
sim.comm.dist = as.matrix(sim.eco.dist)

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

# make the competitors correlograms
library(ecodist)
start_time <- Sys.time()
exp.mgram = mgram(exp.eco.dist, as.dist(phydist), trace = TRUE)
save(exp.mgram, file="../otu_data/dispersal_selection_data/exp_mgram.RData")
end_time <- Sys.time()

start_time <- Sys.time()
sim.mgram = mgram(sim.comm.dist, as.dist(phydist), trace = TRUE)
save(sim.mgram, file="../otu_data/dispersal_selection_data/sim_mgram.RData")
end_time <- Sys.time()
