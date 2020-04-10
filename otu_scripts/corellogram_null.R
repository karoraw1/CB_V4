library(vegan)

env_data_file = "../otu_data/WaterQualityData/matched_cleaned_data/all_mdata_with_habitat.txt"
select_samples <- rownames(read.delim(env_data_file, row.names=1))

# rarefied otu table
tsv.data <- read.delim("../otu_data/dispersal_selection_data/final_rarefied_table.tsv", row.names=1)

# subset 
tsv.data = tsv.data[select_samples, ]

# Remove zeros columns & rows
tsv.data =  tsv.data[rowSums(tsv.data) > 0, colSums(tsv.data) > 0]

# null model community
sim.tsv <- permatswap(tsv.data, "quasiswap", times=10)$perm[[5]]

# measure distances
exp.eco.dist = vegdist(t(tsv.data), method='bray')
exp.comm.dist = as.matrix(exp.eco.dist)
sim.eco.dist = vegdist(t(sim.tsv), method="bray")
sim.comm.dist = as.matrix(sim.eco.dist)

# read in phylo dist mat
phydf <- read.delim("../otu_data/dispersal_selection_data/not_full_tree_distances.tsv", row.names=1)

# subset it to observed abundances and convert to matrix
phydf = phydf[ colnames(tsv.data), colnames(tsv.data) ]
phydist = as.matrix(phydf)

our_breaks = c(0, 0.001, 0.005, 0.015, 0.025, 0.03, 0.04, 0.05, 0.06, 0.07, 0.09, 0.11, 0.13, 0.15, 0.17,
               0.20, 0.22, 0.24, 0.26, 0.28, 0.30, 0.4, 0.5, 0.6, 0.7, 0.9, 1.1, 2.55, 4)

sim_gram_f1 = "../otu_data/dispersal_selection_data/sim_correlog2.RData"
sim.correlog = mantel.correlog(sim.comm.dist, D.geo=phydist, mult="BH", 
                                r.type="spearman", cutoff = F,
                                nperm=999, break.pts=our_breaks)
save(sim.correlog, file=sim_gram_f1)
write(end_time - start_time, stdout())


