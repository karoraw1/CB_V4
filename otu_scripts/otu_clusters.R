library(pvclust)
library(parallel)


env_data_file = "/Volumes/KeithSSD/CB_V4/otu_data/WaterQualityData/matched_cleaned_data/all_mdata_with_habitat.txt"
env_data = read.delim(env_data_file, row.names=1)
env_data = env_data[env_data$CollectionAgency != 'Preheim', ]
select_samples <- rownames(env_data)

tsv.data <- read.delim("/Volumes/KeithSSD/CB_V4/otu_data/dispersal_selection_data/final_rarefied_table.tsv", row.names=1)
tsv.data = tsv.data[select_samples, ]
tsv.data2 = (tsv.data + 0.1)/rowSums(tsv.data + 0.1)
tsv.data3 = t(t(tsv.data2)/rowSums(t(tsv.data2)))

cl_inst <- makeCluster(3, type = "PSOCK")
res.pv <- pvclust(tsv.data3, method.hclust = "ward.D2", method.dist = "euclidean", parallel=cl_inst, iseed = 1)
stopCluster(cl_inst)

# save the data
save(res.pv, "/Volumes/KeithSSD/CB_V4/otu_data/otuclusters.RData")

# view cluster sizes 
cluster_sizes = data.frame('sizes'=c(), 'cutoff'=c())

counter=0
for (cutoff_ in seq(0.94, 0.99, 0.005)){
	clusters_i <- pvpick(res.pv, alpha=cutoff_)$clusters
	clust_lens = unlist(lapply(clusters_i, FUN=length))
	cat(cutoff_, length(clust_lens), "\n")
	counter=counter+length(clust_lens)
	cluster_sizes[(counter+1):(counter+length(clust_lens)), 'sizes'] = clust_lens
	cluster_sizes[(counter+1):(counter+length(clust_lens)), 'cutoff'] = cutoff_
	counter=counter+length(clust_lens)
}

# see how many clusters by each cutoff value
table(cluster_sizes$cutoff)

library(lattice)
# see what size clusters by each cutoff value 
bwplot(sizes~cutoff, data=cluster_sizes)

# pick a cutoff value 











