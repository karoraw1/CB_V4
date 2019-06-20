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


## Test Data is first 10 rows
start_time <- Sys.time()
### Calculate NTD
MNTD_w <- comdistnt.par(tsv.data[1:10,], phydist, abundance.weighted = TRUE, cores=1)
write(MNTD_w, stdout())
save(MNTD_w, file="../otu_data/dispersal_selection_data/test_ntd.RData")
### Calculate NTI
mntd_scores <- ses.comdistnt2(tsv.data[1:10,], phydist, method = "quasiswap", strata = NULL, abundance.weighted = TRUE, runs = 5, cores=1)
write(mntd_scores$comdistnt.obs.z, stdout())
save(mntd_scores, file="../otu_data/dispersal_selection_data/test_ses_nti.RData")
end_time <- Sys.time()
end_time - start_time
## Calculate MPD
mpd_w <- ses.mpd.par(tsv.data[1:10,], phydist, null.model = "taxa.labels", abundance.weighted = TRUE, cores = 1)
save(mpd_w, file="../otu_data/dispersal_selection_data/test_ses_mpd.RData")
end_time <- Sys.time()
end_time - start_time


## Real Data
start_time <- Sys.time()
MNTD_w <- comdistnt.par(tsv.data, phydist, abundance.weighted = TRUE, cores=24)
save(MNTD_w, file="../otu_data/dispersal_selection_data/ntd.RData")
end_time <- Sys.time()
end_time - start_time
mntd_scores <- ses.comdistnt2(tsv.data, phydist, method = "quasiswap", strata = NULL, abundance.weighted = TRUE, runs = 500, cores=24)
save(mntd_scores, file="../otu_data/dispersal_selection_data/ses_nti.RData")
end_time <- Sys.time()
end_time - start_time
mpd_w <- ses.mpd.par(tsv.data, phydist, null.model = "taxa.labels", abundance.weighted = TRUE, cores = 24)
save(mpd_w, file="../otu_data/dispersal_selection_data/ses_mpd.RData")
end_time <- Sys.time()
end_time - start_time


