library(vegan)
library(ape)
library(picante)
library(MicEco)
library(doSNOW)

### Calculate NTR and NTI

## Test Data
start_time <- Sys.time()
tsv.data <- read.delim("../otu_data/clustered_sequences/test_abundances.tsv", row.names=1)
phydf <- read.delim("../otu_data/clustered_sequences/fixed_pplacer_distmat.tsv", row.names=1)
phydist = as.matrix(phydf)
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


