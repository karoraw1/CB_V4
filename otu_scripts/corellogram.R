library(vegan)

tsv.data <- read.delim("../otu_data/clustered_sequences/pplacer_abundances.tsv", row.names=1)
exp.comm =  tsv.data[, colSums(tsv.data != 0) > 0]

sim.comm <- permatswap(exp.comm, "quasiswap", times=3)$perm[[2]]

phydf <- read.delim("../otu_data/clustered_sequences/fixed_pplacer_distmat.tsv", row.names=1)
phydist <- as.matrix(phydf[colnames(exp.comm),colnames(exp.comm)])

exp.comm.dist = as.matrix(vegdist(t(exp.comm), method='bray'))
sim.comm.dist = as.matrix(vegdist(t(sim.comm), method="bray"))

start_time <- Sys.time()
exp.correlog = mantel.correlog(exp.comm.dist, D.geo=phydist, mult="BH", nperm=5)
end_time <- Sys.time()
end_time - start_time
start_time <- Sys.time()
exp.correlog = mantel.correlog(exp.comm.dist, D.geo=phydist, mult="BH", nperm=10)
end_time <- Sys.time()
end_time - start_time
start_time <- Sys.time()
exp.correlog = mantel.correlog(exp.comm.dist, D.geo=phydist, mult="BH", nperm=50)
sim.correlog <- mantel.correlog(sim.comm.dist, D.geo=phydist, mult="BH", nperm=50)
end_time <- Sys.time()
end_time - start_time

summary(exp.correlog)
summary(sim.correlog)
par()
layout(matrix(c(1,2), 1, 2, byrow = TRUE))
plot(exp.correlog, xlab="Phylogenetic distance"); grid();
plot(sim.correlog, xlab="Phylogenetic distance"); grid();


phydist2 <- as.dist(phydf[colnames(exp.comm),colnames(exp.comm)])

exp.comm.dist2 = vegdist(t(exp.comm), method='bray')
sim.comm.dist2 = vegdist(t(sim.comm), method="bray")

exp.mgram = mgram(exp.comm.dist2, phydist2, breaks=(1:34)/10, nperm = 50, nboot = 10, trace = TRUE)
sim.mgram = mgram(sim.comm.dist2, phydist2, breaks=(1:34)/10, nperm = 50, nboot = 10, trace = TRUE)
layout(matrix(c(1,2), 1, 2, byrow = TRUE))
plot(exp.mgram, pval = 0.05, xlab = "Phylo Distance", ylab = "Mantel r"); grid();
plot(sim.mgram, pval = 0.05, xlab = "Phylo Distance", ylab = "Mantel r"); grid();


exp.rc.dist = raupcrick(exp.comm, null = "r1", nsimul = 50, chase = FALSE)
norm.e.rc.d = (exp.rc.dist - .5) * 2

hist(norm.e.rc.d)
sum(norm.e.rc.d > 0.95) / sum(1:370)
sum(norm.e.rc.d < -0.95) / sum(1:370)
sum(abs(norm.e.rc.d) < 0.95) / sum(1:370)

