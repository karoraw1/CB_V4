
setwd("/Volumes/KeithSSD/CB_V4/otu_data/WaterQualityData/cb33_envmodel")
      
mdata <- read.delim("cb33_mdata.tsv", row.names=1, colClasses=c("character", "numeric", "numeric", "numeric", "numeric", "numeric", "factor", "factor", "numeric", "numeric", "numeric"))
comm <- read.delim("cb33comm.tsv", row.names=1, colClasses=c(c('character'), rep('numeric', 1274)))
distmat <- read.delim("cb33distmat.tsv", row.names=1, colClasses=c(c('character'), rep('numeric', 99)))

library(ade4)
library(adegraphics)
library(adespatial)
library(vegan)
library(vegan3d)
library(MASS)
library(ellipse)
library(FactoMineR)
library(rrcov)

hist(mdata$observed_otus, breaks=15)
comm2 = comm[mdata$observed_otus > 100, ]
distmat2 = distmat[mdata$observed_otus > 100, mdata$observed_otus > 100]
mdata2 = mdata[mdata$observed_otus > 100, ]

spe.hel = decostand(comm2, 'hellinger')

mdata_env = mdata2[,colnames(mdata2)[1:5]]
colnames(mdata_env) = c("depth", "do", "wtemp", "pH", "sal")

# basic 
spechem.physio2 <- rda(spe.hel ~ do + wtemp + pH + sal, data = mdata_env)
spechem.physio2 <- rda(distmat2 ~ do + wtemp + pH + sal, data = mdata_env)
anova(spechem.physio2, permutations = how(nperm = 999))
anova(spechem.physio2, permutations = how(nperm = 999), by = "axis")
RsquareAdj(spechem.physio2)$r.squared
RsquareAdj(spechem.physio2)$adj.r.squared
vif.cca(spechem.physio2)
R2a.all = RsquareAdj(spechem.physio2)$adj.r.squared
forward.sel(distmat2, mdata_env, adjR2thresh = R2a.all)

adist2 = as.dist(distmat2)
anv_result = adonis(adist2 ~ do + wtemp + pH + sal, data = mdata_env)
anv_result2 = adonis(spe.hel ~ do + wtemp + pH + sal, data = mdata_env, method='bray')
bray.env.cap <- capscale(comm2 ~ do + wtemp + pH + sal, data = mdata_env, distance = "bray")
uf.env.cap <- capscale(adist2 ~ do + wtemp + pH + sal, data = mdata_env)
anova(bray.env.cap, permutations = how(nperm = 999))
anova(uf.env.cap, permutations = how(nperm = 999))

# interpret this plot? 
plot(spechem.physio2, display=c('wa', 'cn'))

# look at the larger data set with RDA 

# partial out the spatial component ? 

# pg 362 to combine spatial and environmental components



# an LDA of the clusters we found in the previous step? 


