library(vegan)
library(foreach)
library(doParallel)

#setwd("/Volumes/KeithSSD/CB_V4/otu_scripts")

numCores = detectCores(TRUE)

raup_crick_stegen <- function(spXsite, reps) {
	
	registerDoParallel(numCores)
	## count number of sites and total species richness across all plots (gamma)
	n_sites<-nrow(spXsite)
	gamma<-ncol(spXsite)
		
	##make the spXsite matrix into a pres/abs.
	spXsite_pa <- ceiling(spXsite/max(spXsite))

	# make spXsite into relative abundances
	spXsite_ra <- t(scale(t(spXsite), center = FALSE, scale = rowSums(spXsite)))

	##create an occurrence vector- used to give more weight to widely distributed species in the null model:
	occur<-apply(spXsite_pa, MARGIN=2, FUN=sum)

	## determine how many unique species richness values are in the dataset
	alpha_levels<-apply(spXsite_pa, MARGIN=1, FUN=sum)
	# how many individuals to sample per null community
	sample_totals = rowSums(spXsite)
	# the relative abundance of each species across all samples
	species_totals = colSums(spXsite_ra)

	##make_null:

	bootstrap_fxn <- function(pair_idx) {
			##two empty null communities of size gamma:
			a1x = unlist(pair_idx)[[1]]; 
			a2x = unlist(pair_idx)[[2]]; 

			com1<-rep(0,gamma); com2<-rep(0,gamma)
			
			##add alpha1 number of species to com1, weighting by species occurrence frequencies:
			com1_indexes = sample(1:gamma, alpha_levels[a1x], replace=FALSE, prob=occur)
			com1[com1_indexes]<-max(sample_totals)
			## draw species into each selected species in proportion to their RA
			com1_sp_totals = round(species_totals*com1, digits=0)
			## draw until observed population size is obtained
			com1_redraw = rrarefy(com1_sp_totals, sample_totals[[a1x]])

			##same for com2:
			com2_indexes = sample(1:gamma, alpha_levels[a2x], replace=FALSE, prob=occur)
			com2[com2_indexes]<-max(sample_totals)
			## draw species into each selected species in proportion to their RA
			com2_sp_totals = round(species_totals*com2, digits=0)
			## draw until observed population size is obtained
			com2_redraw = rrarefy(com2_sp_totals, sample_totals[[a2x]])

			##bray curtis between random samplings
			this_dist <- vegdist(rbind(com1_redraw, com2_redraw))[1]
			return(this_dist)
			}

	# a table for each element in the distance matrix
	alpha_table<-data.frame(matrix(nrow=sum(dim(spXsite)[1]:1), ncol=2))
	names(alpha_table)<-c("Site1", "Site2")
	##null_array will hold the actual arrays of null distances for every rep for each iteration 
	null_array<-matrix(nrow=dim(alpha_table)[1], ncol=reps)
    
    col_count = 1
	for (a1 in 1:n_sites) {
		for(a2 in a1:n_sites){
			alpha_table[col_count, 'Site1']<-a1
			alpha_table[col_count, 'Site2']<-a2
			#increment the counter for the columns of the alpha table/ elements of the null array
			col_count<-col_count+1
		}
	}

	null_array <- foreach (i = 1:dim(alpha_table)[1], .combine=rbind) %dopar% {
		replicate(reps, bootstrap_fxn(alpha_table[i,]), )
	}
	rownames(null_array) <- 1:dim(null_array)[1]

	##create a new column with both alpha levels to match on:
	alpha_table$matching<-paste(alpha_table[,1], alpha_table[,2], sep="_")

	#####################
	##do the test:

	##build a site by site matrix for the results, with the names of the sites in the row and col names:
	results<-matrix(data=NA, nrow=n_sites, ncol=n_sites, dimnames=list(row.names(spXsite), row.names(spXsite)))

	#for each pair of sites (duplicates effort now to make a full matrix instead of 
	#a half one- but this part should be minimal time as compared to the null model building)
	for(i in 1:n_sites){
		for(j in i:n_sites){
			##how many species are shared between the two sites:
			#n_shared_obs<-sum((spXsite[i,]+spXsite[j,])>1)
			n_shared_obs<-vegdist(rbind(spXsite[i,], spXsite[j,]))[1]
			
			##match against the alpha table- row index identifies which element of the null array contains the correct null distribution for the observed combination of alpha values:
			null_index<-which(alpha_table$matching==paste(i, j, sep="_"))
			
			##how many null observations is the observed value tied with?
			num_exact_matching_in_null<-sum(null_array[null_index,] == n_shared_obs)
			
			##how many null values are bigger than the observed value?
			num_greater_in_null<-sum(null_array[null_index,]>n_shared_obs)
			
			rc<-(reps-num_greater_in_null + (num_exact_matching_in_null/2))/reps
			rc<-(rc-.5)*2
			
			##store the metric in the results matrix:
			results[i,j]<-round(rc, digits=2)
		}
	}
	stopImplicitCluster()	
	print(c("Reps: ", reps))
	print(c("Dim: ", dim(spXsite)[1]))
	print(c("Pairwise: ", sum(dim(spXsite)[1]:1)))
	print(c("Iterations: ", reps*sum(dim(spXsite)[1]:1)))
	return(results)
}

real_spXsite <- read.delim("../otu_data/dispersal_selection_data/final_rarefied_table.tsv", row.names=1)
real_reps <- 999
real_out <- "../otu_data/dispersal_selection_data/raup_crick_data.tsv"

real_res = raup_crick_stegen(real_spXsite, real_reps)
write.table(real_res, file = real_out, sep = "\t")


