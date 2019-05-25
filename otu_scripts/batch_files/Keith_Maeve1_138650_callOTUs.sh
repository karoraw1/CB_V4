#!/bin/bash

#SBATCH
#SBATCH --job-name=Keith_Maeve1_138650_DADA2
#SBATCH --time=5:00:00
#SBATCH --ntasks=48
#SBATCH --cpus-per-task=1
#SBATCH --mem=500G
#SBATCH --partition=lrgmem
#SBATCH --exclusive
#SBATCH --mail-type=END
#SBATCH --mail-user=karoraw1@jhu.edu
#SBATCH --error=DADA2_Keith_Maeve1_138650.err
#SBATCH --output=DADA2_Keith_Maeve1_138650.out

module load R/3.5.1

SEQ_ID=Keith_Maeve1_138650
BASE_IN=/home-3/karoraw1@jhu.edu/work/sprehei1/Keith_Files/Processed_data_group
BASE_OUT=/home-3/karoraw1@jhu.edu/scratch/CB_V3/otu_data/dada2_outputs; mkdir -p $BASE_OUT/$SEQ_ID;
SUFF1=_F_filt
SUFF2=_R_filt
SAMSPLIT=$SUFF1
THREADS=48
SCRIPTS_=~/scratch/CB_V3/otu_scripts/utility_scripts

Rscript $SCRIPTS_/fullDADApipe_PE.R $BASE_IN $BASE_OUT $SEQ_ID $SUFF1 $SUFF2 $SAMSPLIT $THREADS
