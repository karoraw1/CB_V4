#!/bin/bash

#SBATCH
#SBATCH --job-name=sprehei1_149186_DADA2
#SBATCH --time=5:00:00
#SBATCH --nodes=1
#SBATCH --mem=50G
#SBATCH --partition=gpup100
#SBATCH --gres=gpu:2
#SBATCH --ntasks-per-node=2
#SBATCH --cpus-per-task=6
#SBATCH --exclusive
#SBATCH --mail-type=END
#SBATCH --mail-user=karoraw1@jhu.edu
#SBATCH --error=DADA2_sprehei1_149186.err
#SBATCH --output=DADA2_sprehei1_149186.out

module load R/3.5.1

SEQ_ID=sprehei1_149186
BASE_IN=/home-3/karoraw1@jhu.edu/work/sprehei1/Keith_Files/Processed_data_group
BASE_OUT=/home-3/karoraw1@jhu.edu/scratch/CB_V3/otu_data/dada2_outputs; mkdir -p $BASE_OUT/$SEQ_ID;
SUFF1=_F_filt
SUFF2=_R_filt
SAMSPLIT=$SUFF1
THREADS=24
SCRIPTS_=~/scratch/CB_V3/otu_scripts/utility_scripts

Rscript $SCRIPTS_/fullDADApipe_PE.R $BASE_IN $BASE_OUT $SEQ_ID $SUFF1 $SUFF2 $SAMSPLIT $THREADS
