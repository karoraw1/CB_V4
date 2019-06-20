#!/bin/bash

#SBATCH
#SBATCH --job-name=correlg
#SBATCH --time=12:00:00
#SBATCH --nodes=1
#SBATCH --mem=100G
#SBATCH --partition=parallel
#SBATCH --ntasks-per-node=1
#SBATCH --cpus-per-task=24
#SBATCH --mail-type=END
#SBATCH --mail-user=karoraw1@jhu.edu
#SBATCH --error=correlg.err
#SBATCH --output=correlg.out


# takes 4-5 hours on 370 x 1900 otu table

ml R/3.5.1
Rscript corellogram.R
