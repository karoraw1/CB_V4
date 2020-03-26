#!/bin/bash

#SBATCH
#SBATCH --job-name=rc_mod
#SBATCH --time=24:00:00
#SBATCH --nodes=1
#SBATCH --mem=100G
#SBATCH --partition=parallel
#SBATCH --ntasks-per-node=1
#SBATCH --cpus-per-task=24
#SBATCH --mail-type=ALL
#SBATCH --mail-user=karoraw1@jhu.edu
#SBATCH --error=rc_mod_par.err
#SBATCH --output=rc_mod_par.out

# 
ml R/3.5.1

date

Rscript RaupCrick_Final.R

date
