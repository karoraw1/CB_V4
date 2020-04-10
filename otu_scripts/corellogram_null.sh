#!/bin/bash

#SBATCH
#SBATCH --job-name=correlgn
#SBATCH --time=12:00:00
#SBATCH --nodes=1
#SBATCH --mem=100G
#SBATCH --partition=lrgmem,shared,skylake,parallel
#SBATCH --ntasks-per-node=1
#SBATCH --cpus-per-task=1
#SBATCH --mail-type=END
#SBATCH --mail-user=karoraw1@jhu.edu
#SBATCH --error=correlg.n.err
#SBATCH --output=correlg.n.out

ml R/3.5.1
Rscript corellogram_null.R
