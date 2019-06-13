#!/bin/bash

#SBATCH
#SBATCH --job-name=sparcc
#SBATCH --time=12:00:00
#SBATCH --nodes=1
#SBATCH --mem=100G
#SBATCH --partition=shared
#SBATCH --ntasks-per-node=1
#SBATCH --cpus-per-task=24
#SBATCH --mail-type=END
#SBATCH --mail-user=karoraw1@jhu.edu
#SBATCH --error=sparcc_.err
#SBATCH --output=sparcc_.out

#source activate SparCC
#conda install python-dateutil=2.4.2
#conda install numpy=1.9.2 pandas=0.16.2
#conda install libcxx

ml python/2.7-anaconda53
source activate SparCC
source ~/.bash_profile

SPARCC=$(dirname $(pwd))/sandbox/sparcc
OUTDIR=$(dirname $(pwd))/otu_data/sparcc_data
OTU_TABLE=$OUTDIR/filtered_otu_table.txt

#python $SPARCC/SparCC.py $OTU_TABLE --cor_file=$OUTDIR/sparcc_corr.out --cov_file=$OUTDIR/sparcc_cov.out

#N_BOOTS=100
PERM_=$OUTDIR/perms
#mkdir -p $PERM_
#python $SPARCC/MakeBootstraps.py $OTU_TABLE -n $N_BOOTS -t perm_#.txt -p $PERM_/

PCOR_=$OUTDIR/perm_corrs
mkdir -p $PCOR_
BOOT_NUM=99
perl -e 'for(0..99){print "$_\n"}' > itercnt.txt
ml parallel
echo $(which python)
parallel  --jobs 24 "python $SPARCC/SparCC.py $PERM_/perm_{}.txt --cor_file=$PCOR_/perm_cor_{}.txt" :::: itercnt.txt
rm itercnt.txt

python $SPARCC/PseudoPvals.py $OUTDIR/sparcc_corr.out $PCOR_/perm_cor_#.txt 100 -o $OUTDIR/test_pvals.two_sided.txt -t two_sided
