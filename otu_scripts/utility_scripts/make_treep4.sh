#!/bin/bash

#SBATCH --job-name=mt4
#SBATCH --time=18:00:00
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=24
#SBATCH --mem=100G
#SBATCH --partition=parallel
#SBATCH --exclusive
#SBATCH --mail-type=END 
#SBATCH --mail-user=karoraw1@jhu.edu
#SBATCH --error=mt4.err
#SBATCH --output=mt4.out

LOC_DIR=~/scratch/CB_V4/otu_data/tree_data 
NEW_TD=$LOC_DIR/full_tree

mkdir -p $NEW_TD

cd $NEW_TD

QUERY=query_high_abund

source activate otu_caller

tr "[ -%,;\(\):=\.\\\[]\"\']" "_" < ../$QUERY.fasta > $QUERY.fasta

cmalign --cpu 24 --dna -o $QUERY.align.sto --sfile cm_scores.txt --outformat Pfam ../16S_bacteria.cm $QUERY.fasta

seqmagick convert $QUERY.align.sto $QUERY.align.fasta

raxmlHPC-PTHREADS -T 24 -m GTRGAMMA -s $QUERY.align.fasta -n $QUERY.ref.tre -f d -p 12345 
raxmlHPC-PTHREADS -T 24 -m GTRGAMMA -f I -t RAxML_bestTree.$QUERY.ref.tre -n root.$QUERY.ref.tre
raxmlHPC-PTHREADS -T 24 -m GTRGAMMA -f J -p 12345 -t RAxML_rootedTree.root.$QUERY.ref.tre -n conf.root.$QUERY.ref.tre -s $QUERY.align.fasta
raxmlHPC-PTHREADS -T 24 -f x -p 12345 -t RAxML_rootedTree.root.$QUERY.ref.tre -s $QUERY.align.fasta -m GTRGAMMA -n $QUERY.distmat

pplacer -o $QUERY.align.jplace -p -c ../hug_tol.refpkg $QUERY.align.fasta
guppy to_csv --point-mass --pp -o $QUERY.align.csv $QUERY.align.jplace
guppy fat --node-numbers --point-mass --pp -o $QUERY.align.phyloxml $QUERY.align.jplace
guppy distmat -o $QUERY.align.dist.tab $QUERY.align.jplace






