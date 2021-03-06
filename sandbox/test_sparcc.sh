conda deactivate
conda activate sparCC

SPARCC=/Volumes/KeithSSD/CB_V4/sandbox/sparcc
EXAMPLE=$SPARCC/example
PVALS=$EXAMPLE/pvals
BASE_COR=$EXAMPLE/basis_corr
N_BOOTS=3
BOOT_NUM=2

python $SPARCC/SparCC.py $EXAMPLE/fake_data.txt -i 5 --cor_file=$BASE_COR/cor_sparcc_test.out

python $SPARCC/MakeBootstraps.py $EXAMPLE/fake_data.txt -n $N_BOOTS -t test_permutation_#.txt -p $PVALS/

for i in $(seq 0 $BOOT_NUM); do 
    echo $i;
    python $SPARCC/SparCC.py $PVALS/test_permutation_${i}.txt -i 4 --cor_file=$PVALS/test_perm_cor_${i}.txt
done

python $SPARCC/PseudoPvals.py $BASE_COR/cor_sparcc_test.out $PVALS/test_perm_cor_#.txt 3 -o $PVALS/test_pvals.two_sided.txt -t two_sided
