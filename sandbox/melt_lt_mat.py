import pandas as pd
import numpy as np
from statsmodels.stats.multitest import multipletests

corrs_file = "/Volumes/KeithSSD/CB_V4/otu_data/sparcc_data/sparcc_corr.out"
df = pd.read_csv(test_file, sep="\t", index_col=0)

#fake_pvals = np.random.uniform(low=0, high=4, size=df.shape)/4
#p_df = pd.DataFrame(fake_pvals, index=df.index, columns=df.columns)

pval_file = "/Volumes/KeithSSD/CB_V4/otu_data/sparcc_data/test_pvals.two_sided.txt"
p_df = pd.read_csv(pval_file, sep="\t", index_col=0)

def melt_upper_triangle(df_, val_str):
    dfnan = df_.where(np.triu(np.ones(df_.shape)).astype(np.bool))
    melted_df = dfnan.stack().reset_index()
    melted_df.columns = ['OTU_1','OTU_2', val_str]
    melted_df2 = melted_df[melted_df['OTU_1'] != melted_df['OTU_2']]
    return melted_df2.set_index(['OTU_1', 'OTU_2'])

mpdf = melt_upper_triangle(p_df, 'p-value')
mdf = melt_upper_triangle(df, 'correlation')

fulldf = mdf.join(mpdf)

# pull total abundances
# pull taxonomy (order?)

reject, pvals_corrected = multipletests(fulldf['p-value'].values, alpha=0.05, method='fdr_bh')[:2]

thresholded = fulldf.loc[fulldf.index[reject], ['correlation']].reset_index()
corr_cutoff = abs(thresholded.correlation) > 0.5
thresholded_cutoff = thresholded[corr_cutoff]

thresholded_cutoff.to_csv("/Volumes/KeithSSD/CB_V4/otu_data/sparcc_data/test_correlations.txt", sep="\t", index=False)
