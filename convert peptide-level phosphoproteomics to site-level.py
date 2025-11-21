#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Tue Oct 21 20:04:07 2025

@author: harrisbolus
"""
import pandas as pd
import numpy as np

#----------------------------params----------------------------------

input_file = {'workbook': '/Users/harrisbolus/Desktop/Research/Dr. Honkanen/PP5 MCF7 : HEK KO/PP5_summary.xlsx',
              'sheet': 'phos_protein corrected'}

output_file = '.../... site-level phosphoproteomics.xlsx'

uniprot_col = 'protein_accession'
modsite_col = 'modsites'
seq_col = 'centered_sequences'

# list columns here that should not be aggregated. This should include all columns aside from your log2 abundance columns.
non_aggregative_cols = ['protein_name',
'protein_ref',
'gene_id',
'protein_accession',
'protein_description',
'modsites',
'centered_sequences',
'tp',
'up',
'um',
'min_flr']

# list columns here that contain your log2 abundance data
log_cols = ['log2 area(HEK_WT)_1',
'log2 area(HEK_WT)_2',
'log2 area(HEK_WT)_3',
'log2 area(HEK_WT)_4',
'log2 area(HEK_PP5-KO)_1',
'log2 area(HEK_PP5-KO)_2',
'log2 area(HEK_PP5-KO)_3',
'log2 area(HEK_PP5-KO)_4',
'log2 area(MCF7_WT)_1',
'log2 area(MCF7_WT)_2',
'log2 area(MCF7_WT)_3',
'log2 area(MCF7_WT)_4',
'log2 area(MCF7_PP5-KO)_1',
'log2 area(MCF7_PP5-KO)_2',
'log2 area(MCF7_PP5-KO)_3',
'log2 area(MCF7_PP5-KO)_4']

# hopefully your log2 abundance columns have a prefix that can be dropped to make appropriately-named raw abundance columns. Please supply them here
log_col_prefix = 'log2 '

# these will simply be dropped and you will need to recalculate them
stats_cols = ['log2 FC',
'p-value']

#--------------------------------------------------------------------------

# read data in, drop peptide-level statistics
phosdf = pd.read_excel(input_file['workbook'],
                       sheet_name = input_file['sheet']).drop(stats_cols, axis=1)

# split by rows into separate dfs for singly and multiply phosphorylated peptides
singly_phosphorylated_df = phosdf[~phosdf[modsite_col].str.contains(':')]
multiply_phosphorylated_df = phosdf[phosdf[modsite_col].str.contains(':')]

# separate each row with multiple modsites into multiple rows with one modsite each, store in temp
temp = []
for uniprot in set(multiply_phosphorylated_df[uniprot_col]):
    for modsites in set(multiply_phosphorylated_df[multiply_phosphorylated_df[uniprot_col] == uniprot][modsite_col]):
        subset = multiply_phosphorylated_df[(multiply_phosphorylated_df[uniprot_col] == uniprot) & (multiply_phosphorylated_df[modsite_col] == modsites)]

        sequences = subset[seq_col].item().split(',')
        modsites = modsites.split(':')

        for seq, mod in zip(sequences, modsites):
            subset.loc[:, modsite_col] = mod
            subset.loc[:, seq_col] = seq
            temp.append(subset.copy(deep=True))

# add these new rows back to the singly phosphorylated df
singly_phosphorylated_df = pd.concat([singly_phosphorylated_df, *temp], axis=0)

# unlog
singly_phosphorylated_df.loc[:,'unique id'] = singly_phosphorylated_df[uniprot_col] + singly_phosphorylated_df[modsite_col]
abd_cols = [i.replace(log_col_prefix,'') for i in log_cols]
singly_phosphorylated_df[log_cols] = 2**singly_phosphorylated_df[log_cols]
singly_phosphorylated_df = singly_phosphorylated_df.rename(dict(zip(log_cols, abd_cols)), axis=1)

# aggregate
non_aggregative_cols = {i: 'first' for i in non_aggregative_cols}
non_aggregative_cols.update({i:'sum' for i in abd_cols})
collapsed_df = singly_phosphorylated_df.groupby('unique id', as_index=False).agg(non_aggregative_cols)

# recalculate log2 abundances, calling them 0 where abundance is 0
collapsed_df[log_cols] = np.where(collapsed_df[abd_cols] != 0, np.log2(collapsed_df[abd_cols]), 0)
collapsed_df = collapsed_df.drop(abd_cols, axis=1)
collapsed_df.to_excel(output_file, index=False)
