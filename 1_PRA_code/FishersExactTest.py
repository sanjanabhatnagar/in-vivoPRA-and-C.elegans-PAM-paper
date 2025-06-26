import scipy.stats as stats
import pandas as pd
import os
import sys

CL_count_df = sys.argv[1] # This consists the output table from the Splicing_Classifier_PRA_PSIcal.py script where the intronic cis-element, barcode and counts in all categories are present.
CL_file = os.path.basename(CL_count_df)
CL_file_name = os.path.splitext(CL_file)[0] 
output_file = f'{CL_file_name}_FishersExact.csv'

p_val_dic = {}
ACTG_Skip = 2883 # Specifying total count of spliced-out/exon-skipped isform for reference reporter
ACTG_Inc = 201  # Specifying total count of spliced-in/exon-included isform for reference reporter

i = 0
while i < len(CL_count_df):
    test_table = [[ACTG_Inc, ACTG_Skip],
                 [CL_count_df.iloc[i, 1], CL_count_df.iloc[i, 2]]]
    p_value = stats.fisher_exact(test_table)
    p_val_dic.update({CL_count_df.iloc[i, 0]: p_value})
    i = i + 1

CL_count_df['p-value'] = CL_count_df['Unnamed: 0'].map(p_val_dic)
CL_count_df.to_csv(output_file)

