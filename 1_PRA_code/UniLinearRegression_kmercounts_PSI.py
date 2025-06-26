import pandas as pd
import numpy as np

# Open the kmer count matrix file
df_x = pd.read_csv('kmer_count.csv')

# This is only needed for 7mer datasets, otherwise can be commented out
# This is used to change the column names to kmers
import ast
with open('column_names_7merMatrix.txt', 'r') as file:
    content = file.read()
column_names = ast.literal_eval(content)
df_x.rename(columns=column_names, inplace=True)

# Changing the name of the kmer count matrix dataframe and adding PSI column to the matrix dataframe
standardized_kmer_matrix = df_x
df_PRA = pd.read_excel('AllPRA_PSI.xlsx', sheet_name=0)
df_PRA_PSI = df_PRA[['PRA','PSI']]

standardized_df = df_x.merge(df_PRA_PSI, left_on='PRA', right_on='PRA')
std_kmer_counts = np.std(df_x.iloc[:,1:], axis=0)

# Now only considering those kmers for regression analysis which are present in more than 5 PRA elements.
# Count non-zero values in each column
non_zero_counts = (standardized_df != 0).sum(axis=0)
# Filter columns based on the condition (at least 5 non-zero values)
filtered_df = standardized_df.loc[:, non_zero_counts >= 5]

# Dropping kmers with no occurence.
filtered_df.replace(0,np.nan,inplace=True)
filtered_df.dropna(axis=1, thresh=1, inplace=True)
filtered_df.dropna(axis=0, thresh=3, inplace=True)
filtered_df.replace(np.nan,0,inplace=True)
print(filtered_df.shape)

# Storing the filtered count matrix into a new file for UMAP purposes.
filtered_df.to_csv('kmerscountsAug2024_kmers_in>5PRA.csv')

# MAIN CODE THAT RUNS REGRESSION #
from sklearn.linear_model import LinearRegression
import scipy
from scipy import stats

# formatting the dataframe containing all information and converting it into a numpy array/matrix
standardized_df1 = pd.DataFrame(filtered_df.drop(['PSI','PRA'], axis=1))
standardized_df2 = standardized_df1.to_numpy()

# For getting the names of different kmers
kmer = list(standardized_df1.columns.values.tolist())
y = filtered_df['PSI']

model = LinearRegression(copy_X=True, fit_intercept=True, n_jobs=None, normalize=False)
standardized_coefficients=[]

with open('PRA_kmerRegression_kmersin_>5PRA_2024.tsv','a') as h4:
    h4.write('kmer\tIntercept\tCoefficient\tSlope\tRvalue\tR-squared\tpvalue\tStd_Err\tt_stat\tStandardised_Coefficient')

    for i in range(standardized_df2.shape[1]):
        X = standardized_df2[:, i]
        X = X.reshape(-1, 1)
        model.fit(X, y)
        beta = model.coef_[0]
        standardized_coefficient = beta / std_kmer_counts[i]

        standardized_coefficients.append(standardized_coefficient)
        slope, intercept, r_value, p_value, std_err = scipy.stats.linregress(standardized_df1.iloc[:,i], y)

        t_statistic = slope / std_err

        h4.write("%s\t" %kmer[i])
        h4.write("%f\t" %model.intercept_)
        h4.write("%f\t" %model.coef_[0])
        h4.write("%f\t" %slope)
        h4.write("%f\t" %r_value)
        h4.write("%f\t" %(r_value**2))
        h4.write('{:.20e}\t'.format(p_value))
        h4.write('{:.20e}\t'.format(std_err))
        h4.write("%f\t" %t_statistic)
        h4.write("%f\t\n" %standardized_coefficient)

# The output will be a tsv file with t_stat, std_error etc. for each kmer occuring in more than 5 PRA elements.
