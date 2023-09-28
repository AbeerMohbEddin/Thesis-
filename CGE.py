#!/usr/bin/env python
# coding: utf-8

# In[2]:


import os 

# Change the working directory
os.chdir("C:/Users/abeer/Desktop/MScWD")

# Get the current directory and list all files in it
current_directory = os.getcwd()
files = os.listdir(current_directory)

# Print the files to check
print(files)


# In[3]:


import pandas as pd
import numpy as np
from scipy.stats import rankdata
from statsmodels.stats.multitest import multipletests
from itertools import permutations
import random


# In[4]:


gene_expression_data = pd.read_csv("DeAbsent_cpmNorm_Counts_NP11.csv", index_col = 0)
# Rename the index column
gene_expression_data.index.name = 'GeneID'

gene_expression_data


# In[5]:


array = np.array(gene_expression_data)
array


# In[6]:


index_values = gene_expression_data.index.values
index_column_name = gene_expression_data.index.name
mapped_array = np.hstack((index_values.reshape(-1, 1), array))
mapped_array


# In[7]:


# Extract the expression columns from the mapped array
expression_columns = mapped_array[:, 1:].astype(float)

# Rank the expression columns
ranked_columns = np.apply_along_axis(lambda x: np.argsort(-x).argsort() + 1, axis=0, arr=expression_columns)

# Create a new array by combining the index column and ranked expression columns
ranked_array = np.hstack((mapped_array[:, 0].reshape(-1, 1), ranked_columns))

# Get the number of columns in the array
num_columns = ranked_array.shape[1]

ranked_array


# In[8]:


ranked_df = pd.DataFrame(ranked_array)

# Rename the numeric column names
ranked_df = ranked_df.rename(columns={0: 'GeneID'}) 
ranked_df


# In[9]:


ranked_array2 = ranked_array[:, 1:].astype(float)
rank_prod = np.prod(ranked_array2, axis=1)
rank_prod


# In[10]:


rankProductdf = pd.DataFrame(rank_prod)

# Rename the numeric column names
rankProductdf = rankProductdf.rename(columns={0: 'RP value'}) 
rankProductdf


# In[11]:


# Select the column 1 using iloc
GeneID_column = ranked_df.iloc[:, 0]

# Select the column 1 using iloc
RP_column = rankProductdf.iloc[:, 0]

# Create a new DataFrame by concatenating the two columns
new_df = pd.concat([GeneID_column, RP_column], axis=1)
new_df


# In[12]:


# Sort the DataFrame by the 'RP value' column in ascending order
RP_Ascend_sort = new_df.sort_values('RP value', ascending=True)
RP_Ascend_sort


# In[13]:


# Merge the DataFrames based on the 'ID' column
merged_Ascend_df = pd.merge(RP_Ascend_sort, gene_expression_data, on='GeneID')
merged_Ascend_df 


# In[14]:


Ascen_sorted_geneexpressionData = merged_Ascend_df.drop('RP value', axis=1)

# Convert 'GeneID' to the index column
Ascen_sorted_geneexpressionData.set_index('GeneID', inplace=True)
Ascen_sorted_geneexpressionData


# In[15]:


# Calculate the mean expression without considering the 'GeneID' column
merged_Ascend_df['Mean Expression'] = merged_Ascend_df.drop('RP value', axis=1).mean(axis=1)
#Ascend_mean_expression = merged_Ascend_df.drop('RP value', axis=0).mean()
merged_Ascend_df


# In[16]:


# Drop columns ''RP value', and 'Mean Expression' from the DataFrame
merged_Ascend_df_dropped = merged_Ascend_df.drop(['RP value', 'Mean Expression'], axis=1)

# Calculate the standard deviation by gene (row-wise)
merged_Ascend_df_dropped['std_dev'] = merged_Ascend_df_dropped.std(axis=1)
merged_Ascend_df_dropped


# In[17]:


columnRP_Mean = merged_Ascend_df[['GeneID', 'RP value', 'Mean Expression']]
column_StdDev = merged_Ascend_df_dropped[['std_dev']]
# Create a new DataFrame by concatenating the two columns
Ascen_RP_SortedTable1 = pd.concat([columnRP_Mean, column_StdDev], axis=1)

# Set the first column as index
#Ascen_RP_SortedTable1.set_index(Ascen_RP_SortedTable1.columns[0], inplace=True)
Ascen_RP_SortedTable1


# In[18]:


# Select the column 1 using iloc
rp = Ascen_RP_SortedTable1.iloc[:, 1]
rp_array = np.array(rp)
rp_array


# In[19]:


def calculate_rank_product(gene_expression_data):
       
    # Convert DataFrame to array
    array_data = gene_expression_data.values
    index_values = gene_expression_data.index.values
    index_column_name = gene_expression_data.index.name
    mapped_array = np.hstack((index_values.reshape(-1, 1), array_data))

    # Extract the expression columns from the mapped array
    expression_columns = mapped_array[:, 1:].astype(float)

    # Rank the expression columns
    ranked_columns = np.apply_along_axis(lambda x: np.argsort(-x).argsort() + 1, axis=0, arr=expression_columns)

    # Create a new array by combining the index column and ranked expression columns
    ranked_array = np.hstack((mapped_array[:, 0].reshape(-1, 1), ranked_columns))

    # Get the number of columns in the array
    num_columns = ranked_array.shape[1]

    # Create a copy of the original array
    sorted_array = ranked_array.copy()

    sorted_array = sorted_array[:, 1:].astype(float)
    #Rprod = np.prod(sorted_array, axis=1)
    Rank_Product = np.prod(sorted_array, axis=1)

    return Rank_Product


def permutation_test(gene_expression_data, Rank_Product, num_permutations, seed=None):
    if seed is not None:
        np.random.seed(seed)
        random.seed(seed)

    permuted_rps = np.zeros((num_permutations, len(Rank_Product)))

    for i in range(num_permutations):
        permuted_data = gene_expression_data.copy()
        for column in gene_expression_data:
            permuted_data[column] = np.random.permutation(permuted_data[column])

        permuted_rank_product = calculate_rank_product(permuted_data)
        permuted_rps[i] = permuted_rank_product

    return permuted_rps


# In[20]:


Rank_Product = calculate_rank_product(Ascen_sorted_geneexpressionData)
Rank_Product


# In[22]:


# Set the random seed
seed_value = 42

# Call permutation_test with the seed
permutedRP = permutation_test(Ascen_sorted_geneexpressionData, Rank_Product, num_permutations=1000, seed=seed_value)
permutedRP


# In[23]:


def calculate_p_values(rank_product, permuted_values, num_permutations):
    num_genes = len(rank_product)
    p_values = np.zeros(num_genes)

    for i in range(num_genes):
        # Calculate the observed rank product for the current gene
        observed_rank_product = rank_product[i]

        # Calculate the number of permuted rank products that are greater than the observed rank product
        num_greater = np.sum(permuted_values[:, i] > observed_rank_product)

        # Calculate the p-value for the current gene by dividing the number of permuted rank products by the total number of permutations
        p_value = (num_greater + 1) / (num_permutations + 1)     

        # Invert the p-value to represent lower consistency levels as higher values
        inverted_p_value = 1 - p_value

        p_values[i] = inverted_p_value

    return p_values


# In[24]:


num_permutations=1000
p_values = calculate_p_values(Rank_Product, permutedRP, num_permutations)
p_values


# In[25]:


Ascen_pval_df = pd.DataFrame(p_values)

# Rename the column
Ascen_pval_df = Ascen_pval_df.rename(columns={0: 'P value'})
Ascen_pval_df


# In[26]:


# Perform the merge based on the index column
Ascen_RP_Sorted_PvalTable = pd.merge(Ascen_RP_SortedTable1, Ascen_pval_df, left_index=True, right_index=True)

# Convert 'Column1' to the index column
Ascen_RP_Sorted_PvalTable.set_index('GeneID', inplace=True)
#RP_Sorted_PvalTable.to_csv('GW8PE_RP_Sorted_PvalTable.csv')
Ascen_RP_Sorted_PvalTable


# In[27]:


# Set the P-value cutoff and calculate the 75th percentile of Mean Expression
p_value_cutoff = 0.05

# Filter genes based on P-value cutoff
SG_bypval_Ascend = Ascen_RP_Sorted_PvalTable[Ascen_RP_Sorted_PvalTable['P value'] <= p_value_cutoff]
#SG_bypval_df.to_csv('t1.csv')
SG_bypval_Ascend 


# In[28]:


# Calculate the mean expression level for each gene
mean_expression_levels = SG_bypval_Ascend['Mean Expression']
mean_expression_levels = np.array(mean_expression_levels)

# Sort the genes based on their mean expression levels (descending order)
sorted_genes = np.argsort(mean_expression_levels)[::-1]

# Select the top genes with the highest consistent expression levels (e.g., top 10%)
top_genes_percentage = 0.1
num_top_genes = int(len(sorted_genes) * top_genes_percentage)
selected_high_consistent_genes = sorted_genes[:num_top_genes]

# Determine the minimum mean expression level among the selected high consistent genes
max_mean_exp_TopGenes_Ascend = np.max(mean_expression_levels[selected_high_consistent_genes])

min_mean_exp_TopGenes_Ascend = np.min(mean_expression_levels[selected_high_consistent_genes])

# Select the bottom genes with the lowest consistent expression levels (e.g., bottom 10%)
bottom_genes_percentage = 0.1
num_bottom_genes = int(len(sorted_genes) * bottom_genes_percentage)
selected_low_consistent_genes = sorted_genes[-num_bottom_genes:]

# Determine the maximum mean expression level among the selected low consistent genes
min_mean_exp_BottomGenes_Ascend = np.min(mean_expression_levels[selected_low_consistent_genes])

# Determine the maximum mean expression level among the selected low consistent genes
max_mean_exp_BottomGenes_Ascend = np.max(mean_expression_levels[selected_low_consistent_genes])

# Print the minimum and maximum mean expression levels
print("Maximum mean expression level in top genes:", max_mean_exp_TopGenes_Ascend)
print("Minimum mean expression level in top genes:", min_mean_exp_TopGenes_Ascend)
print("Maximum mean expression level in bottom genes:", max_mean_exp_BottomGenes_Ascend)
print("Minimum mean expression level in bottom genes:", min_mean_exp_BottomGenes_Ascend)


# In[29]:


# Set the P-value cutoff and calculate the 75th percentile of Mean Expression
p_value_cutoff = 0.05
mean_Min = max_mean_exp_BottomGenes_Ascend        #7.835853527749999              #7.197171037

# Filter genes based on P-value cutoff
filtered_df1 = SG_bypval_Ascend[SG_bypval_Ascend['P value'] <= p_value_cutoff]

# Categorize high consistent genes based on Mean Expression above the mean_Min
high_consistent_genes = filtered_df1[filtered_df1['Mean Expression'] > mean_Min]
high_consistent_genes


# In[30]:


import pandas as pd
import mygene

def convert_ensembl_to_symbols(ensembl_ids):
    mg = mygene.MyGeneInfo()

    # Retrieve gene annotations for Ensembl IDs
    results = mg.querymany(ensembl_ids, scopes='ensembl.gene', fields='symbol', species='human')

    gene_symbols = []
    for result in results:
        if 'symbol' in result:
            ensembl_id = result['query']
            gene_symbol = result['symbol']
            gene_symbols.append((ensembl_id, gene_symbol))

    return gene_symbols

# Load the Ensemble gene IDs
ensembl_ids = high_consistent_genes.index
gene_symbols = convert_ensembl_to_symbols(ensembl_ids)

# Convert the list of tuples to a DataFrame
df = pd.DataFrame(gene_symbols, columns=['GeneID', 'Gene_Symbol'])
df


# In[31]:


# Convert the index column to a regular column
hcg_df = high_consistent_genes.reset_index()
hcg_df


# In[32]:


# Merge the DataFrames based on the first column (ID)
merged_df1 = pd.merge(df, hcg_df, on='GeneID')
merged_df1


# In[33]:


# Merge the DataFrames based on the first column (ID) and keep the unmerged rows
merged_df2 = pd.merge(df, hcg_df, on='GeneID', how='outer')
merged_df2


# In[34]:


# Sort the DataFrame by the 'RP value' column in ascending order
RP_Descend_sort = new_df.sort_values('RP value', ascending=False) # 'new_df' object from the stpe 4.5
RP_Descend_sort


# In[35]:


# Merge the DataFrames based on the 'ID' column
merged_Descend_df = pd.merge(RP_Descend_sort, gene_expression_data, on='GeneID')
merged_Descend_df 


# In[36]:


Des_sorted_geneexpressionData = merged_Descend_df.drop('RP value', axis=1)

# Convert 'GeneID' to the index column
Des_sorted_geneexpressionData.set_index('GeneID', inplace=True)
Des_sorted_geneexpressionData


# In[37]:


# Calculate the mean expression without considering the 'GeneID' column
    
# Create a copy of the data
Des_sorted_ExpData_Copy1 = Des_sorted_geneexpressionData.copy()

Des_sorted_ExpData_Copy1['Mean Expression'] = Des_sorted_ExpData_Copy1.mean(axis=1)
Des_sorted_ExpData_Copy1


# In[38]:


# Drop columns 'Mean Expression' from the DataFrame
Des_sorted_ExpData_Copy2 = Des_sorted_ExpData_Copy1.drop(['Mean Expression'], axis=1)

# Calculate the standard deviation by gene (row-wise)
Des_sorted_ExpData_Copy2['std_dev'] = Des_sorted_ExpData_Copy2.std(axis=1)
Des_sorted_ExpData_Copy2


# In[39]:


Mean_df = Des_sorted_ExpData_Copy1[['Mean Expression']]
stdv_df = Des_sorted_ExpData_Copy2[['std_dev']]
merged_dataFrame = pd.merge(RP_Descend_sort, Mean_df, on='GeneID')
merged_Des_RPmeanStdv = pd.merge(merged_dataFrame, stdv_df, on='GeneID')
merged_Des_RPmeanStdv


# In[40]:


# Select the column 1 using iloc
rp = merged_Des_RPmeanStdv.iloc[:, 1]
rp_array = np.array(rp)
rp_array


# In[41]:


def calculate_rank_product(gene_expression_data):
       
    # Convert DataFrame to array
    array_data = gene_expression_data.values
    index_values = gene_expression_data.index.values
    index_column_name = gene_expression_data.index.name
    mapped_array = np.hstack((index_values.reshape(-1, 1), array_data))

    # Extract the expression columns from the mapped array
    expression_columns = mapped_array[:, 1:].astype(float)

    # Rank the expression columns
    ranked_columns = np.apply_along_axis(lambda x: np.argsort(-x).argsort() + 1, axis=0, arr=expression_columns)

    # Create a new array by combining the index column and ranked expression columns
    ranked_array = np.hstack((mapped_array[:, 0].reshape(-1, 1), ranked_columns))

    # Get the number of columns in the array
    num_columns = ranked_array.shape[1]

    # Create a copy of the original array
    sorted_array = ranked_array.copy()

    sorted_array = sorted_array[:, 1:].astype(float)
    #Rprod = np.prod(sorted_array, axis=1)
    rank_product = np.prod(sorted_array, axis=1)

    return rank_product

def permutation_test(gene_expression_data, rank_product, num_permutations, seed=None):
    if seed is not None:
        np.random.seed(seed)
        random.seed(seed)

    permuted_rps = np.zeros((num_permutations, len(rank_product)))

    for i in range(num_permutations):
        permuted_data = gene_expression_data.copy()
        for column in gene_expression_data:
            permuted_data[column] = np.random.permutation(permuted_data[column])

        permuted_rank_product = calculate_rank_product(permuted_data)
        permuted_rps[i] = permuted_rank_product

    return permuted_rps


# In[42]:


rank_product = calculate_rank_product(Des_sorted_geneexpressionData) #Des_sorted_geneexpressionData
#rank_product = calculate_rank_product(gene_expression_data)
rank_product


# In[43]:


rank_productDF = pd.DataFrame(rank_product)
rank_productDF = rank_productDF.rename(columns={0: 'RP value'})
rank_productDF


# In[44]:


# Set the random seed
seed_value = 42

# Call permutation_test with the seed
permutedRP = permutation_test(Des_sorted_geneexpressionData, rank_product, num_permutations=1000, seed=seed_value)
#permutedRP = permutation_test(gene_expression_data, rank_product, num_permutations=1000, seed=seed_value)
permutedRP


# In[45]:


def calculate_p_values(rank_product, permuted_values, num_permutations):
    num_genes = len(rank_product)
    p_values = np.zeros(num_genes)

    for i in range(num_genes):
        # Calculate the observed rank product for the current gene
        observed_rank_product = rank_product[i]

        # Calculate the number of permuted rank products that are greater than or equal to the observed rank product
        num_greater = np.sum(permuted_values[:, i] >= observed_rank_product)

        # Calculate the p-value for the current gene by dividing the number of permuted rank products by the total number of permutations
        p_value = (num_greater + 1) / (num_permutations + 1)

        p_values[i] = p_value

    return p_values


# In[46]:


num_permutations=1000
p_values = calculate_p_values(rank_product, permutedRP, num_permutations)
p_values


# In[47]:


pval_df = pd.DataFrame(p_values)

# Rename the column
pval_df = pval_df.rename(columns={0: 'P value'})
pval_df


# In[48]:


# Attach RP value column from the rank_productDF to Des_sorted_geneexpressionData
Des_sorted_geneexpressionData['RP value'] = rank_productDF['RP value'].values

# Attach P value column from the pval_df to Des_sorted_geneexpressionData
Des_sorted_geneexpressionData['P value'] = pval_df['P value'].values

# Attach Mean Expression column from the merged_Des_RPmeanStdv to Des_sorted_geneexpressionData
Des_sorted_geneexpressionData['Mean Expression'] = merged_Des_RPmeanStdv['Mean Expression'].values

# Attach std_dev column from the merged_Des_RPmeanStdv to Des_sorted_geneexpressionData
Des_sorted_geneexpressionData['std_dev'] = merged_Des_RPmeanStdv['std_dev'].values

Descen_RP_Sorted_PvalTable  = Des_sorted_geneexpressionData[['RP value', 'P value', 'Mean Expression', 'std_dev']].copy()
Descen_RP_Sorted_PvalTable 


# In[49]:


# Set the P-value cutoff 
p_value_cutoff = 0.05

# Filter genes based on P-value cutoff
SG_bypval_df = Descen_RP_Sorted_PvalTable[Descen_RP_Sorted_PvalTable['P value'] <= p_value_cutoff]
#SG_bypval_df.to_csv('t1.csv')
SG_bypval_df


# In[50]:


# Sorting based mean expression levels
SG_bypval_Mean = SG_bypval_df.sort_values(by='Mean Expression', ascending=True)
SG_bypval_Mean


# In[51]:


# Calculate the mean expression level for each gene
mean_expression_levels = SG_bypval_df['Mean Expression']
mean_expression_levels = np.array(mean_expression_levels)

# Sort the genes based on their mean expression levels (descending order)
#sorted_genes = np.argsort(mean_expression_levels)[::-1]

# Sort the genes based on their mean expression levels in ascending order
sorted_genes = np.argsort(mean_expression_levels)

# Select the top genes with the highest consistent expression levels (e.g., top 10%)
top_genes_percentage = 0.1
num_top_genes = int(len(sorted_genes) * top_genes_percentage)
selected_high_consistent_genes = sorted_genes[:num_top_genes]

# Determine the minimum mean expression level among the selected high consistent genes
max_mean_exp_TopGenes_Descend = np.max(mean_expression_levels[selected_high_consistent_genes])

min_mean_exp_TopGenes_Descend = np.min(mean_expression_levels[selected_high_consistent_genes])

# Select the bottom genes with the lowest consistent expression levels (e.g., bottom 10%)
bottom_genes_percentage = 0.1
num_bottom_genes = int(len(sorted_genes) * bottom_genes_percentage)
selected_low_consistent_genes = sorted_genes[-num_bottom_genes:]

# Determine the maximum mean expression level among the selected low consistent genes
min_mean_exp_BottomGenes_Descend = np.min(mean_expression_levels[selected_low_consistent_genes])

# Determine the maximum mean expression level among the selected low consistent genes
max_mean_exp_BottomGenes_Descend = np.max(mean_expression_levels[selected_low_consistent_genes])

# Print the minimum and maximum mean expression levels
print("Maximum mean expression level in top genes:", max_mean_exp_TopGenes_Descend)
print("Minimum mean expression level in top genes:", min_mean_exp_TopGenes_Descend)
print("Maximum mean expression level in bottom genes:", max_mean_exp_BottomGenes_Descend)
print("Minimum mean expression level in bottom genes:", min_mean_exp_BottomGenes_Descend)


# In[52]:


# Set the P-value cutoff 
p_value_cutoff = 0.05

# Filter genes based on P-value cutoff
SG_bypval_df = Descen_RP_Sorted_PvalTable[Descen_RP_Sorted_PvalTable['P value'] <= p_value_cutoff]
#SG_bypval_df.to_csv('t1.csv')
SG_bypval_df


# In[53]:


min_mean_exp_BottomGenes_Descend


# In[54]:


# Categorize high consistent genes based on Mean Expression above the mean_Min
#low_consistent_genes = SG_bypval_df[(SG_bypval_df['Mean Expression'] >= 1.3460850563727338) & (SG_bypval_df['Mean Expression'] <= 3.7008423701846147)]
low_consistent_genes = SG_bypval_df[(SG_bypval_df['Mean Expression'] >= min_mean_exp_BottomGenes_Descend) & (SG_bypval_df['Mean Expression'] <= min_mean_exp_BottomGenes_Ascend)]
low_consistent_genes


# In[55]:


import pandas as pd
import mygene

def convert_ensembl_to_symbols(ensembl_ids):
    mg = mygene.MyGeneInfo()

    # Retrieve gene annotations for Ensembl IDs
    results = mg.querymany(ensembl_ids, scopes='ensembl.gene', fields='symbol', species='human')

    gene_symbols = []
    for result in results:
        if 'symbol' in result:
            ensembl_id = result['query']
            gene_symbol = result['symbol']
            gene_symbols.append((ensembl_id, gene_symbol))

    return gene_symbols

# Load the Ensemble gene IDs
ensembl_ids = low_consistent_genes.index
gene_symbols = convert_ensembl_to_symbols(ensembl_ids)

# Convert the list of tuples to a DataFrame
df = pd.DataFrame(gene_symbols, columns=['GeneID', 'Gene_Symbol'])
df


# In[56]:


# Convert the index column to a regular column
lcg_df= low_consistent_genes.reset_index()
lcg_df


# In[57]:


# Merge the DataFrames based on the first column (ID)
merged_df1 = pd.merge(df, lcg_df, on='GeneID')
merged_df1


# In[58]:


# Merge the DataFrames based on the first column (ID) and keep the unmerged rows
merged_df2 = pd.merge(df, lcg_df, on='GeneID', how='outer')
merged_df2


# In[59]:


Ascen_RP_Sorted_PvalTable1 = Ascen_RP_Sorted_PvalTable.rename(columns={'RP value': 'RP value(High)', 
                                                                      'P value': 'P value(High)'})

# Droping the Mean expression and standard deviation column from the data frame
Ascen_RP_Sorted_PvalTable1 = Ascen_RP_Sorted_PvalTable1.drop(['Mean Expression', 'std_dev'], axis=1)
Ascen_RP_Sorted_PvalTable2= Ascen_RP_Sorted_PvalTable1.reset_index()
Ascen_RP_Sorted_PvalTable2


# In[60]:


Descen_RP_Sorted_PvalTable1 = Descen_RP_Sorted_PvalTable.rename(columns={'RP value': 'RP value(Low)', 
                                                                        'P value': 'P value(Low)'})
Descen_RP_Sorted_PvalTable2= Descen_RP_Sorted_PvalTable1.reset_index()
Descen_RP_Sorted_PvalTable2


# In[61]:


FinalTable1 = Ascen_RP_Sorted_PvalTable2.merge(Descen_RP_Sorted_PvalTable2, on='GeneID')

# Set 'GeneID' as the index
FinalTable2 = FinalTable1.set_index('GeneID')
FinalTable2


# In[62]:


high_consistent_genes


# In[63]:


# Convert the index of the DataFrame into a list
HCG_List = high_consistent_genes.index.tolist()

# Extract rows using loc
HCG_DF = FinalTable2.loc[HCG_List]
HCG_DF = HCG_DF[['RP value(High)', 'RP value(Low)', 'P value(High)', 'P value(Low)', 'Mean Expression', 'std_dev']]
HCG_DF


# In[64]:


FinalTable1


# In[65]:


# Create a mask of rows whose Ensembl_ID is not in the list_of_ids
mask1 = ~FinalTable1['GeneID'].isin(HCG_List)

# Apply the mask to the DataFrame to get a DataFrame without the rows with those IDs
HCG_FilteredDF  = FinalTable1[mask1]
HCG_FilteredDF 


# In[66]:


# Convert the index of the DataFrame into a list
LCG_List = low_consistent_genes.index.tolist()

# Extract rows using iloc
LCG_DF = FinalTable2.loc[LCG_List]
LCG_DF = LCG_DF[['RP value(High)', 'RP value(Low)', 'P value(High)', 'P value(Low)', 'Mean Expression', 'std_dev']]
LCG_DF 


# In[67]:


# Create a mask of rows whose Ensembl_ID is not in the list_of_ids
mask2 = ~FinalTable1['GeneID'].isin(LCG_List)

# Apply the mask to the DataFrame to get a DataFrame without the rows with those IDs
InconsistentGenes  = HCG_FilteredDF[mask2]
# Set 'GeneID' as the index
InconsistentGenes = InconsistentGenes.set_index('GeneID')
InconsistentGenes


# In[68]:


# Convert the index of the DataFrame into a list
HCG_List = high_consistent_genes.index.tolist()

# Extract rows using iloc
HCG_gDF = gene_expression_data.loc[HCG_List]
HCG_gDF 


# In[69]:


# Convert the index of the DataFrame into a list
LCG_List = low_consistent_genes.index.tolist()

# Extract rows using iloc
LCG_gDF = gene_expression_data.loc[LCG_List]
LCG_gDF


# In[70]:


InconsistentGenes 


# In[71]:


InconsistentGenes2 = InconsistentGenes.reset_index()
# Convert the index of the DataFrame into a list
ICG_List = InconsistentGenes2['GeneID'].tolist()

# Extract rows using iloc
ICG_gDF = gene_expression_data.loc[ICG_List]
ICG_gDF 


# In[73]:


high_consistent_genes.to_csv('HCG_PE11_pVal0.05.csv')
low_consistent_genes.to_csv('LCG_PE11_pVal0.05.csv')
InconsistentGenes.to_csv('ICG_PE11_pVal0.05.csv') 


# In[74]:


import matplotlib.pyplot as plt
from matplotlib_venn import venn3

# sets
HCG_set = set(HCG_gDF.index)
LCG_set = set(LCG_gDF.index)
ICG_set = set(ICG_gDF.index)

# Create the Venn diagram
venn3([HCG_set, LCG_set, ICG_set], ('PE11_HCG', 'PE11_LCG', 'PE11_ICG'))

# Save the figure
plt.savefig('CGE_Overlap_PE11.png')

# Display the plot
plt.show()


# In[ ]:




