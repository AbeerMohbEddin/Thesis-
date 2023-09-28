#!/usr/bin/env python
# coding: utf-8

# In[1]:


import os
os.chdir("C:/Users/abeer/Desktop/MScWD")


# The os module is imported to interact with the operating system. The os.chdir() function is used to set the current working directory to the specified path where your data files and/or R package are stored.

# # load libraries 

# In[2]:


import pandas as pd
import numpy as np
from scipy.stats import rankdata
from statsmodels.stats.multitest import multipletests
from itertools import permutations
import random


# The script is importing pandas, numpy, scipy.stats, statsmodels, and itertools which are Python libraries commonly used for data processing and statistical analysis.

# # step 3 load the data and data preprocessing for consistent genes analysis process

# # step 3.1: load the data

# In[3]:


gene_expression_data = pd.read_csv("PE_Wk13_MeV.csv", index_col=0)
gene_expression_data


# The pd.read_csv() function is used to load the gene expression data from a CSV file into a pandas DataFrame, with the first column being treated as the index column (GeneID).

# # step 3.2 convert data fram into array data 

# In[4]:


array = np.array(gene_expression_data)
array


# The data from the DataFrame is converted into a numpy array for easier processing.

# # step 3.3 join the gene ID column with array data

# In[5]:


index_values = gene_expression_data.index.values
index_columns_name = gene_expression_data.index.name
mapped_array = np.hstack((index_values.reshape(-1, 1), array))
mapped_array


# The GeneID (which was the index of the DataFrame) is combined with the array data, resulting in an array where the first column is the GeneID and the subsequent columns are gene expression data.

# # step 4.1 Rank the gene expression levels

# In[6]:


# Extract the expression columns from the mapped array
expression_columns = mapped_array[:, 1:].astype(float)

# Rank the expression columns
ranked_columns = np.apply_along_axis(lambda x: np.argsort(-x).argsort() + 1, axis=0, arr=expression_columns)

# Create a new array by combining the index column and ranked expression columns
ranked_array = np.hstack((mapped_array[:, 0].reshape(-1, 1), ranked_columns))

# Get the number of columns in the array
num_columns = ranked_array.shape[1]

ranked_array


# Here, the gene expression levels are ranked within each column (i.e., for each sample), and a new array is created with these ranks.

#  #                                 Identification of high consistent genes

# In[7]:


#step 4.2: convert teh ranked array data into dataframe


# In[8]:


ranked_df = pd.DataFrame(ranked_array)

# Rename the numeric column names
ranked_df = ranked_df.rename(columns={0: 'GeneID'})

ranked_df


# The ranked array is converted back into a DataFrame for easier processing and viewing.

# # step 4.3 calculate the rank product by confirming the values of column 1 (geneID)

# In[9]:


ranked_array2 = ranked_array[:, 1:].astype(float)
rank_product = np.prod(ranked_array2, axis=1)

rank_product


# The rank product is calculated for each gene (row). The rank product is a biologically motivated test for the detection of differentially expressed genes in replicated microarray experiments. It is a simple, non-parametric method that does not depend on specific assumptions about the distribution of the data.

# # step 4.4: convert the RP array data into dataframe

# In[10]:


rankProductdf = pd.DataFrame(rank_product)

#rename the numeric columns
rankProductdf = rankProductdf.rename (columns={0: "RP value"})
rankProductdf


# The resulting rank product values are converted into a DataFrame for easier processing and viewing.

# # step 4.5 : join the gene Id COLUMN with RP values 

# In[11]:


# Select column 1 using iloc
GeneID_column = ranked_df.iloc[:, 0]

# Select column 1 using iloc
RP_column = rankProductdf.iloc[:, 0]

# Create a new DataFrame by concatenating the two columns
new_df = pd.concat([GeneID_column, RP_column], axis=1)

new_df


# The GeneID column and the RP values column are combined into a new DataFrame.

# # step 5: producing sorted RP by ascending sorting

# In[12]:


# Sort the DataFrame by the "RP value" column in ascending order
RP_Ascend_sort = new_df.sort_values('RP value', ascending=True)

RP_Ascend_sort


# The DataFrame is sorted in ascending order based on the RP value. The genes with the smallest RP values will be at the top, indicating they have the most consistent high expression across the samples.

# # step 6: Producing sorted (Ascending)

# In[13]:


# Assuming your DataFrame is named "df"
gene_expression_data.reset_index(inplace=True)
gene_expression_data.rename(columns={'index': 'GeneID'}, inplace=True)


# In[14]:


# Merge the DataFrames based on the ID column
merged_Ascend_df = pd.merge(RP_Ascend_sort, gene_expression_data, on='GeneID')
merged_Ascend_df


# # Step 7: producing sorted (Ascending) gene expression data

# In[15]:


sorted_geneexpressionData = merged_Ascend_df.drop('RP value', axis=1)
sorted_geneexpressionData.set_index('GeneID', inplace=True)

sorted_geneexpressionData


# # step 8: producing sorted (Ascending) RP table with mean Expression and Standard Deviation

# In[16]:


merged_Ascend_df['Mean Expression'] = merged_Ascend_df.drop('RP value', axis=1).mean(axis=1)

merged_Ascend_df


# # step 8.2 Calculating SD 

# In[17]:


# Drop columns 'RP value' and 'Mean Expression' from the DataFrame
merged_Ascend_df_dropped = merged_Ascend_df.drop(['RP value', 'Mean Expression'], axis=1)

# Calculate the standard deviation by gene (row-wise)
merged_Ascend_df_dropped['std_dev'] = merged_Ascend_df_dropped.std(axis=1)

merged_Ascend_df_dropped


# # step 8.3

# In[18]:


# Select columns 1, 2, and 7 using iloc
columnRP_Mean = merged_Ascend_df.iloc[:, [0, 1, 14]]

# Select column 5 using iloc
column_StdDev = merged_Ascend_df_dropped.iloc[:, 13]

# Create a new DataFrame by concatenating the two columns
RP_SortedTable1 = pd.concat([columnRP_Mean, column_StdDev], axis=1)

RP_SortedTable1


# ### step 9

# In[19]:


# Select column 1 using iloc
rp = RP_SortedTable1.iloc[:, 1]

# Convert column 1 to a NumPy array
rp_array = np.array(rp)

rp_array


# # step 9.2: Generate the permuted RP values for p value calculation

# In[20]:


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

   Rprod = np.prod(sorted_array, axis=1)



   return Rprod





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


# In[21]:


rank_product = calculate_rank_product(sorted_geneexpressionData)
rank_product


# In[22]:


seed_value = 42
permutedRP = permutation_test(sorted_geneexpressionData, rank_product, num_permutations=1000, seed=seed_value)
permutedRP


# # step 9.3

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


num_permutations = 1000
p_values = calculate_p_values(rank_product, permutedRP, num_permutations)
p_values


# # step 9.4

# In[25]:


pval_df = pd.DataFrame(p_values)
# Rename the column
pval_df = pval_df.rename(columns={0: 'P value'})
pval_df


# # step 10

# In[26]:


# Perform the merge based on the index column
RP_Sorted_PvalTable = pd.merge(RP_SortedTable1, pval_df, left_index=True, right_index=True)

# Convert 'Column1' to the index column
RP_Sorted_PvalTable.set_index('GeneID', inplace=True)

RP_Sorted_PvalTable


# # step 11

# In[27]:


# Set the P-value cutoff
p_value_cutoff = 0.001

# Filter genes based on P-value cutoff
SG_bypval_df = RP_Sorted_PvalTable[RP_Sorted_PvalTable['P value'] <= p_value_cutoff]

SG_bypval_df.to_csv('HCG.csv', index=False)


# In[28]:


# Calculate the mean expression level for each gene
mean_expression_levels = SG_bypval_df["Mean Expression"]
mean_expression_levels = np.array(mean_expression_levels)

# Sort the genes based on their mean expression levels (descending order)
sorted_genes = np.argsort(mean_expression_levels)[::-1]

# Select the top genes with the highest consistent expression levels (e.g., top 10%)
top_genes_percentage = 0.1
num_top_genes = int(len(sorted_genes) * top_genes_percentage)
selected_high_consistent_genes = sorted_genes[:num_top_genes]

# Determine the minimum mean expression level among the selected high consistent genes
max_mean_expression_level_TopGenes = np.max(mean_expression_levels[selected_high_consistent_genes])
min_mean_expression_level_TopGenes = np.min(mean_expression_levels[selected_high_consistent_genes])

# Select the bottom genes with the lowest consistent expression levels (e.g., bottom 10%)
bottom_genes_percentage = 0.1
num_bottom_genes = int(len(sorted_genes) * bottom_genes_percentage)
selected_low_consistent_genes = sorted_genes[-num_bottom_genes:]

# Determine the maximum mean expression level among the selected low consistent genes
min_mean_expression_level_BottomGenes = np.min(mean_expression_levels[selected_low_consistent_genes])
max_mean_expression_level_BottomGenes = np.max(mean_expression_levels[selected_low_consistent_genes])

# Print the minimum and maximum mean expression levels
print("Maximum mean expression level in top genes:", max_mean_expression_level_TopGenes)
print("Minimum mean expression level in top genes:", min_mean_expression_level_TopGenes)
print("Maximum mean expression level in bottom genes:", max_mean_expression_level_BottomGenes)
print("Minimum mean expression level in bottom genes:", min_mean_expression_level_BottomGenes)


# In[29]:


# Set the P-value cutoff and calculate the 75th percentile of Mean Expression
p_value_cutoff = 0.001
mean_Min = 3.4864588668333325

# Filter genes based on P-value cutoff
filtered_df = SG_bypval_df[SG_bypval_df['P value'] <= p_value_cutoff]

# Categorize high consistent genes based on Mean Expression above the mean_Min
high_consistent_genes = filtered_df[filtered_df['Mean Expression'] > mean_Min]

high_consistent_genes


# In[30]:


get_ipython().system('pip install mygene')


# In[31]:


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

# Load the Ensembl gene IDs
ensembl_ids = high_consistent_genes.index
# Convert Ensembl IDs to symbols
gene_symbols = convert_ensembl_to_symbols(ensembl_ids)

# Convert the list of tuples to a DataFrame
df = pd.DataFrame(gene_symbols, columns=['GeneID', 'Symbol'])
df


# In[32]:


hcg_df = high_consistent_genes.reset_index()
hcg_df


# In[33]:


merged_df1 = pd.merge(df, hcg_df, on='GeneID')
merged_df1


# In[34]:


#Merge the DataFrames based on the first column (ID) and keep the unmerged rows
merged_df2 = pd.merge(df, hcg_df, on='GeneID', how='outer')
merged_df2


# In[35]:


merged_df2.to_csv('WEEK11', index=False)


# In[ ]:




