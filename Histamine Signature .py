#!/usr/bin/env python
# coding: utf-8

# In[27]:


import os
os.chdir("C:/Users/abeer/Desktop/HistamineSignature")


# In[28]:


import pandas as pd

# Load histamine regulated gene files
TxNeg = pd.read_csv("TxNegSig_with_ensembl.csv")
TxPos = pd.read_csv("TxPosSig_with_ensembl.csv")
ConNeg = pd.read_csv("ConNegSig_with_ensembl.csv")
ConPos = pd.read_csv("ConPosSig_with_ensembl.csv")


# In[29]:


import pandas as pd
import matplotlib.pyplot as plt
from matplotlib_venn import venn2

# Load histamine regulated gene files from Excel format
TxNeg = pd.read_csv("TxNegSig_with_ensembl.csv")
ConNeg = pd.read_csv("ConNegSig_with_ensembl.csv")

# Assuming "Gene_Symbol" is the column for gene names
txneg_genes = set(TxNeg['EnsemblID'])
conneg_genes = set(ConNeg['EnsemblID'])

# Create a Venn diagram
venn2([txneg_genes, conneg_genes], ('TxNeg', 'ConNeg'))

plt.show()


# In[30]:


import pandas as pd
import matplotlib.pyplot as plt
from matplotlib_venn import venn2

# Load histamine regulated gene files from Excel format
TxPos = pd.read_csv("TxPosSig_with_ensembl.csv")
ConPos = pd.read_csv("ConPosSig_with_ensembl.csv")

# Assuming "Gene_Symbol" is the column for gene names
txpos_genes = set(TxPos['EnsemblID'])
conpos_genes = set(ConPos['EnsemblID'])

# Create a Venn diagram
venn2([txpos_genes, conpos_genes], ('TxPos', 'ConPos'))

plt.show()


# In[31]:


# Identify the genes unique to TxPos (that are not in ConNeg)
unique_txneg_genes = txneg_genes.difference(conneg_genes)

# Filter the TxNeg DataFrame for these unique genes
unique_txneg_df = TxNeg[TxNeg['EnsemblID'].isin(unique_txneg_genes)]


# In[32]:


print("Number of unique genes specific to TxNeg:", len(unique_txneg_genes))


# In[33]:



# Save the unique genes to a new file
unique_txneg_df.to_csv("Unique_TxNeg.csv", index=False)


# In[34]:


# Identify the genes unique to TxPos (that are not in ConNeg)
unique_txpos_genes = txpos_genes.difference(conpos_genes)

# Filter the TxNeg DataFrame for these unique genes
unique_txpos_df = TxPos[TxPos['EnsemblID'].isin(unique_txpos_genes)]


# In[35]:


print("Number of unique genes specific to TxNeg:", len(unique_txpos_genes))


# In[36]:


# Save the unique genes to a new file
unique_txpos_df.to_csv("Unique_TxPos.csv", index=False)


# In[ ]:




