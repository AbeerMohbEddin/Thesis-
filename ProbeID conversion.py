#!/usr/bin/env python
# coding: utf-8

# In[12]:


import os
os.chdir("C:/Users/abeer/Desktop/HistamineSignature")


# In[14]:


import pandas as pd
from pybiomart import Server

# Step 1: Fetch the conversion mapping from BioMart
def fetch_mapping():
    server = Server(host='http://www.ensembl.org')
    dataset = server.marts['ENSEMBL_MART_ENSEMBL'].datasets['hsapiens_gene_ensembl']

    # Fetch the mapping
    attributes = ['affy_hg_u133_plus_2', 'ensembl_gene_id']  # Adjust the probe type as needed
    mapping = dataset.query(attributes=attributes)

    # Determine the correct column names
    probe_column = [col for col in mapping.columns if 'affy' in col.lower()][0]
    ensembl_column = [col for col in mapping.columns if 'gene stable id' in col.lower()][0]

    # Convert to dictionary for faster look-up
    mapping_dict = dict(zip(mapping[probe_column].values, mapping[ensembl_column].values))
    
    return mapping_dict


# Step 2: Process each file and add EnsemblID column
def process_file(filename, mapping):
    data = pd.read_csv(filename)

    # Map the ProbeID to EnsemblID
    data['EnsemblID'] = data['Probe ID'].map(lambda x: mapping.get(x, "Not_Found"))
    
    # Filter out rows with "Not_Found"
    filtered_data = data[data['EnsemblID'] != "Not_Found"]

    # Save the processed file
    new_filename = filename.replace(".csv", "_with_ensembl.csv")
    filtered_data.to_csv(new_filename, index=False)
    
    # Main Execution
mapping_dict = fetch_mapping()
files = ["TxPosSig.csv", "TxNegSig.csv", "ConNegSig.csv", "ConPosSig.csv"]
for file in files:
    process_file(file, mapping_dict)


# In[ ]:




