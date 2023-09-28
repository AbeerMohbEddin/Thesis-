#!/usr/bin/env python
# coding: utf-8

# In[28]:


import os
os.chdir("C:/Users/abeer/Desktop/track")


# In[16]:


import pandas as pd

def load_genes_from_csv(filename):
    df = pd.read_csv(filename)
    return df['GeneID'].tolist()

def load_data_for_weeks():
    weeks = list(range(8, 17))
    data = {}
    for week in weeks:
        data[f"NP{week}DN"] = load_genes_from_csv(f'NP{week}DN_HisSig_Intersection.csv')
        data[f"NP{week}UP"] = load_genes_from_csv(f'NP{week}UP_HisSig_Intersection.csv')
        data[f"NP{week}NS"] = load_genes_from_csv(f'NP{week}NS_HisSig_Intersection.csv')
        data[f"PE{week}DN"] = load_genes_from_csv(f'PE{week}DN_HisSig_Intersection.csv')
        data[f"PE{week}UP"] = load_genes_from_csv(f'PE{week}UP_HisSig_Intersection.csv')
        data[f"PE{week}NS"] = load_genes_from_csv(f'PE{week}NS_HisSig_Intersection.csv')
    return data

Week8sig = load_genes_from_csv('Week8sig.csv')
all_data = load_data_for_weeks()


# In[17]:


def track_gene_across_weeks(gene, data):
    result = {}
    for category, genes in data.items():
        result[category] = gene in genes
    return result

tracking_results = {gene: track_gene_across_weeks(gene, all_data) for gene in Week8sig}


# In[18]:


tracking_results


# In[19]:


for gene, categories in tracking_results.items():
    print(f"Gene: {gene}")
    for category, present in categories.items():
        if present:
            print(f"    Found in {category}")
    print("-" * 50)


# In[20]:


import pandas as pd

df_results = pd.DataFrame(tracking_results).transpose()
print(df_results)


# In[21]:


import seaborn as sns
import matplotlib.pyplot as plt

plt.figure(figsize=(12, 8))
sns.heatmap(df_results, cmap="coolwarm", cbar=False)
plt.title("Gene Tracking Across Weeks")
plt.ylabel("GeneID")
plt.xlabel("Category")
plt.show()


# In[22]:


df_results.to_csv('tracking_results.csv')


# In[11]:


import pandas as pd

def load_genes_from_csv(filename):
    df = pd.read_csv(filename)
    return df['GeneID'].tolist()

def load_data_for_weeks():
    weeks = list(range(9, 17))
    data = {}
    for week in weeks:
        data[f"NP{week}DN"] = load_genes_from_csv(f'NP{week}DN_HisSig_Intersection.csv')
        data[f"NP{week}UP"] = load_genes_from_csv(f'NP{week}UP_HisSig_Intersection.csv')
        data[f"NP{week}NS"] = load_genes_from_csv(f'NP{week}NS_HisSig_Intersection.csv')
        data[f"PE{week}DN"] = load_genes_from_csv(f'PE{week}DN_HisSig_Intersection.csv')
        data[f"PE{week}UP"] = load_genes_from_csv(f'PE{week}UP_HisSig_Intersection.csv')
        data[f"PE{week}NS"] = load_genes_from_csv(f'PE{week}NS_HisSig_Intersection.csv')
    return data

Week9sig = load_genes_from_csv('Week9sig.csv')
all_data = load_data_for_weeks()


# In[12]:


def track_gene_across_weeks(gene, data):
    result = {}
    for category, genes in data.items():
        result[category] = gene in genes
    return result

tracking_results = {gene: track_gene_across_weeks(gene, all_data) for gene in Week9sig}


# In[13]:


tracking_results


# In[14]:


for gene, categories in tracking_results.items():
    print(f"Gene: {gene}")
    for category, present in categories.items():
        if present:
            print(f"    Found in {category}")
    print("-" * 50)


# In[3]:


import pandas as pd

def load_genes_from_csv(filename):
    df = pd.read_csv(filename)
    return df['GeneID'].tolist()

def load_data_for_weeks():
    weeks = list(range(8, 17))
    data = {}
    for week in weeks:
        data[f"PE{week}DN"] = load_genes_from_csv(f'PE{week}DN_HisSig_Intersection.csv')
        data[f"PE{week}UP"] = load_genes_from_csv(f'PE{week}UP_HisSig_Intersection.csv')
        data[f"PE{week}NS"] = load_genes_from_csv(f'PE{week}NS_HisSig_Intersection.csv')
    return data

Week8sig = load_genes_from_csv('Week8sig.csv')
all_data = load_data_for_weeks()


# In[4]:


def track_gene_across_weeks(gene, data):
    result = {}
    for category, genes in data.items():
        result[category] = gene in genes
    return result

tracking_results = {gene: track_gene_across_weeks(gene, all_data) for gene in Week8sig}


# In[5]:


tracking_results


# In[6]:


for gene, categories in tracking_results.items():
    print(f"Gene: {gene}")
    for category, present in categories.items():
        if present:
            print(f"    Found in {category}")
    print("-" * 50)


# In[58]:


import pandas as pd

def load_genes_from_csv(filename):
    df = pd.read_csv(filename)
    return df['GeneID'].tolist()

def load_data_for_weeks():
    weeks = list(range(10, 17))
    data = {}
    for week in weeks:
        data[f"NP{week}DN"] = load_genes_from_csv(f'NP{week}DN_HisSig_Intersection.csv')
        data[f"NP{week}UP"] = load_genes_from_csv(f'NP{week}UP_HisSig_Intersection.csv')
        data[f"NP{week}NS"] = load_genes_from_csv(f'NP{week}NS_HisSig_Intersection.csv')
        data[f"PE{week}DN"] = load_genes_from_csv(f'PE{week}DN_HisSig_Intersection.csv')
        data[f"PE{week}UP"] = load_genes_from_csv(f'PE{week}UP_HisSig_Intersection.csv')
        data[f"PE{week}NS"] = load_genes_from_csv(f'PE{week}NS_HisSig_Intersection.csv')
    return data

Week10sig = load_genes_from_csv('Week10sig.csv')
all_data = load_data_for_weeks()


# In[39]:


tracking_results


# In[40]:


for gene, categories in tracking_results.items():
    print(f"Gene: {gene}")
    for category, present in categories.items():
        if present:
            print(f"    Found in {category}")
    print("-" * 50)


# In[7]:


import pandas as pd

def load_genes_from_csv(filename):
    df = pd.read_csv(filename)
    return df['GeneID'].tolist()

def load_data_for_weeks():
    weeks = list(range(11, 17))
    data = {}
    for week in weeks:
        data[f"NP{week}DN"] = load_genes_from_csv(f'NP{week}DN_HisSig_Intersection.csv')
        data[f"NP{week}UP"] = load_genes_from_csv(f'NP{week}UP_HisSig_Intersection.csv')
        data[f"NP{week}NS"] = load_genes_from_csv(f'NP{week}NS_HisSig_Intersection.csv')
        data[f"PE{week}DN"] = load_genes_from_csv(f'PE{week}DN_HisSig_Intersection.csv')
        data[f"PE{week}UP"] = load_genes_from_csv(f'PE{week}UP_HisSig_Intersection.csv')
        data[f"PE{week}NS"] = load_genes_from_csv(f'PE{week}NS_HisSig_Intersection.csv')
    return data

Week11sig = load_genes_from_csv('Week11sig.csv')
all_data = load_data_for_weeks()


# In[8]:


def track_gene_across_weeks(gene, data):
    result = {}
    for category, genes in data.items():
        result[category] = gene in genes
    return result

tracking_results = {gene: track_gene_across_weeks(gene, all_data) for gene in Week11sig}


# In[9]:


tracking_results


# In[10]:


for gene, categories in tracking_results.items():
    print(f"Gene: {gene}")
    for category, present in categories.items():
        if present:
            print(f"    Found in {category}")
    print("-" * 50)


# In[34]:


import pandas as pd

def load_genes_from_csv(filename):
    df = pd.read_csv(filename)
    return df['GeneID'].tolist()

def load_data_for_weeks():
    weeks = list(range(10, 17))
    data = {}
    for week in weeks:
        data[f"NP{week}DN"] = load_genes_from_csv(f'NP{week}DN_HisSig_Intersection.csv')
        data[f"NP{week}UP"] = load_genes_from_csv(f'NP{week}UP_HisSig_Intersection.csv')
        data[f"NP{week}NS"] = load_genes_from_csv(f'NP{week}NS_HisSig_Intersection.csv')
    return data

Week11sig = load_genes_from_csv('week11sig.csv')
all_data = load_data_for_weeks()


# In[35]:


def track_gene_across_weeks(gene, data):
    result = {}
    for category, genes in data.items():
        result[category] = gene in genes
    return result

tracking_results = {gene: track_gene_across_weeks(gene, all_data) for gene in Week11sig}


# In[36]:


for gene, categories in tracking_results.items():
    print(f"Gene: {gene}")
    for category, present in categories.items():
        if present:
            print(f"    Found in {category}")
    print("-" * 50)


# In[12]:


import pandas as pd

def load_genes_from_csv(filename):
    df = pd.read_csv(filename)
    return df['GeneID'].tolist()

def load_data_for_weeks():
    weeks = list(range(12, 17))
    data = {}
    for week in weeks:
        data[f"NP{week}DN"] = load_genes_from_csv(f'NP{week}DN_HisSig_Intersection.csv')
        data[f"NP{week}UP"] = load_genes_from_csv(f'NP{week}UP_HisSig_Intersection.csv')
        data[f"NP{week}NS"] = load_genes_from_csv(f'NP{week}NS_HisSig_Intersection.csv')
        data[f"PE{week}DN"] = load_genes_from_csv(f'PE{week}DN_HisSig_Intersection.csv')
        data[f"PE{week}UP"] = load_genes_from_csv(f'PE{week}UP_HisSig_Intersection.csv')
        data[f"PE{week}NS"] = load_genes_from_csv(f'PE{week}NS_HisSig_Intersection.csv')
    return data

Week12sig = load_genes_from_csv('week12sig.csv')
all_data = load_data_for_weeks()


# In[13]:


def track_gene_across_weeks(gene, data):
    result = {}
    for category, genes in data.items():
        result[category] = gene in genes
    return result

tracking_results = {gene: track_gene_across_weeks(gene, all_data) for gene in Week12sig}


# In[14]:


tracking_results


# In[15]:


for gene, categories in tracking_results.items():
    print(f"Gene: {gene}")
    for category, present in categories.items():
        if present:
            print(f"    Found in {category}")
    print("-" * 50)


# In[16]:


import pandas as pd

def load_genes_from_csv(filename):
    df = pd.read_csv(filename)
    return df['GeneID'].tolist()

def load_data_for_weeks():
    weeks = list(range(13, 17))
    data = {}
    for week in weeks:
        data[f"NP{week}DN"] = load_genes_from_csv(f'NP{week}DN_HisSig_Intersection.csv')
        data[f"NP{week}UP"] = load_genes_from_csv(f'NP{week}UP_HisSig_Intersection.csv')
        data[f"NP{week}NS"] = load_genes_from_csv(f'NP{week}NS_HisSig_Intersection.csv')
        data[f"PE{week}DN"] = load_genes_from_csv(f'PE{week}DN_HisSig_Intersection.csv')
        data[f"PE{week}UP"] = load_genes_from_csv(f'PE{week}UP_HisSig_Intersection.csv')
        data[f"PE{week}NS"] = load_genes_from_csv(f'PE{week}NS_HisSig_Intersection.csv')
    return data

Week13sig = load_genes_from_csv('week13sig.csv')
all_data = load_data_for_weeks()


# In[17]:


def track_gene_across_weeks(gene, data):
    result = {}
    for category, genes in data.items():
        result[category] = gene in genes
    return result

tracking_results = {gene: track_gene_across_weeks(gene, all_data) for gene in Week13sig}
tracking_results


# In[18]:


for gene, categories in tracking_results.items():
    print(f"Gene: {gene}")
    for category, present in categories.items():
        if present:
            print(f"    Found in {category}")
    print("-" * 50)


# In[22]:


import pandas as pd

def load_genes_from_csv(filename):
    df = pd.read_csv(filename)
    return df['GeneID'].tolist()

def load_data_for_weeks():
    weeks = list(range(14, 17))
    data = {}
    for week in weeks:
        data[f"NP{week}DN"] = load_genes_from_csv(f'NP{week}DN_HisSig_Intersection.csv')
        data[f"NP{week}UP"] = load_genes_from_csv(f'NP{week}UP_HisSig_Intersection.csv')
        data[f"NP{week}NS"] = load_genes_from_csv(f'NP{week}NS_HisSig_Intersection.csv')
        data[f"PE{week}DN"] = load_genes_from_csv(f'PE{week}DN_HisSig_Intersection.csv')
        data[f"PE{week}UP"] = load_genes_from_csv(f'PE{week}UP_HisSig_Intersection.csv')
        data[f"PE{week}NS"] = load_genes_from_csv(f'PE{week}NS_HisSig_Intersection.csv')
    return data

Week14sig = load_genes_from_csv('week14sig.csv')
all_data = load_data_for_weeks()


# In[23]:


def track_gene_across_weeks(gene, data):
    result = {}
    for category, genes in data.items():
        result[category] = gene in genes
    return result

tracking_results = {gene: track_gene_across_weeks(gene, all_data) for gene in Week14sig}
tracking_results


# In[24]:


for gene, categories in tracking_results.items():
    print(f"Gene: {gene}")
    for category, present in categories.items():
        if present:
            print(f"    Found in {category}")
    print("-" * 50)


# In[33]:


import pandas as pd

def load_genes_from_csv(filename):
    df = pd.read_csv(filename)
    return df['GeneID'].tolist()

def load_data_for_weeks():
    weeks = list(range(15, 17))
    data = {}
    for week in weeks:
        data[f"NP{week}DN"] = load_genes_from_csv(f'NP{week}DN_HisSig_Intersection.csv')
        data[f"NP{week}UP"] = load_genes_from_csv(f'NP{week}UP_HisSig_Intersection.csv')
        data[f"NP{week}NS"] = load_genes_from_csv(f'NP{week}NS_HisSig_Intersection.csv')
        data[f"PE{week}DN"] = load_genes_from_csv(f'PE{week}DN_HisSig_Intersection.csv')
        data[f"PE{week}UP"] = load_genes_from_csv(f'PE{week}UP_HisSig_Intersection.csv')
        data[f"PE{week}NS"] = load_genes_from_csv(f'PE{week}NS_HisSig_Intersection.csv')
    return data

Week15sig = load_genes_from_csv('week15sig.csv')
all_data = load_data_for_weeks()


# In[34]:


def track_gene_across_weeks(gene, data):
    result = {}
    for category, genes in data.items():
        result[category] = gene in genes
    return result

tracking_results = {gene: track_gene_across_weeks(gene, all_data) for gene in Week15sig}
tracking_results


# In[35]:


for gene, categories in tracking_results.items():
    print(f"Gene: {gene}")
    for category, present in categories.items():
        if present:
            print(f"    Found in {category}")
    print("-" * 50)


# In[36]:


import pandas as pd

def load_genes_from_csv(filename):
    df = pd.read_csv(filename)
    return df['GeneID'].tolist()

def load_data_for_weeks():
    weeks = list(range(16, 17))
    data = {}
    for week in weeks:
        data[f"NP{week}DN"] = load_genes_from_csv(f'NP{week}DN_HisSig_Intersection.csv')
        data[f"NP{week}UP"] = load_genes_from_csv(f'NP{week}UP_HisSig_Intersection.csv')
        data[f"NP{week}NS"] = load_genes_from_csv(f'NP{week}NS_HisSig_Intersection.csv')
        data[f"PE{week}DN"] = load_genes_from_csv(f'PE{week}DN_HisSig_Intersection.csv')
        data[f"PE{week}UP"] = load_genes_from_csv(f'PE{week}UP_HisSig_Intersection.csv')
        data[f"PE{week}NS"] = load_genes_from_csv(f'PE{week}NS_HisSig_Intersection.csv')
    return data

Week16sig = load_genes_from_csv('exclusive_genes_week16sig.csv')
all_data = load_data_for_weeks()


# In[37]:


def track_gene_across_weeks(gene, data):
    result = {}
    for category, genes in data.items():
        result[category] = gene in genes
    return result

tracking_results = {gene: track_gene_across_weeks(gene, all_data) for gene in Week16sig}
tracking_results


# In[38]:


for gene, categories in tracking_results.items():
    print(f"Gene: {gene}")
    for category, present in categories.items():
        if present:
            print(f"    Found in {category}")
    print("-" * 50)


# In[ ]:




