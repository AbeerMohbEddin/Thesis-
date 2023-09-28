#!/usr/bin/env python
# coding: utf-8

# In[1]:


import os
os.chdir("C:/Users/abeer/Desktop/week13")


# In[18]:


import os

directory_path = 'C:/Users/abeer/Desktop/week13'  # Replace with your directory path
files = os.listdir(directory_path)

for file in files:
    print(file)


# In[19]:


import pandas as pd

files_to_convert = [
    "HCG_NP13.xlsx",
    "HCG_PE13.xlsx",
    "ICG_NP13.xlsx",
    "ICG_PE13.xlsx",
    "LCG_NP13.xlsx",
    "LCG_PE13.xlsx"
]

for file in files_to_convert:
    # Read the Excel file
    data = pd.read_excel(file)
    
    # Convert to CSV
    csv_name = file.replace(".xlsx", ".csv")
    data.to_csv(csv_name, index=False)


# In[20]:


import pandas as pd

PE13UP = pd.read_csv("HCG_PE13.csv")
PE13DN = pd.read_csv("LCG_PE13.csv")
NP13UP = pd.read_csv("HCG_NP13.csv")
NP13DN = pd.read_csv("LCG_NP13.csv")
PE13NS = pd.read_csv("ICG_PE13.csv")
NP13NS = pd.read_csv("ICG_NP13.csv")
HisSig = pd.read_csv("Histamine_Signature.csv")


# In[21]:


# Create lists from dataframes
PE13UP_List = PE13UP['GeneID'].to_list()
PE13DN_List = PE13DN['GeneID'].to_list()
PE13NS_List = PE13NS['GeneID'].to_list()
NP13UP_List = NP13UP['GeneID'].to_list()
NP13DN_List = NP13DN['GeneID'].to_list()
NP13NS_List = NP13NS['GeneID'].to_list()
HisSig_List = HisSig['GeneID'].to_list()

# Convert lists to sets
PE13UP_Set = set(PE13UP_List)
PE13DN_Set = set(PE13DN_List)
PE13NS_Set = set(PE13NS_List)
NP13UP_Set = set(NP13UP_List)
NP13DN_Set = set(NP13DN_List)
NP13NS_Set = set(NP13NS_List)
HisSig_Set = set(HisSig_List)

# Calculate overlaps
overlap_PE13UP_HisSig = PE13UP_Set.intersection(HisSig_Set)
overlap_PE13DN_HisSig = PE13DN_Set.intersection(HisSig_Set)
overlap_PE13NS_HisSig = PE13NS_Set.intersection(HisSig_Set)
overlap_NP13UP_HisSig = NP13UP_Set.intersection(HisSig_Set)
overlap_NP13DN_HisSig = NP13DN_Set.intersection(HisSig_Set)
overlap_NP13NS_HisSig = NP13NS_Set.intersection(HisSig_Set)

# Print sizes of overlaps
print("Number of overlapping genes between PE13UP and HisSig:", len(overlap_PE13UP_HisSig))
print("Number of overlapping genes between PE13DN and HisSig:", len(overlap_PE13DN_HisSig))
print("Number of overlapping genes between PE13NS and HisSig:", len(overlap_PE13NS_HisSig))
print("Number of overlapping genes between NP13UP and HisSig:", len(overlap_NP13UP_HisSig))
print("Number of overlapping genes between NP13DN and HisSig:", len(overlap_NP13DN_HisSig))
print("Number of overlapping genes between NP13NS and HisSig:", len(overlap_NP13NS_HisSig))


# In[22]:


import matplotlib.pyplot as plt
from matplotlib_venn import venn2

def plot_venn(set1, set2, label1, label2, title):
    plt.figure()
    venn2([set1, set2], (label1, label2))
    plt.title(title)
    plt.show()

# Generate Venn diagrams
plot_venn(PE13UP_Set, HisSig_Set, "PE13UP", "HisSig", "Overlap between PE13UP and HisSig")


# In[23]:



# Generate Venn diagrams
plot_venn(PE13DN_Set, HisSig_Set, "PE13DN", "HisSig", "Overlap between PE13DN and HisSig")


# In[24]:


plot_venn(PE13NS_Set, HisSig_Set, "PE13NS", "HisSig", "Overlap between PE13NS and HisSig")


# In[25]:


plot_venn(NP13UP_Set, HisSig_Set, "NP13UP", "HisSig", "Overlap between NP13UP and HisSig")


# In[26]:


plot_venn(NP13DN_Set, HisSig_Set, "NP13DN", "HisSig", "Overlap between NP13DN and HisSig")


# In[30]:


plot_venn(NP13NS_Set, HisSig_Set, "NP13NS", "HisSig", "Overlap between NP13NS and HisSig")


# In[29]:



def extract_and_save_intersections(df1, df2, column_name, save_filename):
    # Get intersecting GeneIDs
    intersection = set(df1[column_name]).intersection(set(df2[column_name]))
    
    # Filter rows based on intersecting GeneIDs
    result_df = df1[df1[column_name].isin(intersection)]
    
    # Save to CSV
    result_df.to_csv(save_filename, index=False)

# Use the function to extract and save intersections with the HisSig data frame
extract_and_save_intersections(PE13UP, HisSig, 'GeneID', "PE13UP_HisSig_Intersection.csv")
extract_and_save_intersections(PE13DN, HisSig, 'GeneID', "PE13DN_HisSig_Intersection.csv")
extract_and_save_intersections(PE13NS, HisSig, 'GeneID', "PE13NS_HisSig_Intersection.csv")
extract_and_save_intersections(NP13UP, HisSig, 'GeneID', "NP13UP_HisSig_Intersection.csv")
extract_and_save_intersections(NP13DN, HisSig, 'GeneID', "NP13DN_HisSig_Intersection.csv")
extract_and_save_intersections(NP13NS, HisSig, 'GeneID', "NP13NS_HisSig_Intersection.csv")


# In[2]:


import pandas as pd
import matplotlib.pyplot as plt
from matplotlib_venn import venn2

# Load previously saved intersection CSVs
PE13UP_HisSig = pd.read_csv("PE13UP_HisSig_Intersection.csv")
PE13DN_HisSig = pd.read_csv("PE13DN_HisSig_Intersection.csv")
NP13UP_HisSig = pd.read_csv("NP13UP_HisSig_Intersection.csv")
NP13DN_HisSig = pd.read_csv("NP13DN_HisSig_Intersection.csv")

# Overlap PE8DN_HisSig with NP8DN_HisSig and PE8UP_HisSig with NP8UP_HisSig
#extract_and_save_intersections(PE8DN_HisSig, NP8DN_HisSig, 'GeneID', "PE8DN_NP8DN_HisSig_Intersection.csv")
#extract_and_save_intersections(PE8UP_HisSig, NP8UP_HisSig, 'GeneID', "PE8UP_NP8UP_HisSig_Intersection.csv")

# Isolate genes exclusive to PE intersections
def extract_and_save_exclusive(df1, df2, column_name, save_filename):
    exclusive = set(df1[column_name]).difference(set(df2[column_name]))
    result_df = df1[df1[column_name].isin(exclusive)]
    result_df.to_csv(save_filename, index=False)

# Get exclusive genes
extract_and_save_exclusive(PE13UP_HisSig, NP13UP_HisSig, 'GeneID', "Exclusive_PE13UP_HisSig.csv")
extract_and_save_exclusive(PE13DN_HisSig, NP13DN_HisSig, 'GeneID', "Exclusive_PE13DN_HisSig.csv")


def plot_venn(set1, set2, label1, label2, title):
    plt.figure()
    venn2([set1, set2], (label1, label2))
    plt.title(title)
    plt.show()
# Optional: Generate Venn diagrams for the new intersections
plot_venn(set(PE13UP_HisSig['GeneID']), set(NP13UP_HisSig['GeneID']), "PE13UP_HisSig", "NP13UP_HisSig", "Overlap between PE13UP_HisSig and NP13UP_HisSig")
plot_venn(set(PE13DN_HisSig['GeneID']), set(NP13DN_HisSig['GeneID']), "PE13DN_HisSig", "NP13DN_HisSig", "Overlap between PE13DN_HisSig and NP13DN_HisSig")


# In[ ]:




