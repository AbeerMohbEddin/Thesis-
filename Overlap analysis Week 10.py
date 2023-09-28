#!/usr/bin/env python
# coding: utf-8

# In[5]:


import os
os.chdir("C:/Users/abeer/Desktop/week10")


# In[6]:


import pandas as pd

PE10UP = pd.read_csv("HCG_PE10.csv")
PE10DN = pd.read_csv("LCG_PE10.csv")
NP10UP = pd.read_csv("HCG_NP10.csv")
NP10DN = pd.read_csv("LCG_NP10.csv")
PE10NS = pd.read_csv("ICG_PE10.csv")
NP10NS = pd.read_csv("ICG_NP10.csv")
HisSig = pd.read_csv("Histamine_Signature.csv")


# In[7]:


# Create lists from dataframes
PE10UP_List = PE10UP['GeneID'].to_list()
PE10DN_List = PE10DN['GeneID'].to_list()
PE10NS_List = PE10NS['GeneID'].to_list()
NP10UP_List = NP10UP['GeneID'].to_list()
NP10DN_List = NP10DN['GeneID'].to_list()
NP10NS_List = NP10NS['GeneID'].to_list()
HisSig_List = HisSig['GeneID'].to_list()

# Convert lists to sets
PE10UP_Set = set(PE10UP_List)
PE10DN_Set = set(PE10DN_List)
PE10NS_Set = set(PE10NS_List)
NP10UP_Set = set(NP10UP_List)
NP10DN_Set = set(NP10DN_List)
NP10NS_Set = set(NP10NS_List)
HisSig_Set = set(HisSig_List)

# Calculate overlaps
overlap_PE10UP_HisSig = PE10UP_Set.intersection(HisSig_Set)
overlap_PE10DN_HisSig = PE10DN_Set.intersection(HisSig_Set)
overlap_PE10NS_HisSig = PE10NS_Set.intersection(HisSig_Set)
overlap_NP10UP_HisSig = NP10UP_Set.intersection(HisSig_Set)
overlap_NP10DN_HisSig = NP10DN_Set.intersection(HisSig_Set)
overlap_NP10NS_HisSig = NP10NS_Set.intersection(HisSig_Set)

# Print sizes of overlaps
print("Number of overlapping genes between PE10UP and HisSig:", len(overlap_PE10UP_HisSig))
print("Number of overlapping genes between PE10DN and HisSig:", len(overlap_PE10DN_HisSig))
print("Number of overlapping genes between PE10NS and HisSig:", len(overlap_PE10NS_HisSig))
print("Number of overlapping genes between NP10UP and HisSig:", len(overlap_NP10UP_HisSig))
print("Number of overlapping genes between NP10DN and HisSig:", len(overlap_NP10DN_HisSig))
print("Number of overlapping genes between NP10NS and HisSig:", len(overlap_NP10NS_HisSig))


# In[8]:


import matplotlib.pyplot as plt
from matplotlib_venn import venn2

def plot_venn(set1, set2, label1, label2, title):
    plt.figure()
    venn2([set1, set2], (label1, label2))
    plt.title(title)
    plt.show()

# Generate Venn diagrams
plot_venn(PE10UP_Set, HisSig_Set, "PE10UP", "HisSig", "Overlap between PE10UP and HisSig")


# In[9]:



# Generate Venn diagrams
plot_venn(PE10DN_Set, HisSig_Set, "PE10DN", "HisSig", "Overlap between PE10DN and HisSig")


# In[10]:


plot_venn(PE10NS_Set, HisSig_Set, "PE10NS", "HisSig", "Overlap between PE10NS and HisSig")


# In[11]:


plot_venn(NP10UP_Set, HisSig_Set, "NP10UP", "HisSig", "Overlap between NP10UP and HisSig")


# In[12]:


plot_venn(NP10DN_Set, HisSig_Set, "NP10DN", "HisSig", "Overlap between NP10DN and HisSig")


# In[27]:


plot_venn(NP10NS_Set, HisSig_Set, "NP10NS", "HisSig", "Overlap between NP10NS and HisSig")


# In[13]:



def extract_and_save_intersections(df1, df2, column_name, save_filename):
    # Get intersecting GeneIDs
    intersection = set(df1[column_name]).intersection(set(df2[column_name]))
    
    # Filter rows based on intersecting GeneIDs
    result_df = df1[df1[column_name].isin(intersection)]
    
    # Save to CSV
    result_df.to_csv(save_filename, index=False)

# Use the function to extract and save intersections with the HisSig data frame
extract_and_save_intersections(PE10UP, HisSig, 'GeneID', "PE10UP_HisSig_Intersection.csv")
extract_and_save_intersections(PE10DN, HisSig, 'GeneID', "PE10DN_HisSig_Intersection.csv")
extract_and_save_intersections(PE10NS, HisSig, 'GeneID', "PE10NS_HisSig_Intersection.csv")
extract_and_save_intersections(NP10UP, HisSig, 'GeneID', "NP10UP_HisSig_Intersection.csv")
extract_and_save_intersections(NP10DN, HisSig, 'GeneID', "NP10DN_HisSig_Intersection.csv")
extract_and_save_intersections(NP10NS, HisSig, 'GeneID', "NP10NS_HisSig_Intersection.csv")


# In[7]:


import pandas as pd
import matplotlib.pyplot as plt
from matplotlib_venn import venn2

# Load previously saved intersection CSVs
PE10UP_HisSig = pd.read_csv("PE10UP_HisSig_Intersection.csv")
PE10DN_HisSig = pd.read_csv("PE10DN_HisSig_Intersection.csv")
NP10UP_HisSig = pd.read_csv("NP10UP_HisSig_Intersection.csv")
NP10DN_HisSig = pd.read_csv("NP10DN_HisSig_Intersection.csv")

# Overlap PE8DN_HisSig with NP8DN_HisSig and PE8UP_HisSig with NP8UP_HisSig
#extract_and_save_intersections(PE8DN_HisSig, NP8DN_HisSig, 'GeneID', "PE8DN_NP8DN_HisSig_Intersection.csv")
#extract_and_save_intersections(PE8UP_HisSig, NP8UP_HisSig, 'GeneID', "PE8UP_NP8UP_HisSig_Intersection.csv")

# Isolate genes exclusive to PE intersections
def extract_and_save_exclusive(df1, df2, column_name, save_filename):
    exclusive = set(df1[column_name]).difference(set(df2[column_name]))
    result_df = df1[df1[column_name].isin(exclusive)]
    result_df.to_csv(save_filename, index=False)

# Get exclusive genes
extract_and_save_exclusive(PE10UP_HisSig, NP10UP_HisSig, 'GeneID', "Exclusive_PE10UP_HisSig.csv")
extract_and_save_exclusive(PE10DN_HisSig, NP10DN_HisSig, 'GeneID', "Exclusive_PE10DN_HisSig.csv")


def plot_venn(set1, set2, label1, label2, title, save_filename):  # Added save_filename as an argument
    plt.figure()
    venn2([set1, set2], (label1, label2))
    plt.title(title)
    plt.savefig(save_filename)  # Save the figure before displaying
    plt.show()
    
# Optional: Generate Venn diagrams for the new intersections
plot_venn(set(PE10UP_HisSig['GeneID']), set(NP10UP_HisSig['GeneID']), 
          "PE10UP_HisSig", "NP10UP_HisSig", 
          "Overlap between PE10UP_HisSig and NP10UP_HisSig", 
          "PE10UP_NP10UP_HisSig_Venn.png")  # Pass the filename for saving

plot_venn(set(PE10DN_HisSig['GeneID']), set(NP10DN_HisSig['GeneID']), 
          "PE10DN_HisSig", "NP10DN_HisSig", 
          "Overlap between PE10DN_HisSig and NP10DN_HisSig", 
          "PE10DN_NP10DN_HisSig_Venn.png")  # Pass the filename for saving


# In[ ]:




