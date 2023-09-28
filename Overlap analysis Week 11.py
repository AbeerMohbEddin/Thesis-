#!/usr/bin/env python
# coding: utf-8

# In[1]:


import os
os.chdir("C:/Users/abeer/Desktop/week11")


# In[41]:


import pandas as pd

PE11UP = pd.read_csv("HCG_PE11.csv")
PE11DN = pd.read_csv("LCG_PE11.csv")
NP11UP = pd.read_csv("HCG_NP11.csv")
NP11DN = pd.read_csv("LCG_NP11.csv")
PE11NS = pd.read_csv("ICG_PE11.csv")
NP11NS = pd.read_csv("ICG_NP11.csv")
HisSig = pd.read_csv("Histamine_Signature.csv")


# In[45]:


# Create lists from dataframes
PE11UP_List = PE11UP['GeneID'].to_list()
PE11DN_List = PE11DN['GeneID'].to_list()
PE11NS_List = PE11NS['GeneID'].to_list()
NP11UP_List = NP11UP['GeneID'].to_list()
NP11DN_List = NP11DN['GeneID'].to_list()
NP11NS_List = NP11NS['GeneID'].to_list()
HisSig_List = HisSig['GeneID'].to_list()

# Convert lists to sets
PE11UP_Set = set(PE11UP_List)
PE11DN_Set = set(PE11DN_List)
PE110NS_Set = set(PE11NS_List)
NP11UP_Set = set(NP11UP_List)
NP11DN_Set = set(NP11DN_List)
NP1NS_Set = set(NP11NS_List)
HisSig_Set = set(HisSig_List)

# Calculate overlaps
overlap_PE11UP_HisSig = PE11UP_Set.intersection(HisSig_Set)
overlap_PE11DN_HisSig = PE11DN_Set.intersection(HisSig_Set)
overlap_PE11NS_HisSig = PE11NS_Set.intersection(HisSig_Set)
overlap_NP11UP_HisSig = NP11UP_Set.intersection(HisSig_Set)
overlap_NP11DN_HisSig = NP11DN_Set.intersection(HisSig_Set)
overlap_NP11NS_HisSig = NP11NS_Set.intersection(HisSig_Set)

# Print sizes of overlaps
print("Number of overlapping genes between PE11UP and HisSig:", len(overlap_PE11UP_HisSig))
print("Number of overlapping genes between PE11DN and HisSig:", len(overlap_PE11DN_HisSig))
print("Number of overlapping genes between PE11NS and HisSig:", len(overlap_PE11NS_HisSig))
print("Number of overlapping genes between NP11UP and HisSig:", len(overlap_NP11UP_HisSig))
print("Number of overlapping genes between NP11DN and HisSig:", len(overlap_NP11DN_HisSig))
print("Number of overlapping genes between NP11NS and HisSig:", len(overlap_NP11NS_HisSig))


# In[46]:


import matplotlib.pyplot as plt
from matplotlib_venn import venn2

def plot_venn(set1, set2, label1, label2, title):
    plt.figure()
    venn2([set1, set2], (label1, label2))
    plt.title(title)
    plt.show()

# Generate Venn diagrams
plot_venn(PE11UP_Set, HisSig_Set, "PE11UP", "HisSig", "Overlap between PE11UP and HisSig")


# In[47]:



# Generate Venn diagrams
plot_venn(PE11DN_Set, HisSig_Set, "PE11DN", "HisSig", "Overlap between PE11DN and HisSig")


# In[48]:


plot_venn(PE11NS_Set, HisSig_Set, "PE11NS", "HisSig", "Overlap between PE11NS and HisSig")


# In[49]:


plot_venn(NP11UP_Set, HisSig_Set, "NP11UP", "HisSig", "Overlap between NP11UP and HisSig")


# In[50]:


plot_venn(NP11DN_Set, HisSig_Set, "NP11DN", "HisSig", "Overlap between NP11DN and HisSig")


# In[51]:


plot_venn(NP11NS_Set, HisSig_Set, "NP11NS", "HisSig", "Overlap between NP11NS and HisSig")


# In[53]:



def extract_and_save_intersections(df1, df2, column_name, save_filename):
    # Get intersecting GeneIDs
    intersection = set(df1[column_name]).intersection(set(df2[column_name]))
    
    # Filter rows based on intersecting GeneIDs
    result_df = df1[df1[column_name].isin(intersection)]
    
    # Save to CSV
    result_df.to_csv(save_filename, index=False)

# Use the function to extract and save intersections with the HisSig data frame
extract_and_save_intersections(PE11UP, HisSig, 'GeneID', "PE11UP_HisSig_Intersection.csv")
extract_and_save_intersections(PE11DN, HisSig, 'GeneID', "PE11DN_HisSig_Intersection.csv")
extract_and_save_intersections(PE11NS, HisSig, 'GeneID', "PE11NS_HisSig_Intersection.csv")
extract_and_save_intersections(NP11UP, HisSig, 'GeneID', "NP11UP_HisSig_Intersection.csv")
extract_and_save_intersections(NP11DN, HisSig, 'GeneID', "NP11DN_HisSig_Intersection.csv")
extract_and_save_intersections(NP11NS, HisSig, 'GeneID', "NP11NS_HisSig_Intersection.csv")


# In[2]:


import pandas as pd
import matplotlib.pyplot as plt
from matplotlib_venn import venn2

# Load previously saved intersection CSVs
PE11UP_HisSig = pd.read_csv("PE11UP_HisSig_Intersection.csv")
PE11DN_HisSig = pd.read_csv("PE11DN_HisSig_Intersection.csv")
NP11UP_HisSig = pd.read_csv("NP11UP_HisSig_Intersection.csv")
NP11DN_HisSig = pd.read_csv("NP11DN_HisSig_Intersection.csv")

# Overlap PE8DN_HisSig with NP8DN_HisSig and PE8UP_HisSig with NP8UP_HisSig
#extract_and_save_intersections(PE8DN_HisSig, NP8DN_HisSig, 'GeneID', "PE8DN_NP8DN_HisSig_Intersection.csv")
#extract_and_save_intersections(PE8UP_HisSig, NP8UP_HisSig, 'GeneID', "PE8UP_NP8UP_HisSig_Intersection.csv")

# Isolate genes exclusive to PE intersections
def extract_and_save_exclusive(df1, df2, column_name, save_filename):
    exclusive = set(df1[column_name]).difference(set(df2[column_name]))
    result_df = df1[df1[column_name].isin(exclusive)]
    result_df.to_csv(save_filename, index=False)

# Get exclusive genes
extract_and_save_exclusive(PE11UP_HisSig, NP11UP_HisSig, 'GeneID', "Exclusive_PE11UP_HisSig.csv")
extract_and_save_exclusive(PE11DN_HisSig, NP11DN_HisSig, 'GeneID', "Exclusive_PE11DN_HisSig.csv")


def plot_venn(set1, set2, label1, label2, title, save_filename):  # Added save_filename as an argument
    plt.figure()
    venn2([set1, set2], (label1, label2))
    plt.title(title)
    plt.savefig(save_filename)  # Save the figure before displaying
    plt.show()
    
# Optional: Generate Venn diagrams for the new intersections
plot_venn(set(PE11UP_HisSig['GeneID']), set(NP11UP_HisSig['GeneID']), 
          "PE11UP_HisSig", "NP11UP_HisSig", 
          "Overlap between PE11UP_HisSig and NP11UP_HisSig", 
          "PE11UP_NP11UP_HisSig_Venn.png")  # Pass the filename for saving

plot_venn(set(PE11DN_HisSig['GeneID']), set(NP11DN_HisSig['GeneID']), 
          "PE11DN_HisSig", "NP11DN_HisSig", 
          "Overlap between PE11DN_HisSig and NP11DN_HisSig", 
          "PE11DN_NP11DN_HisSig_Venn.png")  # Pass the filename for saving


# In[ ]:




