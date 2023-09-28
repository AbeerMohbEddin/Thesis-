#!/usr/bin/env python
# coding: utf-8

# In[2]:


import os
os.chdir("C:/Users/abeer/Desktop/week8")


# In[21]:


import os

directory_path = 'C:/Users/abeer/Desktop/week8'  # Replace with your directory path
files = os.listdir(directory_path)

for file in files:
    print(file)


# In[22]:


import pandas as pd

files_to_convert = [
    "HCG_NP8.xlsx",
    "HCG_PE8.xlsx",
    "ICG_NP8.xlsx",
    "ICG_PE8.xlsx",
    "LCG_NP8.xlsx",
    "LCG_PE8.xlsx"
]

for file in files_to_convert:
    # Read the Excel file
    data = pd.read_excel(file)
    
    # Convert to CSV
    csv_name = file.replace(".xlsx", ".csv")
    data.to_csv(csv_name, index=False)


# In[2]:


import pandas as pd

PE8UP = pd.read_csv("HCG_PE8.csv")
PE8DN = pd.read_csv("LCG_PE8.csv")
NP8UP = pd.read_csv("HCG_NP8.csv")
NP8DN = pd.read_csv("LCG_NP8.csv")
PE8NS = pd.read_csv("ICG_PE8.csv")
NP8NS = pd.read_csv("ICG_NP8.csv")
HisSig = pd.read_csv("Histamine_Signature.csv")


# In[3]:


# Create lists from dataframes
PE8UP_List = PE8UP['GeneID'].to_list()
PE8DN_List = PE8DN['GeneID'].to_list()
PE8NS_List = PE8NS['GeneID'].to_list()
NP8UP_List = NP8UP['GeneID'].to_list()
NP8DN_List = NP8DN['GeneID'].to_list()
NP8NS_List = NP8NS['GeneID'].to_list()
HisSig_List = HisSig['GeneID'].to_list()

# Convert lists to sets
PE8UP_Set = set(PE8UP_List)
PE8DN_Set = set(PE8DN_List)
PE8NS_Set = set(PE8NS_List)
NP8UP_Set = set(NP8UP_List)
NP8DN_Set = set(NP8DN_List)
NP8NS_Set = set(NP8NS_List)
HisSig_Set = set(HisSig_List)

# Calculate overlaps
overlap_PE8UP_HisSig = PE8UP_Set.intersection(HisSig_Set)
overlap_PE8DN_HisSig = PE8DN_Set.intersection(HisSig_Set)
overlap_PE8NS_HisSig = PE8NS_Set.intersection(HisSig_Set)
overlap_NP8UP_HisSig = NP8UP_Set.intersection(HisSig_Set)
overlap_NP8DN_HisSig = NP8DN_Set.intersection(HisSig_Set)
overlap_NP8NS_HisSig = NP8NS_Set.intersection(HisSig_Set)

# Print sizes of overlaps
print("Number of overlapping genes between PE8UP and HisSig:", len(overlap_PE8UP_HisSig))
print("Number of overlapping genes between PE8DN and HisSig:", len(overlap_PE8DN_HisSig))
print("Number of overlapping genes between PE8NS and HisSig:", len(overlap_PE8NS_HisSig))
print("Number of overlapping genes between NP8UP and HisSig:", len(overlap_NP8UP_HisSig))
print("Number of overlapping genes between NP8DN and HisSig:", len(overlap_NP8DN_HisSig))
print("Number of overlapping genes between NP8NS and HisSig:", len(overlap_NP8NS_HisSig))


# In[4]:


import matplotlib.pyplot as plt
from matplotlib_venn import venn2

def plot_venn(set1, set2, label1, label2, title):
    plt.figure()
    venn2([set1, set2], (label1, label2))
    plt.title(title)
    plt.show()

# Generate Venn diagrams
plot_venn(PE8UP_Set, HisSig_Set, "PE8UP", "HisSig", "Overlap between PE8UP and HisSig")


# In[29]:



# Generate Venn diagrams
plot_venn(PE8DN_Set, HisSig_Set, "PE8DN", "HisSig", "Overlap between PE8DN and HisSig")


# In[30]:


plot_venn(PE8NS_Set, HisSig_Set, "PE8NS", "HisSig", "Overlap between PE8NS and HisSig")


# In[31]:


plot_venn(NP8UP_Set, HisSig_Set, "NP8UP", "HisSig", "Overlap between NP8UP and HisSig")


# In[32]:


plot_venn(NP8DN_Set, HisSig_Set, "NP8DN", "HisSig", "Overlap between NP8DN and HisSig")


# In[33]:


plot_venn(NP8NS_Set, HisSig_Set, "NP8NS", "HisSig", "Overlap between NP8NS and HisSig")


# In[34]:



def extract_and_save_intersections(df1, df2, column_name, save_filename):
    # Get intersecting GeneIDs
    intersection = set(df1[column_name]).intersection(set(df2[column_name]))
    
    # Filter rows based on intersecting GeneIDs
    result_df = df1[df1[column_name].isin(intersection)]
    
    # Save to CSV
    result_df.to_csv(save_filename, index=False)

# Use the function to extract and save intersections with the HisSig data frame
extract_and_save_intersections(PE8UP, HisSig, 'GeneID', "PE8UP_HisSig_Intersection.csv")
extract_and_save_intersections(PE8DN, HisSig, 'GeneID', "PE8DN_HisSig_Intersection.csv")
extract_and_save_intersections(PE8NS, HisSig, 'GeneID', "PE8NS_HisSig_Intersection.csv")
extract_and_save_intersections(NP8UP, HisSig, 'GeneID', "NP8UP_HisSig_Intersection.csv")
extract_and_save_intersections(NP8DN, HisSig, 'GeneID', "NP8DN_HisSig_Intersection.csv")
extract_and_save_intersections(NP8NS, HisSig, 'GeneID', "NP8NS_HisSig_Intersection.csv")


# In[6]:


import pandas as pd
import matplotlib.pyplot as plt
from matplotlib_venn import venn2

# Load previously saved intersection CSVs
PE8UP_HisSig = pd.read_csv("PE8UP_HisSig_Intersection.csv")
PE8DN_HisSig = pd.read_csv("PE8DN_HisSig_Intersection.csv")
NP8UP_HisSig = pd.read_csv("NP8UP_HisSig_Intersection.csv")
NP8DN_HisSig = pd.read_csv("NP8DN_HisSig_Intersection.csv")

# Overlap PE8DN_HisSig with NP8DN_HisSig and PE8UP_HisSig with NP8UP_HisSig
#extract_and_save_intersections(PE8DN_HisSig, NP8DN_HisSig, 'GeneID', "PE8DN_NP8DN_HisSig_Intersection.csv")
#extract_and_save_intersections(PE8UP_HisSig, NP8UP_HisSig, 'GeneID', "PE8UP_NP8UP_HisSig_Intersection.csv")

# Isolate genes exclusive to PE intersections
def extract_and_save_exclusive(df1, df2, column_name, save_filename):
    exclusive = set(df1[column_name]).difference(set(df2[column_name]))
    result_df = df1[df1[column_name].isin(exclusive)]
    result_df.to_csv(save_filename, index=False)

# Get exclusive genes
extract_and_save_exclusive(PE8UP_HisSig, NP8UP_HisSig, 'GeneID', "Exclusive_PE8UP_HisSig.csv")
extract_and_save_exclusive(PE8DN_HisSig, NP8DN_HisSig, 'GeneID', "Exclusive_PE8DN_HisSig.csv")


def plot_venn(set1, set2, label1, label2, title, save_filename):
    plt.figure()
    venn2([set1, set2], (label1, label2))
    plt.title(title)
    plt.savefig(save_filename)  # Save the figure before displaying
    plt.show()

# Optional: Generate Venn diagrams for the new intersections
plot_venn(set(PE8UP_HisSig['GeneID']), set(NP8UP_HisSig['GeneID']), 
          "PE8UP_HisSig", "NP8UP_HisSig", 
          "Overlap between PE8UP_HisSig and NP8UP_HisSig", 
          "PE8UP_NP8UP_HisSig_Venn.png")  # Filename for saving

plot_venn(set(PE8DN_HisSig['GeneID']), set(NP8DN_HisSig['GeneID']), 
          "PE8DN_HisSig", "NP8DN_HisSig", 
          "Overlap between PE8DN_HisSig and NP8DN_HisSig", 
          "PE8DN_NP8DN_HisSig_Venn.png")  # Filename for saving


# In[ ]:




