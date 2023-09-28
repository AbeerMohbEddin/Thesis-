#!/usr/bin/env python
# coding: utf-8

# In[5]:


import os
os.chdir("C:/Users/abeer/Desktop/week16")


# In[2]:


import os

directory_path = 'C:/Users/abeer/Desktop/week16'  # Replace with your directory path
files = os.listdir(directory_path)

for file in files:
    print(file)


# In[3]:


import pandas as pd

files_to_convert = [
    "HCG_NP16.xlsx",
    "HCG_PE16.xlsx",
    "ICG_NP16.xlsx",
    "ICG_PE16.xlsx",
    "LCG_NP16.xlsx",
    "LCG_PE16.xlsx"
]

for file in files_to_convert:
    # Read the Excel file
    data = pd.read_excel(file)
    
    # Convert to CSV
    csv_name = file.replace(".xlsx", ".csv")
    data.to_csv(csv_name, index=False)


# In[5]:


import pandas as pd

PE16UP = pd.read_csv("HCG_PE16.csv")
PE16DN = pd.read_csv("LCG_PE16.csv")
NP16UP = pd.read_csv("HCG_NP16.csv")
NP16DN = pd.read_csv("LCG_NP16.csv")
PE16NS = pd.read_csv("ICG_PE16.csv")
NP16NS = pd.read_csv("ICG_NP16.csv")
HisSig = pd.read_csv("Histamine_Signature.csv")


# In[6]:


# Create lists from dataframes
PE16UP_List = PE16UP['GeneID'].to_list()
PE16DN_List = PE16DN['GeneID'].to_list()
PE16NS_List = PE16NS['GeneID'].to_list()
NP16UP_List = NP16UP['GeneID'].to_list()
NP16DN_List = NP16DN['GeneID'].to_list()
NP16NS_List = NP16NS['GeneID'].to_list()
HisSig_List = HisSig['GeneID'].to_list()

# Convert lists to sets
PE16UP_Set = set(PE16UP_List)
PE16DN_Set = set(PE16DN_List)
PE16NS_Set = set(PE16NS_List)
NP16UP_Set = set(NP16UP_List)
NP16DN_Set = set(NP16DN_List)
NP16NS_Set = set(NP16NS_List)
HisSig_Set = set(HisSig_List)

# Calculate overlaps
overlap_PE16UP_HisSig = PE16UP_Set.intersection(HisSig_Set)
overlap_PE16DN_HisSig = PE16DN_Set.intersection(HisSig_Set)
overlap_PE16NS_HisSig = PE16NS_Set.intersection(HisSig_Set)
overlap_NP16UP_HisSig = NP16UP_Set.intersection(HisSig_Set)
overlap_NP16DN_HisSig = NP16DN_Set.intersection(HisSig_Set)
overlap_NP16NS_HisSig = NP16NS_Set.intersection(HisSig_Set)

# Print sizes of overlaps
print("Number of overlapping genes between PE16UP and HisSig:", len(overlap_PE16UP_HisSig))
print("Number of overlapping genes between PE16DN and HisSig:", len(overlap_PE16DN_HisSig))
print("Number of overlapping genes between PE16NS and HisSig:", len(overlap_PE16NS_HisSig))
print("Number of overlapping genes between NP16UP and HisSig:", len(overlap_NP16UP_HisSig))
print("Number of overlapping genes between NP16DN and HisSig:", len(overlap_NP16DN_HisSig))
print("Number of overlapping genes between NP16NS and HisSig:", len(overlap_NP16NS_HisSig))


# In[7]:


import matplotlib.pyplot as plt
from matplotlib_venn import venn2

def plot_venn(set1, set2, label1, label2, title):
    plt.figure()
    venn2([set1, set2], (label1, label2))
    plt.title(title)
    plt.show()

# Generate Venn diagrams
plot_venn(PE16UP_Set, HisSig_Set, "PE16UP", "HisSig", "Overlap between PE16UP and HisSig")


# In[7]:



# Generate Venn diagrams
plot_venn(PE15DN_Set, HisSig_Set, "PE15DN", "HisSig", "Overlap between PE15DN and HisSig")


# In[8]:


plot_venn(PE16NS_Set, HisSig_Set, "PE16NS", "HisSig", "Overlap between PE16NS and HisSig")


# In[9]:


plot_venn(NP16UP_Set, HisSig_Set, "NP16UP", "HisSig", "Overlap between NP16UP and HisSig")


# In[10]:


plot_venn(NP16DN_Set, HisSig_Set, "NP16DN", "HisSig", "Overlap between NP16DN and HisSig")


# In[11]:


plot_venn(NP16NS_Set, HisSig_Set, "NP16NS", "HisSig", "Overlap between NP16NS and HisSig")


# In[12]:


def extract_and_save_intersections(df1, df2, column_name, save_filename):
    # Get intersecting GeneIDs
    intersection = set(df1[column_name]).intersection(set(df2[column_name]))
    
    # Filter rows based on intersecting GeneIDs
    result_df = df1[df1[column_name].isin(intersection)]
    
    # Save to CSV
    result_df.to_csv(save_filename, index=False)

# Use the function to extract and save intersections with the HisSig data frame
extract_and_save_intersections(PE16UP, HisSig, 'GeneID', "PE16UP_HisSig_Intersection.csv")
extract_and_save_intersections(PE16DN, HisSig, 'GeneID', "PE16DN_HisSig_Intersection.csv")
extract_and_save_intersections(PE16NS, HisSig, 'GeneID', "PE16NS_HisSig_Intersection.csv")
extract_and_save_intersections(NP16UP, HisSig, 'GeneID', "NP16UP_HisSig_Intersection.csv")
extract_and_save_intersections(NP16DN, HisSig, 'GeneID', "NP16DN_HisSig_Intersection.csv")
extract_and_save_intersections(NP16NS, HisSig, 'GeneID', "NP16NS_HisSig_Intersection.csv")


# In[6]:


import pandas as pd
import matplotlib.pyplot as plt
from matplotlib_venn import venn2

# Load previously saved intersection CSVs
PE16UP_HisSig = pd.read_csv("PE16UP_HisSig_Intersection.csv")
PE16DN_HisSig = pd.read_csv("PE16DN_HisSig_Intersection.csv")
NP16UP_HisSig = pd.read_csv("NP16UP_HisSig_Intersection.csv")
NP16DN_HisSig = pd.read_csv("NP16DN_HisSig_Intersection.csv")


# Isolate genes exclusive to PE intersections
def extract_and_save_exclusive(df1, df2, column_name, save_filename):
    exclusive = set(df1[column_name]).difference(set(df2[column_name]))
    result_df = df1[df1[column_name].isin(exclusive)]
    result_df.to_csv(save_filename, index=False)

# Get exclusive genes
extract_and_save_exclusive(PE16UP_HisSig, NP16UP_HisSig, 'GeneID', "Exclusive_PE16UP_HisSig.csv")
extract_and_save_exclusive(PE16DN_HisSig, NP16DN_HisSig, 'GeneID', "Exclusive_PE16DN_HisSig.csv")


def plot_venn(set1, set2, label1, label2, title, save_filename):  # Added save_filename as an argument
    plt.figure()
    venn2([set1, set2], (label1, label2))
    plt.title(title)
    plt.savefig(save_filename)  # Save the figure before displaying
    plt.show()
    
# Optional: Generate Venn diagrams for the new intersections
plot_venn(set(PE16UP_HisSig['GeneID']), set(NP16UP_HisSig['GeneID']), 
          "PE16UP_HisSig", "NP16UP_HisSig", 
          "Overlap between PE16UP_HisSig and NP16UP_HisSig", 
          "PE16UP_NP16UP_HisSig_Venn.png")  # Pass the filename for saving

plot_venn(set(PE16DN_HisSig['GeneID']), set(NP16DN_HisSig['GeneID']), 
          "PE16DN_HisSig", "NP16DN_HisSig", 
          "Overlap between PE16DN_HisSig and NP16DN_HisSig", 
          "PE16DN_NP16DN_HisSig_Venn.png")  # Pass the filename for saving


# In[ ]:




