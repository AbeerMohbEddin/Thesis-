#!/usr/bin/env python
# coding: utf-8

# In[1]:


import os
os.chdir("C:/Users/abeer/Desktop/week14")


# In[2]:


import os

directory_path = 'C:/Users/abeer/Desktop/week14'  # Replace with your directory path
files = os.listdir(directory_path)

for file in files:
    print(file)


# In[3]:


import pandas as pd

files_to_convert = [
    "HCG_NP14.xlsx",
    "HCG_PE14.xlsx",
    "ICG_NP14.xlsx",
    "ICG_PE14.xlsx",
    "LCG_NP14.xlsx",
    "LCG_PE14.xlsx"
]

for file in files_to_convert:
    # Read the Excel file
    data = pd.read_excel(file)
    
    # Convert to CSV
    csv_name = file.replace(".xlsx", ".csv")
    data.to_csv(csv_name, index=False)


# In[4]:


import pandas as pd

PE14UP = pd.read_csv("HCG_PE14.csv")
PE14DN = pd.read_csv("LCG_PE14.csv")
NP14UP = pd.read_csv("HCG_NP14.csv")
NP14DN = pd.read_csv("LCG_NP14.csv")
PE14NS = pd.read_csv("ICG_PE14.csv")
NP14NS = pd.read_csv("ICG_NP14.csv")
HisSig = pd.read_csv("Histamine_Signature.csv")


# In[5]:


# Create lists from dataframes
PE14UP_List = PE14UP['GeneID'].to_list()
PE14DN_List = PE14DN['GeneID'].to_list()
PE14NS_List = PE14NS['GeneID'].to_list()
NP14UP_List = NP14UP['GeneID'].to_list()
NP14DN_List = NP14DN['GeneID'].to_list()
NP14NS_List = NP14NS['GeneID'].to_list()
HisSig_List = HisSig['GeneID'].to_list()

# Convert lists to sets
PE14UP_Set = set(PE14UP_List)
PE14DN_Set = set(PE14DN_List)
PE14NS_Set = set(PE14NS_List)
NP14UP_Set = set(NP14UP_List)
NP14DN_Set = set(NP14DN_List)
NP14NS_Set = set(NP14NS_List)
HisSig_Set = set(HisSig_List)

# Calculate overlaps
overlap_PE14UP_HisSig = PE14UP_Set.intersection(HisSig_Set)
overlap_PE14DN_HisSig = PE14DN_Set.intersection(HisSig_Set)
overlap_PE14NS_HisSig = PE14NS_Set.intersection(HisSig_Set)
overlap_NP14UP_HisSig = NP14UP_Set.intersection(HisSig_Set)
overlap_NP14DN_HisSig = NP14DN_Set.intersection(HisSig_Set)
overlap_NP14NS_HisSig = NP14NS_Set.intersection(HisSig_Set)

# Print sizes of overlaps
print("Number of overlapping genes between PE14UP and HisSig:", len(overlap_PE14UP_HisSig))
print("Number of overlapping genes between PE14DN and HisSig:", len(overlap_PE14DN_HisSig))
print("Number of overlapping genes between PE14NS and HisSig:", len(overlap_PE14NS_HisSig))
print("Number of overlapping genes between NP14UP and HisSig:", len(overlap_NP14UP_HisSig))
print("Number of overlapping genes between NP14DN and HisSig:", len(overlap_NP14DN_HisSig))
print("Number of overlapping genes between NP14NS and HisSig:", len(overlap_NP14NS_HisSig))


# In[6]:


import matplotlib.pyplot as plt
from matplotlib_venn import venn2

def plot_venn(set1, set2, label1, label2, title):
    plt.figure()
    venn2([set1, set2], (label1, label2))
    plt.title(title)
    plt.show()

# Generate Venn diagrams
plot_venn(PE14UP_Set, HisSig_Set, "PE14UP", "HisSig", "Overlap between PE14UP and HisSig")


# In[7]:



# Generate Venn diagrams
plot_venn(PE14DN_Set, HisSig_Set, "PE14DN", "HisSig", "Overlap between PE14DN and HisSig")


# In[8]:


plot_venn(PE14NS_Set, HisSig_Set, "PE14NS", "HisSig", "Overlap between PE14NS and HisSig")


# In[9]:


plot_venn(NP14UP_Set, HisSig_Set, "NP14UP", "HisSig", "Overlap between NP14UP and HisSig")


# In[10]:


plot_venn(NP14DN_Set, HisSig_Set, "NP14DN", "HisSig", "Overlap between NP14DN and HisSig")


# In[11]:


plot_venn(NP14NS_Set, HisSig_Set, "NP14NS", "HisSig", "Overlap between NP14NS and HisSig")


# In[12]:


def extract_and_save_intersections(df1, df2, column_name, save_filename):
    # Get intersecting GeneIDs
    intersection = set(df1[column_name]).intersection(set(df2[column_name]))
    
    # Filter rows based on intersecting GeneIDs
    result_df = df1[df1[column_name].isin(intersection)]
    
    # Save to CSV
    result_df.to_csv(save_filename, index=False)

# Use the function to extract and save intersections with the HisSig data frame
extract_and_save_intersections(PE14UP, HisSig, 'GeneID', "PE14UP_HisSig_Intersection.csv")
extract_and_save_intersections(PE14DN, HisSig, 'GeneID', "PE14DN_HisSig_Intersection.csv")
extract_and_save_intersections(PE14NS, HisSig, 'GeneID', "PE14NS_HisSig_Intersection.csv")
extract_and_save_intersections(NP14UP, HisSig, 'GeneID', "NP14UP_HisSig_Intersection.csv")
extract_and_save_intersections(NP14DN, HisSig, 'GeneID', "NP14DN_HisSig_Intersection.csv")
extract_and_save_intersections(NP14NS, HisSig, 'GeneID', "NP14NS_HisSig_Intersection.csv")


# In[2]:


import pandas as pd
import matplotlib.pyplot as plt
from matplotlib_venn import venn2

# Load previously saved intersection CSVs
PE14UP_HisSig = pd.read_csv("PE14UP_HisSig_Intersection.csv")
PE14DN_HisSig = pd.read_csv("PE14DN_HisSig_Intersection.csv")
NP14UP_HisSig = pd.read_csv("NP14UP_HisSig_Intersection.csv")
NP14DN_HisSig = pd.read_csv("NP14DN_HisSig_Intersection.csv")

# Overlap PE8DN_HisSig with NP8DN_HisSig and PE8UP_HisSig with NP8UP_HisSig
#extract_and_save_intersections(PE8DN_HisSig, NP8DN_HisSig, 'GeneID', "PE8DN_NP8DN_HisSig_Intersection.csv")
#extract_and_save_intersections(PE8UP_HisSig, NP8UP_HisSig, 'GeneID', "PE8UP_NP8UP_HisSig_Intersection.csv")

# Isolate genes exclusive to PE intersections
def extract_and_save_exclusive(df1, df2, column_name, save_filename):
    exclusive = set(df1[column_name]).difference(set(df2[column_name]))
    result_df = df1[df1[column_name].isin(exclusive)]
    result_df.to_csv(save_filename, index=False)

# Get exclusive genes
extract_and_save_exclusive(PE14UP_HisSig, NP14UP_HisSig, 'GeneID', "Exclusive_PE14UP_HisSig.csv")
extract_and_save_exclusive(PE14DN_HisSig, NP14DN_HisSig, 'GeneID', "Exclusive_PE14DN_HisSig.csv")


def plot_venn(set1, set2, label1, label2, title, save_filename):  # Added save_filename as an argument
    plt.figure()
    venn2([set1, set2], (label1, label2))
    plt.title(title)
    plt.savefig(save_filename)  # Save the figure before displaying
    plt.show()
    
# Optional: Generate Venn diagrams for the new intersections
plot_venn(set(PE14UP_HisSig['GeneID']), set(NP14UP_HisSig['GeneID']), 
          "PE14UP_HisSig", "NP14UP_HisSig", 
          "Overlap between PE14UP_HisSig and NP14UP_HisSig", 
          "PE14UP_NP14UP_HisSig_Venn.png")  # Pass the filename for saving

plot_venn(set(PE14DN_HisSig['GeneID']), set(NP14DN_HisSig['GeneID']), 
          "PE14DN_HisSig", "NP14DN_HisSig", 
          "Overlap between PE14DN_HisSig and NP14DN_HisSig", 
          "PE14DN_NP14DN_HisSig_Venn.png")  # Pass the filename for saving


# In[ ]:




