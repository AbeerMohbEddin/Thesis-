#!/usr/bin/env python
# coding: utf-8

# In[1]:


import os
os.chdir("C:/Users/abeer/Desktop/week15")


# In[2]:


import os

directory_path = 'C:/Users/abeer/Desktop/week15'  # Replace with your directory path
files = os.listdir(directory_path)

for file in files:
    print(file)


# In[3]:


import pandas as pd

files_to_convert = [
    "HCG_NP15.xlsx",
    "HCG_PE15.xlsx",
    "ICG_NP15.xlsx",
    "ICG_PE15.xlsx",
    "LCG_NP15.xlsx",
    "LCG_PE15.xlsx"
]

for file in files_to_convert:
    # Read the Excel file
    data = pd.read_excel(file)
    
    # Convert to CSV
    csv_name = file.replace(".xlsx", ".csv")
    data.to_csv(csv_name, index=False)


# In[4]:


import pandas as pd

PE15UP = pd.read_csv("HCG_PE15.csv")
PE15DN = pd.read_csv("LCG_PE15.csv")
NP15UP = pd.read_csv("HCG_NP15.csv")
NP15DN = pd.read_csv("LCG_NP15.csv")
PE15NS = pd.read_csv("ICG_PE15.csv")
NP15NS = pd.read_csv("ICG_NP15.csv")
HisSig = pd.read_csv("Histamine_Signature.csv")


# In[5]:


# Create lists from dataframes
PE15UP_List = PE15UP['GeneID'].to_list()
PE15DN_List = PE15DN['GeneID'].to_list()
PE15NS_List = PE15NS['GeneID'].to_list()
NP15UP_List = NP15UP['GeneID'].to_list()
NP15DN_List = NP15DN['GeneID'].to_list()
NP15NS_List = NP15NS['GeneID'].to_list()
HisSig_List = HisSig['GeneID'].to_list()

# Convert lists to sets
PE15UP_Set = set(PE15UP_List)
PE15DN_Set = set(PE15DN_List)
PE15NS_Set = set(PE15NS_List)
NP15UP_Set = set(NP15UP_List)
NP15DN_Set = set(NP15DN_List)
NP15NS_Set = set(NP15NS_List)
HisSig_Set = set(HisSig_List)

# Calculate overlaps
overlap_PE15UP_HisSig = PE15UP_Set.intersection(HisSig_Set)
overlap_PE15DN_HisSig = PE15DN_Set.intersection(HisSig_Set)
overlap_PE15NS_HisSig = PE15NS_Set.intersection(HisSig_Set)
overlap_NP15UP_HisSig = NP15UP_Set.intersection(HisSig_Set)
overlap_NP15DN_HisSig = NP15DN_Set.intersection(HisSig_Set)
overlap_NP15NS_HisSig = NP15NS_Set.intersection(HisSig_Set)

# Print sizes of overlaps
print("Number of overlapping genes between PE15UP and HisSig:", len(overlap_PE15UP_HisSig))
print("Number of overlapping genes between PE15DN and HisSig:", len(overlap_PE15DN_HisSig))
print("Number of overlapping genes between PE15NS and HisSig:", len(overlap_PE15NS_HisSig))
print("Number of overlapping genes between NP15UP and HisSig:", len(overlap_NP15UP_HisSig))
print("Number of overlapping genes between NP15DN and HisSig:", len(overlap_NP15DN_HisSig))
print("Number of overlapping genes between NP15NS and HisSig:", len(overlap_NP15NS_HisSig))


# In[6]:


import matplotlib.pyplot as plt
from matplotlib_venn import venn2

def plot_venn(set1, set2, label1, label2, title):
    plt.figure()
    venn2([set1, set2], (label1, label2))
    plt.title(title)
    plt.show()

# Generate Venn diagrams
plot_venn(PE15UP_Set, HisSig_Set, "PE15UP", "HisSig", "Overlap between PE15UP and HisSig")


# In[7]:



# Generate Venn diagrams
plot_venn(PE15DN_Set, HisSig_Set, "PE15DN", "HisSig", "Overlap between PE15DN and HisSig")


# In[8]:


plot_venn(PE15NS_Set, HisSig_Set, "PE15NS", "HisSig", "Overlap between PE15NS and HisSig")


# In[9]:


plot_venn(NP15UP_Set, HisSig_Set, "NP15UP", "HisSig", "Overlap between NP15UP and HisSig")


# In[10]:


plot_venn(NP15DN_Set, HisSig_Set, "NP15DN", "HisSig", "Overlap between NP15DN and HisSig")


# In[11]:


plot_venn(NP14NS_Set, HisSig_Set, "NP14NS", "HisSig", "Overlap between NP14NS and HisSig")


# In[11]:


def extract_and_save_intersections(df1, df2, column_name, save_filename):
    # Get intersecting GeneIDs
    intersection = set(df1[column_name]).intersection(set(df2[column_name]))
    
    # Filter rows based on intersecting GeneIDs
    result_df = df1[df1[column_name].isin(intersection)]
    
    # Save to CSV
    result_df.to_csv(save_filename, index=False)

# Use the function to extract and save intersections with the HisSig data frame
extract_and_save_intersections(PE15UP, HisSig, 'GeneID', "PE15UP_HisSig_Intersection.csv")
extract_and_save_intersections(PE15DN, HisSig, 'GeneID', "PE15DN_HisSig_Intersection.csv")
extract_and_save_intersections(PE15NS, HisSig, 'GeneID', "PE15NS_HisSig_Intersection.csv")
extract_and_save_intersections(NP15UP, HisSig, 'GeneID', "NP15UP_HisSig_Intersection.csv")
extract_and_save_intersections(NP15DN, HisSig, 'GeneID', "NP15DN_HisSig_Intersection.csv")
extract_and_save_intersections(NP15NS, HisSig, 'GeneID', "NP15NS_HisSig_Intersection.csv")


# In[2]:


import pandas as pd
import matplotlib.pyplot as plt
from matplotlib_venn import venn2

# Load previously saved intersection CSVs
PE15UP_HisSig = pd.read_csv("PE15UP_HisSig_Intersection.csv")
PE15DN_HisSig = pd.read_csv("PE15DN_HisSig_Intersection.csv")
NP15UP_HisSig = pd.read_csv("NP15UP_HisSig_Intersection.csv")
NP15DN_HisSig = pd.read_csv("NP15DN_HisSig_Intersection.csv")

# Overlap PE8DN_HisSig with NP8DN_HisSig and PE8UP_HisSig with NP8UP_HisSig
#extract_and_save_intersections(PE8DN_HisSig, NP8DN_HisSig, 'GeneID', "PE8DN_NP8DN_HisSig_Intersection.csv")
#extract_and_save_intersections(PE8UP_HisSig, NP8UP_HisSig, 'GeneID', "PE8UP_NP8UP_HisSig_Intersection.csv")

# Isolate genes exclusive to PE intersections
def extract_and_save_exclusive(df1, df2, column_name, save_filename):
    exclusive = set(df1[column_name]).difference(set(df2[column_name]))
    result_df = df1[df1[column_name].isin(exclusive)]
    result_df.to_csv(save_filename, index=False)

# Get exclusive genes
extract_and_save_exclusive(PE15UP_HisSig, NP15UP_HisSig, 'GeneID', "Exclusive_PE15UP_HisSig.csv")
extract_and_save_exclusive(PE15DN_HisSig, NP15DN_HisSig, 'GeneID', "Exclusive_PE15DN_HisSig.csv")


def plot_venn(set1, set2, label1, label2, title, save_filename):  # Added save_filename as an argument
    plt.figure()
    venn2([set1, set2], (label1, label2))
    plt.title(title)
    plt.savefig(save_filename)  # Save the figure before displaying
    plt.show()
    
# Optional: Generate Venn diagrams for the new intersections
plot_venn(set(PE15UP_HisSig['GeneID']), set(NP15UP_HisSig['GeneID']), 
          "PE15UP_HisSig", "NP15UP_HisSig", 
          "Overlap between PE15UP_HisSig and NP15UP_HisSig", 
          "PE15UP_NP15UP_HisSig_Venn.png")  # Pass the filename for saving

plot_venn(set(PE15DN_HisSig['GeneID']), set(NP15DN_HisSig['GeneID']), 
          "PE15DN_HisSig", "NP15DN_HisSig", 
          "Overlap between PE15DN_HisSig and NP15DN_HisSig", 
          "PE15DN_NP15DN_HisSig_Venn.png")  # Pass the filename for saving


# In[ ]:




