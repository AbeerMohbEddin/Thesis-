#!/usr/bin/env python
# coding: utf-8

# In[1]:


import os
os.chdir("C:/Users/abeer/Desktop/week12")


# In[13]:


import os

directory_path = 'C:/Users/abeer/Desktop/week12'  # Replace with your directory path
files = os.listdir(directory_path)

for file in files:
    print(file)


# In[4]:


import pandas as pd

files_to_convert = [
    "HCG_NP12.xlsx",
    "HCG_PE12.xlsx",
    "ICG_NP12.xlsx",
    "ICG_PE12.xlsx",
    "LCG_NP12.xlsx",
    "LCG_PE12.xlsx"
]

for file in files_to_convert:
    # Read the Excel file
    data = pd.read_excel(file)
    
    # Convert to CSV
    csv_name = file.replace(".xlsx", ".csv")
    data.to_csv(csv_name, index=False)


# In[15]:


import pandas as pd

PE12UP = pd.read_csv("HCG_PE12.csv")
PE12DN = pd.read_csv("LCG_PE12.csv")
NP12UP = pd.read_csv("HCG_NP12.csv")
NP12DN = pd.read_csv("LCG_NP12.csv")
PE12NS = pd.read_csv("ICG_PE12.csv")
NP12NS = pd.read_csv("ICG_NP12.csv")
HisSig = pd.read_csv("Histamine_Signature.csv")


# In[16]:


# Create lists from dataframes
PE12UP_List = PE12UP['GeneID'].to_list()
PE12DN_List = PE12DN['GeneID'].to_list()
PE12NS_List = PE12NS['GeneID'].to_list()
NP12UP_List = NP12UP['GeneID'].to_list()
NP12DN_List = NP12DN['GeneID'].to_list()
NP12NS_List = NP12NS['GeneID'].to_list()
HisSig_List = HisSig['GeneID'].to_list()

# Convert lists to sets
PE12UP_Set = set(PE12UP_List)
PE12DN_Set = set(PE12DN_List)
PE12NS_Set = set(PE12NS_List)
NP12UP_Set = set(NP12UP_List)
NP12DN_Set = set(NP12DN_List)
NP12NS_Set = set(NP12NS_List)
HisSig_Set = set(HisSig_List)

# Calculate overlaps
overlap_PE12UP_HisSig = PE12UP_Set.intersection(HisSig_Set)
overlap_PE12DN_HisSig = PE12DN_Set.intersection(HisSig_Set)
overlap_PE12NS_HisSig = PE12NS_Set.intersection(HisSig_Set)
overlap_NP12UP_HisSig = NP12UP_Set.intersection(HisSig_Set)
overlap_NP12DN_HisSig = NP12DN_Set.intersection(HisSig_Set)
overlap_NP12NS_HisSig = NP12NS_Set.intersection(HisSig_Set)

# Print sizes of overlaps
print("Number of overlapping genes between PE12UP and HisSig:", len(overlap_PE12UP_HisSig))
print("Number of overlapping genes between PE12DN and HisSig:", len(overlap_PE12DN_HisSig))
print("Number of overlapping genes between PE12NS and HisSig:", len(overlap_PE12NS_HisSig))
print("Number of overlapping genes between NP12UP and HisSig:", len(overlap_NP12UP_HisSig))
print("Number of overlapping genes between NP12DN and HisSig:", len(overlap_NP12DN_HisSig))
print("Number of overlapping genes between NP12NS and HisSig:", len(overlap_NP12NS_HisSig))


# In[17]:


import matplotlib.pyplot as plt
from matplotlib_venn import venn2

def plot_venn(set1, set2, label1, label2, title):
    plt.figure()
    venn2([set1, set2], (label1, label2))
    plt.title(title)
    plt.show()

# Generate Venn diagrams
plot_venn(PE12UP_Set, HisSig_Set, "PE12UP", "HisSig", "Overlap between PE12UP and HisSig")


# In[18]:



# Generate Venn diagrams
plot_venn(PE12DN_Set, HisSig_Set, "PE12DN", "HisSig", "Overlap between PE12DN and HisSig")


# In[19]:


plot_venn(PE12NS_Set, HisSig_Set, "PE12NS", "HisSig", "Overlap between PE12NS and HisSig")


# In[20]:


plot_venn(NP12UP_Set, HisSig_Set, "NP12UP", "HisSig", "Overlap between NP12UP and HisSig")


# In[21]:


plot_venn(NP12DN_Set, HisSig_Set, "NP12DN", "HisSig", "Overlap between NP12DN and HisSig")


# In[23]:


plot_venn(NP12NS_Set, HisSig_Set, "NP12NS", "HisSig", "Overlap between NP12NS and HisSig")


# In[24]:



def extract_and_save_intersections(df1, df2, column_name, save_filename):
    # Get intersecting GeneIDs
    intersection = set(df1[column_name]).intersection(set(df2[column_name]))
    
    # Filter rows based on intersecting GeneIDs
    result_df = df1[df1[column_name].isin(intersection)]
    
    # Save to CSV
    result_df.to_csv(save_filename, index=False)

# Use the function to extract and save intersections with the HisSig data frame
extract_and_save_intersections(PE12UP, HisSig, 'GeneID', "PE12UP_HisSig_Intersection.csv")
extract_and_save_intersections(PE12DN, HisSig, 'GeneID', "PE12DN_HisSig_Intersection.csv")
extract_and_save_intersections(PE12NS, HisSig, 'GeneID', "PE12NS_HisSig_Intersection.csv")
extract_and_save_intersections(NP12UP, HisSig, 'GeneID', "NP12UP_HisSig_Intersection.csv")
extract_and_save_intersections(NP12DN, HisSig, 'GeneID', "NP12DN_HisSig_Intersection.csv")
extract_and_save_intersections(NP12NS, HisSig, 'GeneID', "NP12NS_HisSig_Intersection.csv")


# In[2]:


import pandas as pd
import matplotlib.pyplot as plt
from matplotlib_venn import venn2

# Load previously saved intersection CSVs
PE12UP_HisSig = pd.read_csv("PE12UP_HisSig_Intersection.csv")
PE12DN_HisSig = pd.read_csv("PE12DN_HisSig_Intersection.csv")
NP12UP_HisSig = pd.read_csv("NP12UP_HisSig_Intersection.csv")
NP12DN_HisSig = pd.read_csv("NP12DN_HisSig_Intersection.csv")

# Overlap PE8DN_HisSig with NP8DN_HisSig and PE8UP_HisSig with NP8UP_HisSig
#extract_and_save_intersections(PE8DN_HisSig, NP8DN_HisSig, 'GeneID', "PE8DN_NP8DN_HisSig_Intersection.csv")
#extract_and_save_intersections(PE8UP_HisSig, NP8UP_HisSig, 'GeneID', "PE8UP_NP8UP_HisSig_Intersection.csv")

# Isolate genes exclusive to PE intersections
def extract_and_save_exclusive(df1, df2, column_name, save_filename):
    exclusive = set(df1[column_name]).difference(set(df2[column_name]))
    result_df = df1[df1[column_name].isin(exclusive)]
    result_df.to_csv(save_filename, index=False)

# Get exclusive genes
extract_and_save_exclusive(PE12UP_HisSig, NP12UP_HisSig, 'GeneID', "Exclusive_PE12UP_HisSig.csv")
extract_and_save_exclusive(PE12DN_HisSig, NP12DN_HisSig, 'GeneID', "Exclusive_PE12DN_HisSig.csv")


def plot_venn(set1, set2, label1, label2, title, save_filename):  # Added save_filename as an argument
    plt.figure()
    venn2([set1, set2], (label1, label2))
    plt.title(title)
    plt.savefig(save_filename)  # Save the figure before displaying
    plt.show()
    
# Optional: Generate Venn diagrams for the new intersections
plot_venn(set(PE12UP_HisSig['GeneID']), set(NP12UP_HisSig['GeneID']), 
          "PE12UP_HisSig", "NP12UP_HisSig", 
          "Overlap between PE12UP_HisSig and NP12UP_HisSig", 
          "PE12UP_NP12UP_HisSig_Venn.png")  # Pass the filename for saving

plot_venn(set(PE12DN_HisSig['GeneID']), set(NP12DN_HisSig['GeneID']), 
          "PE12DN_HisSig", "NP12DN_HisSig", 
          "Overlap between PE12DN_HisSig and NP12DN_HisSig", 
          "PE12DN_NP12DN_HisSig_Venn.png")  # Pass the filename for saving


# In[ ]:




