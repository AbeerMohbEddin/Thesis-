#!/usr/bin/env python
# coding: utf-8

# In[1]:


import os
os.chdir("C:/Users/abeer/Desktop/week9")


# In[21]:


import os

directory_path = 'C:/Users/abeer/Desktop/week9'  # Replace with your directory path
files = os.listdir(directory_path)

for file in files:
    print(file)


# In[5]:


import pandas as pd

files_to_convert = [
    "HCG_NP9.xlsx",
    "HCG_PE9.xlsx",
    "ICG_NP9.xlsx",
    "ICG_PE9.xlsx",
    "LCG_NP9.xlsx",
    "LCG_PE9.xlsx"
]

for file in files_to_convert:
    # Read the Excel file
    data = pd.read_excel(file)
    
    # Convert to CSV
    csv_name = file.replace(".xlsx", ".csv")
    data.to_csv(csv_name, index=False)


# In[2]:


import pandas as pd

PE9UP = pd.read_csv("HCG_PE9.csv")
PE9DN = pd.read_csv("LCG_PE9.csv")
NP9UP = pd.read_csv("HCG_NP9.csv")
NP9DN = pd.read_csv("LCG_NP9.csv")
PE9NS = pd.read_csv("ICG_PE9.csv")
NP9NS = pd.read_csv("ICG_NP9.csv")
HisSig = pd.read_csv("Histamine_Signature.csv")


# In[3]:


# Create lists from dataframes
PE9UP_List = PE9UP['GeneID'].to_list()
PE9DN_List = PE9DN['GeneID'].to_list()
PE9NS_List = PE9NS['GeneID'].to_list()
NP9UP_List = NP9UP['GeneID'].to_list()
NP9DN_List = NP9DN['GeneID'].to_list()
NP9NS_List = NP9NS['GeneID'].to_list()
HisSig_List = HisSig['GeneID'].to_list()

# Convert lists to sets
PE9UP_Set = set(PE9UP_List)
PE9DN_Set = set(PE9DN_List)
PE9NS_Set = set(PE9NS_List)
NP9UP_Set = set(NP9UP_List)
NP9DN_Set = set(NP9DN_List)
NP9NS_Set = set(NP9NS_List)
HisSig_Set = set(HisSig_List)

# Calculate overlaps
overlap_PE9UP_HisSig = PE9UP_Set.intersection(HisSig_Set)
overlap_PE9DN_HisSig = PE9DN_Set.intersection(HisSig_Set)
overlap_PE9NS_HisSig = PE9NS_Set.intersection(HisSig_Set)
overlap_NP9UP_HisSig = NP9UP_Set.intersection(HisSig_Set)
overlap_NP9DN_HisSig = NP9DN_Set.intersection(HisSig_Set)
overlap_NP9NS_HisSig = NP9NS_Set.intersection(HisSig_Set)

# Print sizes of overlaps
print("Number of overlapping genes between PE9UP and HisSig:", len(overlap_PE9UP_HisSig))
print("Number of overlapping genes between PE9DN and HisSig:", len(overlap_PE9DN_HisSig))
print("Number of overlapping genes between PE9NS and HisSig:", len(overlap_PE9NS_HisSig))
print("Number of overlapping genes between NP9UP and HisSig:", len(overlap_NP9UP_HisSig))
print("Number of overlapping genes between NP9DN and HisSig:", len(overlap_NP9DN_HisSig))
print("Number of overlapping genes between NP9NS and HisSig:", len(overlap_NP9NS_HisSig))


# In[4]:


import matplotlib.pyplot as plt
from matplotlib_venn import venn2

def plot_venn(set1, set2, label1, label2, title):
    plt.figure()
    venn2([set1, set2], (label1, label2))
    plt.title(title)
    plt.show()

# Generate Venn diagrams
plot_venn(PE9UP_Set, HisSig_Set, "PE9UP", "HisSig", "Overlap between PE9UP and HisSig")


# In[25]:



# Generate Venn diagrams
plot_venn(PE9DN_Set, HisSig_Set, "PE9DN", "HisSig", "Overlap between PE9DN and HisSig")


# In[26]:


plot_venn(PE9NS_Set, HisSig_Set, "PE9NS", "HisSig", "Overlap between PE9NS and HisSig")


# In[27]:


plot_venn(NP9UP_Set, HisSig_Set, "NP9UP", "HisSig", "Overlap between NP9UP and HisSig")


# In[28]:


plot_venn(NP9DN_Set, HisSig_Set, "NP9DN", "HisSig", "Overlap between NP9DN and HisSig")


# In[29]:


plot_venn(NP9NS_Set, HisSig_Set, "NP9NS", "HisSig", "Overlap between NP9NS and HisSig")


# In[30]:



def extract_and_save_intersections(df1, df2, column_name, save_filename):
    # Get intersecting GeneIDs
    intersection = set(df1[column_name]).intersection(set(df2[column_name]))
    
    # Filter rows based on intersecting GeneIDs
    result_df = df1[df1[column_name].isin(intersection)]
    
    # Save to CSV
    result_df.to_csv(save_filename, index=False)

# Use the function to extract and save intersections with the HisSig data frame
extract_and_save_intersections(PE9UP, HisSig, 'GeneID', "PE9UP_HisSig_Intersection.csv")
extract_and_save_intersections(PE9DN, HisSig, 'GeneID', "PE9DN_HisSig_Intersection.csv")
extract_and_save_intersections(PE9NS, HisSig, 'GeneID', "PE9NS_HisSig_Intersection.csv")
extract_and_save_intersections(NP9UP, HisSig, 'GeneID', "NP9UP_HisSig_Intersection.csv")
extract_and_save_intersections(NP9DN, HisSig, 'GeneID', "NP9DN_HisSig_Intersection.csv")
extract_and_save_intersections(NP9NS, HisSig, 'GeneID', "NP9NS_HisSig_Intersection.csv")


# In[2]:


import pandas as pd
import matplotlib.pyplot as plt
from matplotlib_venn import venn2

# Load previously saved intersection CSVs
PE9UP_HisSig = pd.read_csv("PE9UP_HisSig_Intersection.csv")
PE9DN_HisSig = pd.read_csv("PE9DN_HisSig_Intersection.csv")
NP9UP_HisSig = pd.read_csv("NP9UP_HisSig_Intersection.csv")
NP9DN_HisSig = pd.read_csv("NP9DN_HisSig_Intersection.csv")

# Overlap PE8DN_HisSig with NP8DN_HisSig and PE8UP_HisSig with NP8UP_HisSig
#extract_and_save_intersections(PE8DN_HisSig, NP8DN_HisSig, 'GeneID', "PE8DN_NP8DN_HisSig_Intersection.csv")
#extract_and_save_intersections(PE8UP_HisSig, NP8UP_HisSig, 'GeneID', "PE8UP_NP8UP_HisSig_Intersection.csv")

# Isolate genes exclusive to PE intersections
def extract_and_save_exclusive(df1, df2, column_name, save_filename):
    exclusive = set(df1[column_name]).difference(set(df2[column_name]))
    result_df = df1[df1[column_name].isin(exclusive)]
    result_df.to_csv(save_filename, index=False)

# Get exclusive genes
extract_and_save_exclusive(PE9UP_HisSig, NP9UP_HisSig, 'GeneID', "Exclusive_PE9UP_HisSig.csv")
extract_and_save_exclusive(PE9DN_HisSig, NP9DN_HisSig, 'GeneID', "Exclusive_PE9DN_HisSig.csv")


def plot_venn(set1, set2, label1, label2, title, save_filename):
    plt.figure()
    venn2([set1, set2], (label1, label2))
    plt.title(title)
    plt.savefig(save_filename)  # Save the figure before displaying
    plt.show()

# Optional: Generate Venn diagrams for the new intersections
plot_venn(set(PE9UP_HisSig['GeneID']), set(NP9UP_HisSig['GeneID']), 
          "PE9UP_HisSig", "NP9UP_HisSig", 
          "Overlap between PE9UP_HisSig and NP9UP_HisSig", 
          "PE9UP_NP9UP_HisSig_Venn.png")  # Filename for saving

plot_venn(set(PE9DN_HisSig['GeneID']), set(NP9DN_HisSig['GeneID']), 
          "PE9DN_HisSig", "NP9DN_HisSig", 
          "Overlap between PE9DN_HisSig and NP9DN_HisSig", 
          "PE9DN_NP9DN_HisSig_Venn.png")  # Filename for saving


# In[ ]:




