#!/usr/bin/env python
# coding: utf-8

# ## Library Import

# In[1]:


import matplotlib.pyplot as plt
import matplotlib
import seaborn as sns
import pandas as pd

from sklearn.decomposition import PCA

import scanpy as sc
import scipy.stats as sps

import scipy.spatial
import scipy.cluster.hierarchy as sch
from scipy.io import mmread
import numpy as np


# ## Dendrogram Plotting

# In[2]:


df = pd.read_csv('DCM_HCM_UMAP_V1.txt',delimiter="\t",low_memory=False)
#Reading the data file which is tab delimited
display(df)


# In[3]:


df_1 = df.drop(labels = [0],axis = 0)
# Cleaning the dataframe to remove the unnecessary components
df_1.drop('NAME', inplace=True, axis=1)
# This column contains teh barcodes which are used to identify the RNA seq expression data. Since it is not needed
#for our analysis, it is dropped along axis 1(rows)
display(df_1)
columns = list(df_1.columns.values)
#This stores the column header values in an array
print(columns)


# In[4]:


unique = df_1['Category'].unique().tolist()
#Stores the unique values of the categories in a list. This will be used later to create the labels of the dendrogram
unique.sort()
unique_1 = []
for n in range(len(unique)):
    unique_1.append(unique[n][3:])
    #This for loop removes the numbers at the front of the labels for formatting.

print(unique_1)


# In[5]:


df_columns = df_1['Category'].copy()
display(df_columns)
df_2 = df_1.apply(pd.to_numeric, errors='coerce')
#This is done to convert numbers which are represented as strings to float values so that we can find the mean of
#the points later
df_2.pop('Category')
#Since the column category only contains strings which could not be converted into float, they are removed.
display(df_2)
df_3 = pd.concat([df_2,df_columns],axis = 1)
display(df_3)
#This is the dataframe upon which we are going to find the center of the clusters.


# In[6]:


# df_1.groupby(by='Category')[['X','Y']].mean()
# new_df['Category'] = ['00.Cardiomyocyte_I', '09.Adipocyte', '17.Proliferating_macrophage', '16.Cardiomyocyte_III', '04.Macrophage', '08.Endocardial', '01.Fibroblast_I', '15.Endothelial_III', '10.Neuronal', '12.Cardiomyocyte_II', '11.Lymphatic_endothelial', '05.VSMC', '02.Endothelial_I', '06.Endothelial_II', '03.Pericyte_I', '13.Activated_fibroblast', '18.Pericyte_II', '14.Mast_cell', '07.Lymphocyte', '19.Epicardial', '20.Fibroblast_II']

new_df = pd.DataFrame(columns = ['X','Y','Category'])



for i in unique:
    #This entire for loop iterates through the list and gets the categories which are clustered together, extracts
    #the rows of the categories. It then finds the center of the clusters by finding out the mean of all the points
    #inside the columns. The resulting data is stored in a new dataframe.
    df_indv = df_3[df_3['Category'] == i]
    x_mean = df_indv['X'].mean()
    new_df.loc[i,['X']] = [x_mean]
    y_mean = df_indv['Y'].mean()
    new_df.loc[i,['Y']] = [y_mean]
    
    
new_df.pop('Category')
new_df.reset_index(inplace = True)
#Since the index contains the names of the categories we are resetting the index to numbers and moving the category
#names to a new column
new_df.rename(columns = {'index':'Category'},inplace = True)
#The column containing the category is renamed to our preference.
display(new_df)


# In[7]:


change_df = new_df.copy()
nparray = change_df.to_numpy(copy = True)
#A numpy array is created from the dataframe which is used to build the dendrogram. This is done because the scipy
#linkage and dendrogram functions only take arrays as inputs.
print(nparray)
indices = nparray[:,1:3]
#Slicing through the array is necessary to remove the category names from the arrays.
print("\n")
print(indices)


# ## Dendrogram

# In[8]:


cluster = sch.linkage(indices, method = 'average', metric = 'euclidean')
#This creates the linkage matrix upon which the dendrogram is generated
#Method describes the clustering method ued to create the linkage matrix,
#Metric describes how the distance between the points is measured
dendrogram = sch.dendrogram(cluster,labels = unique_1, leaf_rotation = 90)
#This method takes the linkage matrix as the input
#Labels takes the cluster names as the input
#Leaf_rotation rotates the leaves(Clusternames) to 90 degrees to align the names properly
plt.title(label = "Average Method", loc = 'center')
plt.show()


# ## Different methods to obtain Dendrogram Clustering

# In[9]:


# methods = ['single','complete','average','weighted','centroid','median','ward']

# for met in methods:
#     cluster = sch.linkage(indices, method = met, metric = 'euclidean')
#     dendrogram = sch.dendrogram(cluster,labels = unique_1, leaf_rotation = 90)
#     plt.title(label = "%s method"%met , loc = 'center')
#     plt.show()  


# ## Heatmap with Hierarchial Clustering

# In[10]:


heatmap_df = pd.read_excel('/Users/yamuna/Documents/Bioinformatics/Project/Source Data Extended Data Fig.3.xlsx',sheet_name = 'PanelB.HeatMap')
#This creates a dataframe of the data provided for the image.
# heatmap_df[]
heatmap_df = heatmap_df.rename({'Unnamed: 0': 'Genes'}, axis=1)
#Renaming one column to our required name
Gene_names = heatmap_df['Genes']
#This stores the gene names in a list
column_headers = list(heatmap_df.columns.values)
column_headers = column_headers[1:]
#These 2 lines of code are used to store the names of the various cell clusters in a list
# print(column_headers)
# print(Gene_names)
display(heatmap_df)


# In[11]:


heatmap_df.pop('Genes')
#In order to convert the given data into a numpy array which can be taken as an input by the clustemap function of
#sns.clustermap we have to remove the gene column
vectors = heatmap_df.values
#Here we are converting the data frame into a numpy array
vectors


# In[12]:


sns.set_context("paper",font_scale = 1)
#This sets the font size of the letterings on the image


# ## Best Method for Clustermap Recreation

# In[13]:


Z_columns = sch.linkage(np.transpose(vectors),method = 'ward' ,metric='euclidean')
#This is done to create a dendrogram which is identical to the one produced in the clustermap.
#The vector data has to be transposed in order to create a linkage matrix upon which the dendrogram is generated
dendrogram = sch.dendrogram(Z_columns,labels = column_headers, leaf_rotation = 90)
#This method generates the dendrograms.

sns_plot = sns.clustermap(vectors,method = 'ward', metric = 'euclidean', 
                          xticklabels = column_headers, yticklabels = False, 
                          row_cluster = False, col_linkage = Z_columns, 
                          cmap="seismic", cbar_pos = [1.1, 0.5, 0.02, 0.3], 
                          vmin = -3.2, vmax = 3.2, figsize = (10,15))
#This line of code generates the heatmap with the clustering of the variables. Each of the individual parameters
#allow us to modify the plot to generate an image similar to the one in the paper.
#xticklabels and yticklabels takes lists which contain the column and row headers as input and display them in the
#plot
#col_linkage allows us to use a linkage matrix which has already been generated and use it in the clustermap
#cmap allows us to use a predefined colour scale in the clustermap
#cbar_pos is used to define the position of the colorbar
#vmin and vmax defines the minimum and maximum limit of the colorbar
#figsize determines the size of the clustermap image
plt.title("Z-score Expression")
plt.show()


# ## Different methods for Clustermap generation

# In[14]:


# for met in methods:
#     sns_plot = sns.clustermap(vectors,method = met, metric = 'euclidean', 
#                               xticklabels = column_headers, yticklabels = False, 
#                               row_cluster = False, cmap="seismic", cbar_pos = [1.1, 0.5, 0.02, 0.3], 
#                               vmin = -3, vmax = 3)


# ## Clustermap with clustering done on Y axis(Gene Names)

# In[15]:


# sns_plot = sns.clustermap(vectors,method = 'ward', metric = 'euclidean', 
#                           xticklabels = column_headers, yticklabels = False, 
#                           row_cluster = True, col_linkage = Z_columns, 
#                           cmap="seismic", cbar_pos = [1.1, 0.5, 0.02, 0.3], 
#                           vmin = -3.2, vmax = 3.2, figsize = (10,15))
# plt.title("Z-score Expression")
# plt.show()

