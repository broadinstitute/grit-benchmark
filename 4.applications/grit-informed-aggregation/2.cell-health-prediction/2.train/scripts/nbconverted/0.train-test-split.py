#!/usr/bin/env python
# coding: utf-8

# # Train Test Split

# In[1]:


import os
import numpy as np
import pandas as pd
from pathlib import Path

from sklearn.model_selection import train_test_split
from sklearn.preprocessing import StandardScaler

from pycytominer import get_na_columns


# ### Load the MODZ consensus profiles

# In[2]:


consensus_folder = '../1.generate-profiles/data/consensus/'
method='median'
# method = 'weighted'


# cell painting profliles
x_median_df = (pd.read_csv(Path("{}{}_agg_cell_painting_modz.tsv.gz".format(consensus_folder, method)), sep='\t'))

# the cell health labels
y_median_df = pd.read_csv(Path("{}{}_agg_cell_health_modz.tsv.gz".format(consensus_folder, method)), sep='\t')


# ### Perform the train-test split stratified by cell line

# In[3]:


x_train_df, x_test_df, y_train_df, y_test_df = train_test_split(x_median_df, y_median_df, test_size = 0.15, stratify = x_median_df[['Metadata_cell_line']])


# In[4]:


print(x_train_df.shape, x_test_df.shape)


# In[5]:


training_samples = x_train_df.Metadata_profile_id.tolist()
testing_samples = x_test_df.Metadata_profile_id.tolist()
training_samples == training_samples


# In[6]:


get_ipython().run_cell_magic('time', '', '# save these splits\nsplits_folder = "data/train_test/{}_agg/".format(method)\nif not os.path.exists(splits_folder):\n    os.makedirs(splits_folder)\n\nfile = f"{splits_folder}x_train_modz.tsv.gz"\nx_train_df.to_csv(Path(file), sep="\\t", index=False)\n\nfile = f"{splits_folder}y_train_modz.tsv.gz"\ny_train_df.to_csv(Path(file), sep="\\t", index=False)\n\nfile = f"{splits_folder}x_test_modz.tsv.gz"\nx_test_df.to_csv(Path(file), sep="\\t", index=False)\n\nfile = f"{splits_folder}y_test_modz.tsv.gz"\ny_test_df.to_csv(Path(file), sep="\\t", index=False)')

