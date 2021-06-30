#!/usr/bin/env python
# coding: utf-8

# # Build consensus signatures (1 signature for each CRISPR guide) from bulk profiles for Cell Health prediction pipeline
# **Generate consensus signatures with median + moderated-z-score**
# 
# Consensus-profiles are generated via:
# MODZ (moderated z-score) 
# 
# *reference: cell-health/1.generate-profiles/2.build-consensus-signatures*

# In[1]:


import os
import glob
import gzip
from pathlib import Path
import pickle
import re

import numpy as np
import pandas as pd

from pycytominer.consensus import modz
from pycytominer import get_na_columns, aggregate

# from scipy.special import softmax 

from pycytominer import aggregate
from pycytominer.cyto_utils import infer_cp_features
from scripts.utils import calculate_weighted_agg


# ## Load Cell Painting Data
# These are individual df's pf well-level profiles (level 3) that are concatatenated into a single file per aggregation method.

# In[2]:


input_folder = 'data/processed/'
output_folder = 'data/profiles/'
method_list = list(set([(x.split("_")[1:][0]).split('.')[0] for x in glob.glob(input_folder+'*.tsv')]))

# since single-cell grit was not calculated for EMPTY wells, we will use median-aggregated well-level profiles 
# for EMPTY perturbations to form the EMPTY consensus profile for the cell health prediction pipeline
empty_list = []
for file in glob.glob(input_folder+'*.tsv'):
    file_cell_line = file.split('/')[-1].split('.')[0].split('_')[0]
    file_method = file.split('/')[-1].split('.')[0].split('_')[1]
    if "EMPTY" in file:
        print(f"adding {file} to list")
        empties_df = (pd.read_csv(file,sep='\t')
                      .assign(Metadata_cellline = file_cell_line, 
                              Metadata_aggmethod = file_method)
                           )
        empty_list.append(empties_df)
#     print(file.split('/')[-1].split('.')[0].split('_')[1])
empty_profiles = pd.concat(empty_list)
empty_profiles = empty_profiles[sorted(empty_profiles, key = lambda x: x not in empty_profiles.filter(like="Metadata").columns)]
print("total shape: ", empty_profiles.shape)
display(empty_profiles.head())


# perform for both well-level aggregation methods (median and grit-informed)
for method in method_list:# ['weighted']: # 
    print(f"for method is: {method}")
    df_list = []
    for file in glob.glob(input_folder+'*.tsv'):
        file_cell_line = file.split('/')[-1].split('.')[0].split('_')[0]
        file_method = file.split('/')[-1].split('.')[0].split('_')[1]
        if method in file and "EMPTY" not in file:
            print(f"adding {file} to {method} df")
            cell_line_df = (pd.read_csv(file,sep='\t')
                            .assign(Metadata_cellline = file_cell_line, Metadata_aggmethod = file_method)
                           )
            df_list.append(cell_line_df)        
    level3profiles = pd.concat(df_list, axis='rows')
    # add in the EMPTY wells
    level3profiles = pd.concat([level3profiles, empty_profiles], axis='rows')
    # reorder the columns
    level3profiles = level3profiles[sorted(level3profiles, key = lambda x: x not in level3profiles.filter(like="Metadata").columns)]
    print(level3profiles.shape)
    display(level3profiles.head())
    print(infer_cp_features(level3profiles, metadata=True))
    
    # Output final merged file (for all cell lines)
    filename = Path(f"{output_folder}cell_health_profiles_{method}_merged.tsv.gz")
    print(f"filename will be: {filename}")
    level3profiles.to_csv(filename, index=False, sep='\t')


# ## Build Consensus Signatures for aggregation methods
# The remainder of this script generates consensus signatures (1 signature for each CRISPR guide perturbation). The remaining cells are
# 1. ...run once with `method='weighted'` to generate consensus signatures using grit-weighted aggregation of single-cell profiles into well-level profiles
# 2. ...run again with `method='median'` to generate consensus signatures using standard median aggregation of single-cells profiles into well-level profiles

# ### Read in well-level profiles

# In[23]:


folder = 'data/profiles/'
method='weighted'
# method='median'


# In[24]:


x_df = pd.read_csv(Path(f"{folder}cell_health_profiles_{method}_merged.tsv.gz"), sep='\t', low_memory=False)
print(x_df.shape)
display(x_df.head())
x_df.groupby(['Metadata_cell_line']).apply(lambda x: len(get_na_columns(x)))


# ### Load Cell Health labels from cell-health/ project

# In[25]:


commit = "8244680d6e6db1a2bc1f709b9dabf7783c4a9670"
base_url = f"https://github.com/broadinstitute/cell-health/raw/{commit}"
url = f"{base_url}/1.generate-profiles/data/labels/normalized_cell_health_labels.tsv"

y_df = pd.read_csv(url, sep='\t').drop(["plate_name", "well_col", "well_row"], axis="columns")

print(y_df.shape)
y_df.head(3)


# ## Determine how many Cell Painting profiles have Cell Health status labels

# In[26]:


x_groupby_cols = ["Metadata_gene_name", "Metadata_pert_name", "Metadata_cell_line"]

x_metacount_df = (
    x_df
    .loc[:, x_groupby_cols]
    .assign(n_measurements=1)
    .groupby(x_groupby_cols)
    .count()
    .reset_index()
    .assign(data_type="cell_painting")
    .merge(x_df.loc[:, x_groupby_cols + ["Metadata_Well", "Metadata_Plate"]],
           how="left",
           on=x_groupby_cols)
)

print(x_metacount_df.shape)
x_metacount_df.head(2)


# In[27]:


# cell health labels
y_groupby_cols = ["guide", "cell_id"]

y_metacount_df = (
    y_df
    .loc[:, y_groupby_cols]
    .assign(n_measurements=1)
    .groupby(y_groupby_cols)
    .count()
    .reset_index()
    .assign(data_type="cell_health")
)

print(y_metacount_df.shape)
y_metacount_df.head(2)


# In[28]:


all_measurements_df = (
    x_metacount_df
    .merge(
        y_metacount_df,
        left_on=["Metadata_pert_name", "Metadata_cell_line"],
        right_on=["guide", "cell_id"],
        suffixes=["_paint", "_health"],
        how="inner")
    .sort_values(by=["Metadata_cell_line", "Metadata_pert_name"])
    .reset_index(drop=True)
    .drop(["Metadata_Well", "guide", "cell_id"], axis="columns")
)

file = os.path.join("results", "all_profile_metadata.tsv")
# all_measurements_df.to_csv(file, sep='\t', index=False)

print(all_measurements_df.shape)
all_measurements_df.head()


# In[29]:


[len(all_measurements_df[x].unique()) for x in all_measurements_df.columns]


# In[30]:


[len(all_measurements_df[x].unique()) for x in all_measurements_df.columns]


# # apply median consensus aggregation...
# since the modz didnt work initially

# ### 1. to Cell Painting Profiles

# In[31]:


x_median_df = aggregate(
    x_df,
    strata=["Metadata_cell_line", "Metadata_pert_name"],
    features="infer",
    operation="median"
)

x_median_df = (
    x_median_df
    .query("Metadata_pert_name in @all_measurements_df.Metadata_pert_name.unique()")
    .query("Metadata_cell_line in @all_measurements_df.Metadata_cell_line.unique()")
    .reset_index(drop=True)
    .reset_index()
    .rename({"index": "Metadata_profile_id"}, axis='columns')
)
x_median_df.Metadata_profile_id = ["profile_{}".format(x) for x in x_median_df.Metadata_profile_id]

print(x_median_df.shape)
x_median_df.head()


# In[32]:


# Output Profile Mapping for Downstream Analysis
profile_id_mapping_df = x_median_df.loc[:, x_median_df.columns.str.startswith("Metadata")]
file = os.path.join("data", "{}_profile_id_metadata_mapping.tsv".format(method))
print(file)
profile_id_mapping_df.to_csv(file, sep='\t', index=False)

print(profile_id_mapping_df.shape)
profile_id_mapping_df.head()


# ### 2. to Cell Health Panel readouts

# In[33]:


cell_health_meta_features = ["cell_id", "guide"]
cell_health_features = y_df.drop(cell_health_meta_features, axis="columns").columns.tolist()
y_meta_merge_cols = ["Metadata_profile_id", "Metadata_pert_name", "Metadata_cell_line"]


# In[34]:


y_median_df = aggregate(
    y_df,
    strata=cell_health_meta_features,
    features=cell_health_features,
    operation="median"
)

print(y_median_df.shape)
y_median_df.head()


# In[35]:


y_median_df = (
    y_median_df
    .reset_index(drop=True)
    .merge(
        x_median_df.loc[:, y_meta_merge_cols],
        left_on=["guide", "cell_id"],
        right_on=["Metadata_pert_name", "Metadata_cell_line"],
        how="right"
    )
)

# Get columns in correct order
y_columns = (
    y_meta_merge_cols +
    y_median_df
    .loc[:, ~y_median_df.columns.str.startswith("Metadata_")]
    .columns
    .tolist()
)

y_median_df = (
    y_median_df
    .loc[:, y_columns]
    .drop(["guide", "cell_id"], axis="columns")
)

print(y_median_df.shape)
y_median_df.head(5)


# In[36]:


# Confirm that matrices are aligned

pd.testing.assert_series_equal(
    x_median_df.Metadata_profile_id,
    y_median_df.Metadata_profile_id,
    check_names=True
)

# Are the guides aligned?
pd.testing.assert_series_equal(
    x_median_df.Metadata_pert_name,
    y_median_df.Metadata_pert_name,
    check_names=True
)

# Are the cells aligned?
pd.testing.assert_series_equal(
    x_median_df.Metadata_cell_line,
    y_median_df.Metadata_cell_line,
    check_names=True
)


# # apply MODZ consensus aggregation

# ### ...to Cell Painting Profiles

# In[37]:


get_ipython().run_cell_magic('time', '', '\nx_consensus_df = modz(\n    x_df,\n    replicate_columns=["Metadata_cell_line", "Metadata_pert_name"],\n    precision=5\n)\n\nx_consensus_df.head()')


# In[38]:


x_consensus_df = (
    x_consensus_df
    .reset_index()
    .query("Metadata_pert_name in @all_measurements_df.Metadata_pert_name.unique()")
    .query("Metadata_cell_line in @all_measurements_df.Metadata_cell_line.unique()")
    .reset_index(drop=True)
    .reset_index()
    .rename(
        {
            "index": "Metadata_profile_id"
        },
        axis='columns'
    )
)
x_consensus_df.Metadata_profile_id = ["profile_{}".format(x) for x in x_consensus_df.Metadata_profile_id]

print(x_consensus_df.shape)
x_consensus_df.head(5)


# ### Cell health assays data

# In[39]:


get_ipython().run_cell_magic('time', '', '\ny_consensus_df = modz(\n    y_df,\n    features=cell_health_features,\n    replicate_columns=cell_health_meta_features,\n    precision=5\n)\n\nprint(y_consensus_df.shape)\ny_consensus_df.head()')


# In[40]:


y_consensus_df = (
    y_consensus_df
    .reset_index()
    .reset_index(drop=True)
    .merge(
        x_consensus_df.loc[:, y_meta_merge_cols],
        left_on=["guide", "cell_id"],
        right_on=["Metadata_pert_name", "Metadata_cell_line"],
        how="right"
    )
    .loc[:, y_columns]
    .drop(["guide", "cell_id"], axis="columns")
)

print(y_consensus_df.shape)
y_consensus_df.head(5)


# In[41]:


# Confirm that matrices are aligned
pd.testing.assert_series_equal(
    x_consensus_df.Metadata_profile_id,
    y_consensus_df.Metadata_profile_id,
    check_names=True
)

# Are the guides aligned?
pd.testing.assert_series_equal(
    x_consensus_df.Metadata_pert_name,
    y_consensus_df.Metadata_pert_name,
    check_names=True
)

# Are the cells aligned?
pd.testing.assert_series_equal(
    x_consensus_df.Metadata_cell_line,
    y_consensus_df.Metadata_cell_line,
    check_names=True
)


# In[42]:


get_ipython().run_cell_magic('time', '', 'consensus_folder = \'data/consensus/\'\n\nfile = Path(consensus_folder, "{}_agg_cell_painting_median.tsv.gz".format(method))\nx_median_df.to_csv(file, sep="\\t", index=False)\n\nfile = Path(consensus_folder, "{}_agg_cell_health_median.tsv.gz".format(method))\ny_median_df.to_csv(file, sep="\\t", index=False)\n\nfile = Path(consensus_folder, "{}_agg_cell_painting_modz.tsv.gz".format(method))\nx_consensus_df.to_csv(file, sep="\\t", index=False)\n\nfile = Path(consensus_folder, "{}_agg_cell_health_modz.tsv.gz".format(method))\ny_consensus_df.to_csv(file, sep="\\t", index=False)')

