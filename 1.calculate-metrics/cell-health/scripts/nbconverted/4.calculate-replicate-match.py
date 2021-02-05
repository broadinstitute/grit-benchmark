#!/usr/bin/env python
# coding: utf-8

# ## Calculate replicate matching

# In[1]:


import pathlib
import numpy as np
import pandas as pd

from pycytominer.cyto_utils import infer_cp_features

from cytominer_eval import evaluate
from cytominer_eval.transform import metric_melt
from cytominer_eval.operations.util import assign_replicates


# In[2]:


output_dir = "results"

file = pathlib.Path("data/cell_health_merged_feature_select.csv.gz")
cell_health_df = pd.read_csv(file)

print(cell_health_df.shape)
cell_health_df.head()


# In[3]:


features = infer_cp_features(cell_health_df)
meta_features = infer_cp_features(cell_health_df, metadata=True)

similarity_metric = "pearson"
operation = "percent_strong"

replicate_groups = ["Metadata_cell_line", "Metadata_gene_name", "Metadata_pert_name"]

control_ids = ["Chr2", "Luc", "LacZ"]


# In[4]:


# Melt the input profiles to long format
similarity_melted_df = metric_melt(
    df=cell_health_df,
    features=features,
    metadata_features=meta_features,
    similarity_metric=similarity_metric,
    eval_metric=operation,
)

similarity_melted_df = assign_replicates(
    similarity_melted_df=similarity_melted_df, replicate_groups=replicate_groups
)

print(similarity_melted_df.shape)
similarity_melted_df.head()


# In[5]:


non_replicate_cor_95th = (
    similarity_melted_df
    .query("not group_replicate")
    .groupby("Metadata_cell_line_pair_a")["similarity_metric"]
    .quantile(0.95)
    .reset_index()
    .rename({"Metadata_cell_line_pair_a": "cell_line"}, axis="columns")
)

# Output results
output_file = pathlib.Path(f"{output_dir}/cell_health_nonreplicate_95thpercentile.tsv")

non_replicate_cor_95th.to_csv(output_file, sep="\t", index=False)
non_replicate_cor_95th


# In[6]:


# Capture median replicate correlation
median_replicate_correlation_df = (
    similarity_melted_df
    .query("group_replicate")
    .groupby(
        [
            "Metadata_cell_line_pair_a",
            "Metadata_gene_name_pair_a",
            "Metadata_pert_name_pair_a"
        ]
    )
    ["similarity_metric"]
    .median()
    .reset_index()
    .rename(
        {
            "similarity_metric": "median_replicate_correlation",
            "Metadata_pert_name_pair_a": "perturbation",
            "Metadata_gene_name_pair_a": "group",
            "Metadata_cell_line_pair_a": "cell_line"
        }, axis="columns"
    )
)

print(median_replicate_correlation_df.shape)
median_replicate_correlation_df.head()


# In[7]:


median_empty_correlation_df = (
    similarity_melted_df
    .query("Metadata_gene_name_pair_b in @control_ids")
    .groupby(
        [
            "Metadata_cell_line_pair_a",
            "Metadata_gene_name_pair_a",
            "Metadata_pert_name_pair_a"
        ]
    )
    ["similarity_metric"]
    .median()
    .reset_index()
    .rename(
        {
            "similarity_metric": "median_control_correlation",
            "Metadata_pert_name_pair_a": "perturbation",
            "Metadata_gene_name_pair_a": "group",
            "Metadata_cell_line_pair_a": "cell_line"
        }, axis="columns"
    )
)

print(median_empty_correlation_df.shape)
median_empty_correlation_df.head()


# In[8]:


full_correlation_results_df = (
    median_replicate_correlation_df
    .merge(
        median_empty_correlation_df,
        on=[
            "cell_line",
            "group",
            "perturbation"
        ],
        how="inner"
    )
)

print(full_correlation_results_df.shape)
full_correlation_results_df.head()


# In[9]:


# Output results
output_file = pathlib.Path(f"{output_dir}/cell_health_replicate_reproducibility.tsv")

full_correlation_results_df.to_csv(output_file, sep="\t", index=False)

