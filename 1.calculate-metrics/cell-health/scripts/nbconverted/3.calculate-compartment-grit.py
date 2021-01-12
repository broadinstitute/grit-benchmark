#!/usr/bin/env python
# coding: utf-8

# # Calculate grit for different Cell Profiler feature groups
# 
# * Per compartment (Cells, Cytoplasm, Nuclei)
# * Per compartment feature group (Cells x AreaShape, Cells x Correlation, Nuclei x Texture, etc.)
# * Per channel (DNA, RNA, AGP, ER, Mito)
#     * Use all features that include any information from one of these channels

# In[1]:


import pathlib
import numpy as np
import pandas as pd

from pycytominer.cyto_utils import infer_cp_features

from cytominer_eval import evaluate


# In[2]:


compartments = ["Cells", "Cytoplasm", "Nuclei"]


# In[3]:


# Load Cell Health data
data_file = pathlib.Path("data/cell_health_merged_feature_select.csv.gz")

df = pd.read_csv(data_file, sep=",")

print(df.shape)
df.head(2)


# In[4]:


# Define cell health constants
barcode_col = "Metadata_pert_name"
gene_col = "Metadata_gene_name"

replicate_group_grit = {
    "replicate_id": barcode_col,
    "group_id": gene_col
}

control_group_cut = ["Chr2", "Luc", "LacZ"]

control_barcodes = (
    df.loc[
        df[replicate_group_grit["group_id"]].isin(control_group_cut),
        replicate_group_grit["replicate_id"]
    ]
    .unique()
    .tolist()
)

control_barcodes


# In[5]:


all_features = infer_cp_features(df, compartments=compartments)
meta_features = infer_cp_features(df, metadata=True)

meta_features


# In[6]:


grit_compartment_results = []
for cell_line in df.Metadata_cell_line.unique():
    for compartment in compartments:
        compartment_features = infer_cp_features(df, compartments=compartment)
        subset_df = df.loc[:, meta_features + compartment_features]

        result = evaluate(
            profiles=subset_df.query("Metadata_cell_line == @cell_line"),
            features=compartment_features,
            meta_features=[barcode_col, gene_col],
            replicate_groups=replicate_group_grit,
            operation="grit",
            similarity_metric="pearson",
            grit_control_perts=control_barcodes
        ).assign(
            cell_line=cell_line,
            barcode_control="cutting_control",
            cor_method="pearson",
            compartment=compartment,
            channel="all",
            feature_group="all",
            num_features=len(compartment_features)
        )
        
        grit_compartment_results.append(result)
        
grit_compartment_results = pd.concat(grit_compartment_results).reset_index(drop=True)

print(grit_compartment_results.shape)
grit_compartment_results.head()


# ## Calculate grit for feature groups

# In[7]:


feature_group_compartments = list(set(["_".join(x.split("_")[0:2]) for x in all_features]))

grit_subcompartment_results = []
for cell_line in df.Metadata_cell_line.unique():
    for compartment_group in feature_group_compartments:
        compartment_features = df.loc[:, df.columns.str.startswith(compartment_group)].columns.tolist()
        subset_df = df.loc[:, meta_features + compartment_features]
        
        compartment, feature_group = compartment_group.split("_")

        result = evaluate(
            profiles=subset_df.query("Metadata_cell_line == @cell_line"),
            features=compartment_features,
            meta_features=[barcode_col, gene_col],
            replicate_groups=replicate_group_grit,
            operation="grit",
            similarity_metric="pearson",
            grit_control_perts=control_barcodes
        ).assign(
            cell_line=cell_line,
            barcode_control="cutting_control",
            cor_method="pearson",
            compartment=compartment,
            channel="all",
            feature_group=feature_group,
            num_features=len(compartment_features)
        )
        
        grit_subcompartment_results.append(result)
        
grit_subcompartment_results = pd.concat(grit_subcompartment_results).reset_index(drop=True)

print(grit_subcompartment_results.shape)
grit_subcompartment_results.head()


# ## Calculate grit for channels

# In[8]:


channels = ["DNA", "RNA", "Mito", "AGP", "ER"]

grit_channel_results = []
for cell_line in df.Metadata_cell_line.unique():
    for channel in channels:
        channel_features = df.loc[:, df.columns.str.contains(channel)].columns.tolist()

        subset_df = df.loc[:, meta_features + compartment_features]
        
        compartment, feature_group = compartment_group.split("_")

        result = evaluate(
            profiles=subset_df.query("Metadata_cell_line == @cell_line"),
            features=compartment_features,
            meta_features=[barcode_col, gene_col],
            replicate_groups=replicate_group_grit,
            operation="grit",
            similarity_metric="pearson",
            grit_control_perts=control_barcodes
        ).assign(
            cell_line=cell_line,
            barcode_control="cutting_control",
            cor_method="pearson",
            compartment="all",
            channel=channel,
            feature_group="all",
            num_features=len(channel_features)
        )
        
        grit_channel_results.append(result)
        
grit_channel_results = pd.concat(grit_channel_results).reset_index(drop=True)

print(grit_channel_results.shape)
grit_channel_results.head()


# ## Concatenate results together

# In[9]:


full_grit_results = pd.concat(
    [
        grit_compartment_results,
        grit_subcompartment_results,
        grit_channel_results,
    ],
    axis="rows"
).reset_index(drop=True)

print(full_grit_results.shape)
full_grit_results.head()


# In[10]:


# Output results
output_dir = "results"
output_file = pathlib.Path(f"{output_dir}/cell_health_grit_compartments.tsv.gz")

full_grit_results.to_csv(
    output_file,
    sep="\t",
    compression={"method": "gzip", "mtime": 1},
    index=False
)

