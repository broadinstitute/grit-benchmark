#!/usr/bin/env python
# coding: utf-8

# ## Calculate grit in different simulated scenarios
# 
# 1. Calculate shuffled grit
# 2. Change the number of controls used for grit calculation

# In[1]:


import pathlib
import numpy as np
import pandas as pd

from pycytominer.cyto_utils import infer_cp_features

from cytominer_eval import evaluate


# In[2]:


np.random.seed(123)


# In[3]:


# Load Cell Health data
plate = "SQ00014613"
data_file = pathlib.Path("data/cell_health_merged_feature_select.csv.gz")

df = pd.read_csv(data_file, sep=",").query("Metadata_Plate == @plate").reset_index(drop=True)

print(df.shape)
df.head(2)


# ## 1. Calculate grit with shuffled permutations

# In[4]:


compartments = ["Cells", "Cytoplasm", "Nuclei"]
cor_method = "pearson"
num_shuffle_permutations = 5


# In[5]:


# Define grit 
barcode_col = "Metadata_pert_name"
gene_col = "Metadata_gene_name"

replicate_group_grit = {
    "replicate_id": barcode_col,
    "group_id": gene_col
}

control_group_cut = ["Chr2", "Luc", "LacZ"]
control_group_pert = ["EMPTY"]

control_barcodes_cut = (
    df.loc[
        df[replicate_group_grit["group_id"]].isin(control_group_cut),
        replicate_group_grit["replicate_id"]
    ]
    .unique()
    .tolist()
)

control_barcodes_pert = (
    df.loc[
        df[replicate_group_grit["group_id"]].isin(control_group_pert),
        replicate_group_grit["replicate_id"]
    ]
    .unique()
    .tolist()
)

control_barcodes = {
    "cutting_control": control_barcodes_cut,
    "perturbation_control": control_barcodes_pert
}

control_barcodes


# In[6]:


morph_features = infer_cp_features(df, compartments=compartments)
meta_features = infer_cp_features(df, metadata=True)

meta_features


# In[7]:


get_ipython().run_cell_magic('time', '', 'grit_results = []\nfor i in range(0, num_shuffle_permutations):\n    for cell_line in df.Metadata_cell_line.unique():\n    \n        profiles = df.query("Metadata_cell_line == @cell_line").reset_index(drop=True)\n        meta_df = profiles.loc[:, meta_features]\n\n        feature_df = profiles.drop(meta_features, axis="columns")\n\n        shuffle_df = pd.concat(\n            [\n                meta_df,\n                feature_df.apply(\n                    lambda x: np.random.permutation(x.values)\n                ).reset_index(drop=True),\n            ],\n            axis="columns"\n        )\n\n        result = evaluate(\n            profiles=shuffle_df,\n            features=morph_features,\n            meta_features=[barcode_col, gene_col],\n            replicate_groups=replicate_group_grit,\n            operation="grit",\n            similarity_metric=cor_method,\n            grit_control_perts=control_barcodes["cutting_control"]\n        ).assign(\n            cell_line=cell_line,\n            barcode_control="cutting_control",\n            cor_method=cor_method,\n            random_iteration=i,\n            shuffle_method="independent_column"\n        )\n        grit_results.append(result)\n        \n\n        shuffle_df = df.copy()\n        shuffle_df.loc[:, ["Metadata_gene_name", "Metadata_pert_name"]] = (\n            df.loc[:, ["Metadata_gene_name", "Metadata_pert_name"]]\n            .sample(n=df.shape[0])\n            .reset_index(drop=True)\n        )\n        \n        result = evaluate(\n            profiles=shuffle_df,\n            features=morph_features,\n            meta_features=[barcode_col, gene_col],\n            replicate_groups=replicate_group_grit,\n            operation="grit",\n            similarity_metric=cor_method,\n            grit_control_perts=control_barcodes["cutting_control"]\n        ).assign(\n            cell_line=cell_line,\n            barcode_control="cutting_control",\n            cor_method=cor_method,\n            random_iteration=i,\n            shuffle_method="metadata_shuffle"\n        )\n        grit_results.append(result)\n        \n\n# Calculate plate wise grit using real data\nresult = evaluate(\n    profiles=df,\n    features=morph_features,\n    meta_features=[barcode_col, gene_col],\n    replicate_groups=replicate_group_grit,\n    operation="grit",\n    similarity_metric=cor_method,\n    grit_control_perts=control_barcodes["cutting_control"]\n).assign(\n    cell_line=cell_line,\n    barcode_control="cutting_control",\n    cor_method=cor_method,\n    random_iteration="real",\n    shuffle_method="real",\n)\n\ngrit_results.append(result)\n\ngrit_results = pd.concat(grit_results).reset_index(drop=True)\n\nprint(grit_results.shape)\ngrit_results.head()')


# In[8]:


# Output results
output_dir = "results"
output_file = pathlib.Path(f"{output_dir}/cell_health_grit_randomshuffled_{plate}.tsv")

grit_results.to_csv(output_file, sep="\t", index=False)


# ## 2. Titrate the amount of control perturbations present when calculating grit

# In[9]:


control_ns = [56, 50, 40, 30, 20, 15, 10, 7, 5, 4, 3, 2]


# In[10]:


get_ipython().run_cell_magic('time', '', 'grit_results = []\nfor i in range(0, num_shuffle_permutations):\n    for n in control_ns:\n        for cell_line in df.Metadata_cell_line.unique():\n\n            profiles = df.query("Metadata_cell_line == @cell_line").reset_index(drop=True)\n\n            control_df = (\n                profiles\n                .query("Metadata_pert_name in @control_barcodes[\'perturbation_control\']")\n            )\n            treatment_df = (\n                profiles\n                .query("Metadata_pert_name not in @control_barcodes[\'perturbation_control\']")\n            )\n\n            control_dropped_profiles_df = pd.concat(\n                [\n                    treatment_df,\n                    control_df.sample(n)\n                ], axis="rows"\n            )\n            \n            result = evaluate(\n                profiles=control_dropped_profiles_df,\n                features=morph_features,\n                meta_features=[barcode_col, gene_col],\n                replicate_groups=replicate_group_grit,\n                operation="grit",\n                similarity_metric=cor_method,\n                grit_control_perts=control_barcodes["perturbation_control"]\n            ).assign(\n                cell_line=cell_line,\n                barcode_control="perturbation_control",\n                cor_method=cor_method,\n                random_iteration=i,\n                num_controls=n\n            )\n\n            grit_results.append(result)\n\ngrit_results = pd.concat(grit_results).reset_index(drop=True)\n\nprint(grit_results.shape)\ngrit_results.head()')


# In[11]:


# Output results
output_dir = "results"
output_file = pathlib.Path(f"{output_dir}/cell_health_grit_control_titration_{plate}.tsv")

grit_results.to_csv(output_file, sep="\t", index=False)

