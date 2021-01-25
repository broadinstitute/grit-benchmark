#!/usr/bin/env python
# coding: utf-8

# # Calculate grit for bulk profile data

# In[1]:


import pathlib
import numpy as np
import pandas as pd

from pycytominer.cyto_utils import infer_cp_features, output
from pycytominer import feature_select

from cytominer_eval import evaluate
from cytominer_eval.transform import metric_melt
from cytominer_eval.operations.util import assign_replicates


# In[2]:


# Load Cell Health data
commit = "07e4b40c39dd27084be36fbef4d64c5654b2960f"
base_url = f"https://github.com/broadinstitute/cell-health/raw/{commit}"
url = f"{base_url}/1.generate-profiles/data/processed/cell_health_profiles_merged.tsv.gz"

df = pd.read_csv(url, sep="\t")

print(df.shape)
df.head(2)


# In[3]:


# Perform feature selection
feature_select_ops = [
    "variance_threshold",
    "correlation_threshold",
    "drop_na_columns",
    "blocklist",
    "drop_outliers",
]

df = feature_select(
    profiles=df,
    operation=feature_select_ops,
    na_cutoff=0
)

features = infer_cp_features(df)
meta_features = infer_cp_features(df, metadata=True)

print(df.shape)
df.head(2)


# In[4]:


# Output feature selected file
output_file = pathlib.Path("data/cell_health_merged_feature_select.csv.gz")

output(
    df=df,
    output_filename=output_file,
    sep=",",
    compression_options={"method": "gzip", "mtime": 1}
)


# In[5]:


# Define cell health constants
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


get_ipython().run_cell_magic('time', '', 'grit_results = []\nfor cell_line in df.Metadata_cell_line.unique():\n    for control_barcode in control_barcodes:\n        for cor_method in ["pearson", "spearman"]:\n            result = evaluate(\n                profiles=df.query("Metadata_cell_line == @cell_line"),\n                features=features,\n                meta_features=[barcode_col, gene_col],\n                replicate_groups=replicate_group_grit,\n                operation="grit",\n                similarity_metric=cor_method,\n                grit_control_perts=control_barcodes[control_barcode]\n            ).assign(\n                cell_line=cell_line,\n                barcode_control=control_barcode,\n                cor_method=cor_method\n            )\n\n            grit_results.append(result)\n    \ngrit_results = pd.concat(grit_results).reset_index(drop=True)\n\nprint(grit_results.shape)\ngrit_results.head()')


# In[7]:


# Some perturbations have only one guide per gene, these cannot have grit scores
print(grit_results.grit.isna().sum())
grit_results.loc[grit_results.grit.isna(), :].reset_index(drop=True).head(5)


# In[8]:


# Output results
output_dir = "results"
output_file = pathlib.Path(f"{output_dir}/cell_health_grit.tsv")

grit_results.to_csv(output_file, sep="\t", index=False)


# ## Calculate mp-value

# In[9]:


get_ipython().run_cell_magic('time', '', 'mp_results = []\n\nfor cell_line in df.Metadata_cell_line.unique():\n    for num_permutations in [10, 100, 1000, 5000]:\n\n        mp_value_params = {"nb_permutations": num_permutations}\n\n        result = evaluate(\n            profiles=df.query("Metadata_cell_line == @cell_line"),\n            features=features,\n            meta_features=[barcode_col, gene_col],\n            replicate_groups="Metadata_pert_name",\n            operation="mp_value",\n            grit_control_perts=control_barcodes["cutting_control"],\n            mp_value_params=mp_value_params\n        ).assign(\n            cell_line=cell_line,\n            barcode_control="cutting_control",\n            num_permutations=num_permutations\n        )\n\n        mp_results.append(result)\n    \nmp_results = pd.concat(mp_results).reset_index(drop=True)\n\nprint(mp_results.shape)\nmp_results.head()')


# In[10]:


# Output results
output_file = pathlib.Path(f"{output_dir}/cell_health_mpvalue.tsv")

mp_results.to_csv(output_file, sep="\t", index=False)

