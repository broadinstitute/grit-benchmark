#!/usr/bin/env python
# coding: utf-8

# ## Calculate grit with different normalization schemes
#
# We compare whole well vs. control-based normalization in grit calculations.

# In[1]:


import pathlib
import pandas as pd

from pycytominer.cyto_utils import infer_cp_features, output
from pycytominer.operations import get_na_columns

from cytominer_eval import evaluate
from cytominer_eval.transform import metric_melt
from cytominer_eval.operations.util import assign_replicates


# In[2]:


output_dir = "results"


# In[3]:


# Load different normalized data
data_dir = pathlib.Path("../../0.download-data/data/cell-health/profiles")
plate_file = pathlib.Path(
    f"{data_dir}/cell_health_profiles_merged_wholeplate_normalized_featureselected.tsv.gz"
)

profile_df = pd.read_csv(plate_file, sep="\t")

features = infer_cp_features(profile_df)
meta_features = infer_cp_features(profile_df, metadata=True)

print(profile_df.shape)
profile_df.head()


# In[4]:


# Define cell health constants
barcode_col = "Metadata_pert_name"
gene_col = "Metadata_gene_name"

replicate_group_grit = {"replicate_id": barcode_col, "group_id": gene_col}

control_group_cut = ["Chr2", "Luc", "LacZ"]
control_group_pert = ["EMPTY"]

control_barcodes_cut = (
    profile_df.loc[
        profile_df[replicate_group_grit["group_id"]].isin(control_group_cut),
        replicate_group_grit["replicate_id"],
    ]
    .unique()
    .tolist()
)

control_barcodes_pert = (
    profile_df.loc[
        profile_df[replicate_group_grit["group_id"]].isin(control_group_pert),
        replicate_group_grit["replicate_id"],
    ]
    .unique()
    .tolist()
)

control_barcodes = {
    "cutting_control": control_barcodes_cut,
    "perturbation_control": control_barcodes_pert,
}

control_barcodes


# In[5]:


get_ipython().run_cell_magic(
    "time",
    "",
    'grit_results = []\nfor cell_line in profile_df.Metadata_cell_line.unique():\n    for control_barcode in control_barcodes:\n        for cor_method in ["pearson", "spearman"]:\n            result = evaluate(\n                profiles=profile_df.query("Metadata_cell_line == @cell_line"),\n                features=features,\n                meta_features=[barcode_col, gene_col],\n                replicate_groups=replicate_group_grit,\n                operation="grit",\n                similarity_metric=cor_method,\n                grit_control_perts=control_barcodes[control_barcode]\n            ).assign(\n                cell_line=cell_line,\n                barcode_control=control_barcode,\n                cor_method=cor_method\n            )\n\n            grit_results.append(result)\n    \ngrit_results = pd.concat(grit_results).reset_index(drop=True)\n\nprint(grit_results.shape)\ngrit_results.head()',
)


# In[6]:


# Output results
output_dir = "results"
output_file = pathlib.Path(f"{output_dir}/cell_health_grit_wholeplatenormalized.tsv")

grit_results.to_csv(output_file, sep="\t", index=False)
