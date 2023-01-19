#!/usr/bin/env python
# coding: utf-8

# ## Finalize the single cell perturbseq dataset
#
# The output of the Seurat pipeline in 2.process-perturbseq.ipynb is not easily compatible with our downstream tasks.
# This notebook finalizes the input perturbseq (CRISPRi) dataset.
#
# There are four basic steps:
#
# 1. Load, transpose, and clean the gene expression matrix
# 2. Load the single cell file identities (barcodes and metadata)
# 3. Merge
# 4. Output
#
# Lastly, we use pycytominer.aggregate to form bulk profiles from the single cell readouts

# In[1]:


import pathlib
import numpy as np
import pandas as pd

from pycytominer import aggregate


# In[2]:


gse_id = "GSE132080"
perturbseq_data_dir = pathlib.Path("data/perturbseq/")

output_file = pathlib.Path(f"{perturbseq_data_dir}/{gse_id}_final_analytical.tsv.gz")
output_bulk_file = pathlib.Path(
    f"{perturbseq_data_dir}/{gse_id}_bulk_final_analytical.tsv.gz"
)


# In[3]:


# Load and process gene expression data
gene_exp_file = pathlib.Path(f"{perturbseq_data_dir}/{gse_id}_processed_matrix.tsv.gz")

gene_exp_df = (
    pd.read_csv(gene_exp_file, sep="\t", index_col=0)
    .transpose()
    .reset_index()
    .rename({"index": "Metadata_barcode"}, axis="columns")
)

# Pull out the measured genes
gene_features = gene_exp_df.columns.tolist()
gene_features.remove("Metadata_barcode")

gene_exp_df = gene_exp_df.assign(
    Metadata_sequence=[x.split("-")[0] for x in gene_exp_df.Metadata_barcode]
)
gene_exp_df.columns.name = ""

meta_features = ["Metadata_barcode", "Metadata_sequence"]
gene_features = sorted(gene_exp_df.drop(meta_features, axis="columns").columns.tolist())

gene_exp_df = gene_exp_df.reindex(meta_features + gene_features, axis="columns")

print(gene_exp_df.shape)
gene_exp_df.head()


# In[4]:


# Load cell identities
identity_file = pathlib.Path(f"{perturbseq_data_dir}/{gse_id}_cell_identities.csv.gz")
cell_id_df = pd.read_csv(identity_file, sep=",")

cell_id_df.columns = [f"Metadata_{x}" for x in cell_id_df.columns]
cell_id_df = cell_id_df.assign(
    Metadata_gene_identity=[
        str(x).split("_")[0] for x in cell_id_df.Metadata_guide_identity
    ]
)

print(cell_id_df.shape)
cell_id_df.head()


# In[5]:


# Merge single cells with identifiers
sc_df = cell_id_df.merge(
    gene_exp_df,
    how="right",
    right_on="Metadata_barcode",
    left_on="Metadata_cell_barcode",
)

sc_df = sc_df.reset_index().rename({"index": "Metadata_cell_identity"}, axis="columns")
sc_df.Metadata_cell_identity = [f"sc_profile_{x}" for x in sc_df.Metadata_cell_identity]

print(sc_df.shape)
sc_df.head()


# In[6]:


# Write the file to disk
sc_df.to_csv(
    output_file, index=False, sep="\t", compression={"method": "gzip", "mtime": 1}
)


# ## Calculate bulk perturbseq data

# In[7]:


# Perform single cell aggregation into bulk
bulk_df = aggregate(
    population_df=sc_df,
    strata=["Metadata_guide_identity"],
    features=gene_features,
    operation="median",
)

# remove one row with NaN value
bulk_df = bulk_df[~bulk_df["Metadata_guide_identity"].isnull()]

# create a column for the gene
bulk_df = bulk_df.assign(
    Metadata_gene_identity=[x.split("_")[0] for x in bulk_df.Metadata_guide_identity]
).query("Metadata_gene_identity != '*'")

bulk_df = bulk_df.reindex(
    ["Metadata_guide_identity", "Metadata_gene_identity"] + gene_features,
    axis="columns",
)

print(bulk_df.shape)
bulk_df.head()


# In[8]:


# Write the file to disk
bulk_df.to_csv(
    output_bulk_file, index=False, sep="\t", compression={"method": "gzip", "mtime": 1}
)
