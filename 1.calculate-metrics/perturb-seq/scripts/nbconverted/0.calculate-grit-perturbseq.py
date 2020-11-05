#!/usr/bin/env python
# coding: utf-8

# In[1]:


import pathlib
import numpy as np
import pandas as pd

from pycytominer import aggregate
from pycytominer.cyto_utils import infer_cp_features

from cytominer_eval import evaluate

import umap


# In[2]:


# Load perturbseq data
perturbseq_data_dir = pathlib.Path("../../0.download-data/data/perturbseq")
gse_id = "GSE132080"

file = pathlib.Path(perturbseq_data_dir / f"{gse_id}_processed_matrix.tsv.gz")
df = (
    pd.read_csv(file, sep="\t", index_col=0)
    .transpose()
    .reset_index()
    .rename({"index": "barcode"}, axis="columns")
)

# Pull out the measured genes
gene_features = df.columns.tolist()
gene_features.remove("barcode")

df = df.assign(sequence=[x.split("-")[0] for x in df.barcode])

print(df.shape)
df.head()


# In[3]:


# Load activities results
file = pathlib.Path("supplementary/Table_S16_perturb-seq_screen_phenotypes.txt")
activity_df = pd.read_csv(file, sep="\t").rename({"Unnamed: 0": "id"}, axis="columns")

# Create a perturbation column to match with other IDs
activity_df = activity_df.assign(perturbation=activity_df.gene + "_" + activity_df.id)

print(activity_df.shape)
activity_df.head()


# In[4]:


# Load Cell Identities
cell_id_file = pathlib.Path(f"{perturbseq_data_dir}/{gse_id}_cell_identities.csv.gz")
cell_id_df = pd.read_csv(cell_id_file, sep=",")

print(cell_id_df.shape)
cell_id_df.head()


# In[5]:


# Merge single cells with identifiers
sc_df = cell_id_df.merge(df, how="right", right_on="barcode", left_on="cell_barcode")

print(sc_df.shape)
sc_df.head()


# In[6]:


# Perform single cell aggregation into bulk
bulk_df = aggregate(
    population_df=sc_df,
    strata=["guide_identity"],
    features=gene_features,
    operation="median"
)

# Some genes have very small variance still, remove these!
genes_to_retain = (
    pd.DataFrame(bulk_df.var() > 0.001)
    .reset_index()
    .rename({"index": "gene", 0: "keep"}, axis="columns")
    .query("keep")
    .gene
    .tolist()
)

bulk_df = bulk_df.loc[:, ["guide_identity"] + genes_to_retain]

# create a column for the gene
bulk_df = (
    bulk_df
    .assign(gene_identity=[x.split("_")[0] for x in bulk_df.guide_identity])
    .query("gene_identity != '*'")
)

print(bulk_df.shape)
bulk_df.head()


# ## Calculate Grit

# In[7]:


neg_controls = [x for x in bulk_df.guide_identity if "neg_ctrl" in x]

barcode_col = "guide_identity"
gene_col = "gene_identity"

replicate_group_grit = {
    "replicate_id": barcode_col,
    "group_id": gene_col
}

result = evaluate(
    profiles=bulk_df,
    features=genes_to_retain,
    meta_features=[barcode_col, gene_col],
    replicate_groups=replicate_group_grit,
    operation="grit",
    grit_control_perts=neg_controls
)

result = result.dropna().sort_values(by="grit", ascending=False).reset_index(drop=True)

print(result.shape)
result.head(3)


# In[8]:


# Merge with activity results and output file
output_results_file = pathlib.Path(f"results/{gse_id}_grit.tsv")

result = result.merge(activity_df, left_on="perturbation", right_on="perturbation")

result.to_csv(output_results_file, sep="\t", index=False)

print(result.shape)
result.head(3)


# ## Single Cell Grit

# In[9]:


# Prepare single cell data frame for grit calculation
sc_df = sc_df.assign(gene_identity=[str(x).split("_")[0] for x in sc_df.guide_identity])
sc_df = sc_df.reset_index().rename({"index": "cell_identity"}, axis="columns")
neg_controls_df = sc_df.query("guide_identity in @neg_controls").sample(frac=0.2)

sc_neg_controls = (
    neg_controls_df
    .query("guide_identity in @neg_controls")
    .cell_identity
    .tolist()
)

replicate_group_grit = {'replicate_id': 'cell_identity', 'group_id': 'guide_identity'}


# In[10]:


all_sc_grit_results = []

genes = sc_df.gene_identity.unique()
for gene in genes:
    if gene not in ["neg", "*", "nan"]:
        print(f"Now analyzing {gene}...")
        subset_sc_df = sc_df.query("gene_identity in @gene")
        guides = subset_sc_df.guide_identity.unique()
        for guide in guides:
            subset_guide_df = pd.concat(
                [
                    subset_sc_df.query("guide_identity == @guide"),
                    neg_controls_df
                ]
            ).reset_index(drop=True)
            
            sc_grit_result = evaluate(
                profiles=subset_guide_df,
                features=genes_to_retain,
                meta_features=["guide_identity", "cell_identity"],
                replicate_groups=replicate_group_grit,
                operation="grit",
                grit_control_perts=[str(x) for x in sc_neg_controls]
            )
            
            all_sc_grit_results.append(
                sc_grit_result.assign(grit_gene=gene, grit_guide=guide)
            )


# In[11]:


all_sc_grit_results = pd.concat(all_sc_grit_results).reset_index(drop=True)

print(all_sc_grit_results.shape)
all_sc_grit_results.head()


# In[12]:


# Apply UMAP to single cell profiles
embedding = umap.UMAP().fit_transform(sc_df.loc[:, genes_to_retain])


# In[13]:


# Combine results with single cell dataframe
embedding_df = pd.concat(
    [
        sc_df.drop(gene_features, axis="columns").reset_index(drop=True),
        pd.DataFrame(embedding, columns=["umap_0", "umap_1"])
    ],
    axis="columns"
)

embedding_df.cell_identity = embedding_df.cell_identity.astype(str)

embedding_df = embedding_df.merge(
    all_sc_grit_results,
    left_on="cell_identity",
    right_on="perturbation",
    how="right"
).merge(
    activity_df,
    left_on="guide_identity",
    right_on="perturbation",
    how="outer",
    suffixes=["", "_activity"]
)

print(embedding_df.shape)
embedding_df.head()


# In[14]:


# Output file
output_results_file = pathlib.Path(f"results/{gse_id}_single_cell_grit.tsv.gz")
embedding_df.to_csv(output_results_file, sep="\t", compression="gzip", index=False)

