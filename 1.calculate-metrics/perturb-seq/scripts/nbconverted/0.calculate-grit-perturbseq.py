#!/usr/bin/env python
# coding: utf-8

# ## Calculate Grit for bulk and single cell perturbseq data
# 
# Also calculate UMAP embeddings for each perturbation at the same time.

# In[1]:


import pathlib
import numpy as np
import pandas as pd

from pycytominer import aggregate
from pycytominer.cyto_utils import infer_cp_features

from cytominer_eval import evaluate

import umap


# In[2]:


np.random.seed(2021)


# In[3]:


gse_id = "GSE132080"
perturbseq_data_dir = pathlib.Path("../../0.download-data/data/perturbseq/")
perturbseq_screen_phenotypes = "paper_supplement/Table_S16_perturb-seq_screen_phenotypes.txt"


# In[4]:


# Load finalized single cell perturbseq data
gene_exp_file = pathlib.Path(f"{perturbseq_data_dir}/{gse_id}_final_analytical.tsv.gz")

sc_gene_exp_df = pd.read_csv(gene_exp_file, sep="\t")
gene_features = [x for x in sc_gene_exp_df if not x.startswith("Metadata_")]

print(sc_gene_exp_df.shape)
sc_gene_exp_df.head()


# In[5]:


print(len(gene_features))


# In[6]:


# Load activities results (bulk)
file = perturbseq_data_dir / perturbseq_screen_phenotypes
activity_df = pd.read_csv(file, sep="\t").rename({"Unnamed: 0": "id"}, axis="columns")

# Create a perturbation column to match with other IDs
activity_df = activity_df.assign(perturbation=activity_df.gene + "_" + activity_df.id)

print(activity_df.shape)
activity_df.head()


# In[7]:


# Load bulk perturbseq data
bulk_file = pathlib.Path(f"{perturbseq_data_dir}/{gse_id}_bulk_final_analytical.tsv.gz")
bulk_df = pd.read_csv(bulk_file, sep="\t")

# Some genes have very small variance still, remove these!
genes_to_retain = (
    pd.DataFrame(bulk_df.var() > 0.001)
    .reset_index()
    .rename({"index": "gene", 0: "keep"}, axis="columns")
    .query("keep")
    .gene
    .tolist()
)

bulk_subset_df = bulk_df.loc[:, ["Metadata_guide_identity"] + genes_to_retain]

# create a column for the gene
bulk_subset_df = (
    bulk_df
    .assign(Metadata_gene_identity=[x.split("_")[0] for x in bulk_subset_df.Metadata_guide_identity])
    .query("Metadata_gene_identity != '*'")
)

print(bulk_subset_df.shape)
bulk_subset_df.head()


# ## Calculate Grit
# 
# ### Bulk profiles

# In[8]:


barcode_col = "Metadata_guide_identity"
gene_col = "Metadata_gene_identity"

replicate_group_grit = {
    "profile_col": barcode_col,
    "replicate_group_col": gene_col
}

neg_controls = [x for x in bulk_subset_df.Metadata_guide_identity if "neg_ctrl" in x]
neg_controls


# In[9]:


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


# In[10]:


# Merge with activity results and output file
output_results_file = pathlib.Path(f"results/{gse_id}_grit.tsv")

result = result.merge(activity_df, left_on="perturbation", right_on="perturbation")

result.to_csv(output_results_file, sep="\t", index=False)

print(result.shape)
result.head(3)


# ### Single cells

# In[11]:


# Determine a proportion of negative control guide population
sc_neg_controls_df = sc_gene_exp_df.query("Metadata_guide_identity in @neg_controls").sample(frac=0.2)

sc_neg_controls = (
    sc_neg_controls_df
    .query("Metadata_guide_identity in @neg_controls")
    .Metadata_cell_identity
    .tolist()
)

replicate_group_grit = {
    "profile_col": "Metadata_cell_identity",
    "replicate_group_col": "Metadata_guide_identity"
}


# In[12]:


all_sc_grit_results = []
all_sc_umap_embeddings = []

genes = sc_gene_exp_df.Metadata_gene_identity.unique()
for gene in genes:
    if gene not in ["neg", "*", "nan", np.nan]:
        print(f"Now analyzing {gene}...")
        subset_sc_df = sc_gene_exp_df.query("Metadata_gene_identity in @gene")
        
        # There are a certain number of guides targeting each gene
        guides = subset_sc_df.Metadata_guide_identity.unique()

        # Use the same controls in every experiment
        subset_sc_df = pd.concat([subset_sc_df, sc_neg_controls_df]).reset_index(drop=True)

        # Apply UMAP to single cell profiles (all profiles of one gene + neg controls)
        embedding = umap.UMAP().fit_transform(subset_sc_df.loc[:, genes_to_retain])
        
        # Combine results with single cell dataframe
        embedding_df = pd.concat(
            [
                subset_sc_df.drop(gene_features, axis="columns").reset_index(drop=True),
                pd.DataFrame(embedding, columns=["umap_0", "umap_1"])
            ],
            axis="columns"
        )
        
        # Append to list
        all_sc_umap_embeddings.append(embedding_df.assign(grit_gene=gene))
        
        # Now calculate sc-Grit per guide
        for guide in guides:
            subset_guide_df = pd.concat(
                [
                    subset_sc_df.query("Metadata_guide_identity == @guide"),
                    sc_neg_controls_df
                ]
            ).reset_index(drop=True)
            
            # Calculate Grit
            # Note, every negative control single cell will recieve MULTIPLE grit scores
            # depending on the replicate group information (replicate_group_col)!
            sc_grit_result = evaluate(
                profiles=subset_guide_df,
                features=genes_to_retain,
                meta_features=["Metadata_guide_identity", "Metadata_cell_identity"],
                replicate_groups=replicate_group_grit,
                operation="grit",
                grit_control_perts=[str(x) for x in sc_neg_controls]
            )

            all_sc_grit_results.append(
                sc_grit_result.assign(grit_gene=gene, grit_guide=guide)
            )


# In[13]:


all_sc_umap_embeddings = pd.concat(all_sc_umap_embeddings).reset_index(drop=True)

# Output file
output_results_file = pathlib.Path(f"results/{gse_id}_single_cell_umap_embeddings.tsv.gz")
all_sc_umap_embeddings.to_csv(output_results_file, sep="\t", compression="gzip", index=False)

print(all_sc_umap_embeddings.shape)
all_sc_umap_embeddings.head()


# In[14]:


all_sc_grit_results = pd.concat(all_sc_grit_results).reset_index(drop=True)

# Output file
output_results_file = pathlib.Path(f"results/{gse_id}_single_cell_grit.tsv.gz")
all_sc_grit_results.to_csv(output_results_file, sep="\t", compression="gzip", index=False)

print(all_sc_grit_results.shape)
all_sc_grit_results.head()

