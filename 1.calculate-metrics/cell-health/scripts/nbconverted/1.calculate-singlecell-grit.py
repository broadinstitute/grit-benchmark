#!/usr/bin/env python
# coding: utf-8

# ## Calculate Grit for Single Cell Profiles in Cell Health Data
# 
# We only calculate grit for ES2 cells.

# In[1]:


import pathlib
import numpy as np
import pandas as pd

import umap
import phate

from pycytominer.cyto_utils import infer_cp_features
from pycytominer import feature_select

from cytominer_eval import evaluate
from cytominer_eval.transform import metric_melt
from cytominer_eval.operations.util import assign_replicates


# In[2]:


np.random.seed(2021)


# In[3]:


# Set grit constants
sample_frac = 0.02
control_group_genes_cut = ["Chr2", "Luc", "LacZ"]

exclude_grit_genes = control_group_genes_cut + ["EMPTY"]

replicate_group_grit = {
    "replicate_id": "Metadata_cell_identity",
    "group_id": "Metadata_pert_name"
}


# In[4]:


# Only process ES2 plates
sc_dir = pathlib.Path("../../0.download-data/data/cell_health/normalized/")

plates = [
    "SQ00014613",
    "SQ00014614",
    "SQ00014615",
]

plate_files = {plate: pathlib.Path(f"{sc_dir}/{plate}_normalized_featureselected.csv.gz") for plate in plates}
plate_files


# In[5]:


for plate in plate_files:
    print(f"Now processing {plate}...")
    sc_df = (
        pd.read_csv(plate_files[plate], low_memory=False)
        .reset_index()
        .rename({"index": "Metadata_cell_identity"}, axis="columns")
    )

    sc_df.loc[:, "Metadata_cell_identity"] = [f"cell_{x}" for x in sc_df.Metadata_cell_identity]

    print(sc_df.shape)
    
    morph_features = infer_cp_features(sc_df)
    
    neg_controls_df = (
        sc_df
        .query("Metadata_gene_name in @control_group_genes_cut")
        .sample(frac=sample_frac)
        .reset_index(drop=True)
    )

    control_group_guides_cut = neg_controls_df.Metadata_pert_name.unique()
    sc_neg_control_cells = neg_controls_df.Metadata_cell_identity.tolist()
    
    # Prepare variables for results storage
    all_sc_grit_results = []
    all_sc_umap_embeddings = []
    all_sc_phate_embeddings = []

    genes = sc_df.Metadata_gene_name.unique()
    for gene in genes:
        if gene in exclude_grit_genes:
            continue
        
        print(f"Now analyzing {gene}...")
            
        subset_sc_df = sc_df.query("Metadata_gene_name in @gene")
        subset_sc_df = pd.concat([subset_sc_df, neg_controls_df]).reset_index(drop=True)
        
        guides = subset_sc_df.Metadata_pert_name.unique()

        # Apply UMAP
        embedding = umap.UMAP().fit_transform(subset_sc_df.loc[:, morph_features])

        # Combine results with single cell dataframe
        embedding_df = pd.concat(
            [
                subset_sc_df.drop(morph_features, axis="columns").reset_index(drop=True),
                pd.DataFrame(embedding, columns=["umap_0", "umap_1"])
            ],
            axis="columns"
        )
        all_sc_umap_embeddings.append(embedding_df.assign(grit_gene=gene))
        
        # Apply PHATE
        phate_operator = phate.PHATE(n_jobs=-2)
        phate_operator.set_params(decay=20, t="auto", gamma=0, verbose=0)

        Y_phate = phate_operator.fit_transform(subset_sc_df.loc[:, morph_features])
        
        # Combine results with single cell dataframe
        phate_embedding_df = pd.concat(
            [
                subset_sc_df.drop(morph_features, axis="columns").reset_index(drop=True),
                pd.DataFrame(Y_phate, columns=["phate_0", "phate_1"])
            ],
            axis="columns"
        )
        all_sc_phate_embeddings.append(phate_embedding_df.assign(grit_gene=gene))

        # Now calculate sc-Grit per guide
        for guide in guides:
            if guide in control_group_guides_cut:
                continue
    
            subset_guide_df = pd.concat(
                [
                    subset_sc_df.query("Metadata_pert_name == @guide"),
                    neg_controls_df
                ]
            ).reset_index(drop=True)
            
            # Calculate Grit
            sc_grit_result = evaluate(
                profiles=subset_guide_df,
                features=morph_features,
                meta_features=["Metadata_pert_name", "Metadata_cell_identity"],
                replicate_groups=replicate_group_grit,
                operation="grit",
                grit_control_perts=sc_neg_control_cells
            ).assign(gene=gene, guide=guide)
            
            all_sc_grit_results.append(
                sc_grit_result.assign(grit_gene=gene, grit_guide=guide)
            )
            
    # Output results
    all_sc_umap_embeddings = pd.concat(all_sc_umap_embeddings).reset_index(drop=True)
    output_results_file = pathlib.Path(f"results/cellhealth_single_cell_umap_embeddings_{plate}.tsv.gz")
    all_sc_umap_embeddings.to_csv(output_results_file, sep="\t", compression="gzip", index=False)

    all_sc_phate_embeddings = pd.concat(all_sc_phate_embeddings).reset_index(drop=True)
    output_results_file = pathlib.Path(f"results/cellhealth_single_cell_phate_embeddings_{plate}.tsv.gz")
    all_sc_phate_embeddings.to_csv(output_results_file, sep="\t", compression="gzip", index=False)

    all_sc_grit_results = pd.concat(all_sc_grit_results).reset_index(drop=True)
    output_results_file = pathlib.Path(f"results/cellhealth_single_cell_grit_{plate}.tsv.gz")
    all_sc_grit_results.to_csv(output_results_file, sep="\t", compression="gzip", index=False)

