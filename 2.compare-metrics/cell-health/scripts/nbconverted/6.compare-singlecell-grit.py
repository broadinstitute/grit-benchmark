#!/usr/bin/env python
# coding: utf-8

# ## Visualize the results of single cell grit
# 
# Note, this is only for the ES2 cells

# In[1]:


import pathlib
import pandas as pd

import plotnine as gg


# In[2]:


figure_dir = pathlib.Path("figures/single_cell")

control_group_genes_cut = ["Chr2", "Luc", "LacZ"]

plates = [
    "SQ00014613",
    "SQ00014614",
    "SQ00014615",
]

results_dir = pathlib.Path("../../1.calculate-metrics/cell-health/results")
results_prefixes = {
    "grit": "cellhealth_single_cell_grit_",
    "phate": "cellhealth_single_cell_phate_embeddings_",
    "umap": "cellhealth_single_cell_umap_embeddings_"
}

files = {prefix: {plate: []} for prefix in results_prefixes for plate in plates}
for plate in plates:
    for prefix in results_prefixes:
        file_prefix = results_prefixes[prefix]
        files[prefix][plate] = pathlib.Path(f"{results_dir}/{file_prefix}{plate}.tsv.gz")

files


# In[3]:


for plate in plates:
    grit_df = pd.read_csv(files["grit"][plate], sep="\t")
    umap_df = pd.read_csv(files["umap"][plate], sep="\t")
    phate_df = pd.read_csv(files["phate"][plate], sep="\t")

    for gene in umap_df.grit_gene.unique():

        umap_gene_df = umap_df.query("grit_gene == @gene").merge(
            grit_df.query("grit_gene == @gene"),
            left_on=["Metadata_cell_identity"],
            right_on=["perturbation"]
        )
        
        phate_gene_df = phate_df.query("grit_gene == @gene").merge(
            grit_df.query("grit_gene == @gene"),
            left_on=["Metadata_cell_identity"],
            right_on=["perturbation"]
        )


        control_perts = (
            umap_gene_df
            .query("Metadata_gene_name in @control_group_genes_cut")
            .Metadata_pert_name
            .unique()
            .tolist()
        )

        test_perts = (
            umap_gene_df
            .query("Metadata_pert_name not in @control_perts")
            .Metadata_pert_name
            .unique()
            .tolist()
        )

        pert_order = test_perts + control_perts
        umap_gene_df.Metadata_pert_name = pd.Categorical(umap_gene_df.Metadata_pert_name, categories=pert_order)
        phate_gene_df.Metadata_pert_name = pd.Categorical(phate_gene_df.Metadata_pert_name, categories=pert_order)

        gene_umap_gg = (
            gg.ggplot(umap_gene_df, gg.aes(x="umap_0", y="umap_1"))
            + gg.geom_point(gg.aes(fill="grit"), size=0.6, stroke=0.1, alpha=0.7)
            + gg.facet_wrap("~Metadata_pert_name")
            + gg.theme_bw()
            + gg.theme(strip_background=gg.element_rect(color="black", fill="#fdfff4"))
        )
        
        fig_file = pathlib.Path(f"{figure_dir}/{gene}_{plate}_umap.png")
        gene_umap_gg.save(fig_file, dpi=500, height=5, width=6)
        
        gene_phate_gg = (
            gg.ggplot(phate_gene_df, gg.aes(x="phate_0", y="phate_1"))
            + gg.geom_point(gg.aes(fill="grit"), size=0.6, stroke=0.1, alpha=0.7)
            + gg.facet_wrap("~Metadata_pert_name")
            + gg.theme_bw()
            + gg.theme(strip_background=gg.element_rect(color="black", fill="#fdfff4"))
        )
        
        fig_file = pathlib.Path(f"{figure_dir}/{gene}_{plate}_phate.png")
        gene_phate_gg.save(fig_file, dpi=500, height=5, width=6)

