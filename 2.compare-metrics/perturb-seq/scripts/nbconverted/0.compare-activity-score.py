#!/usr/bin/env python
# coding: utf-8

# ## Compare Grit to perturbseq Relative Activity Score
# 
# The perturbseq data come from the CRISPRi experiment from:
# 
# > Jost, M., Santos, D.A., Saunders, R.A. et al. Titrating gene expression using libraries of systematically attenuated CRISPR guide RNAs. Nat Biotechnol 38, 355–364 (2020). https://doi.org/10.1038/s41587-019-0387-5
# 
# and relative activity is defined as:
# 
# The fold-knockdown of each mismatched variant divided by the fold-knockdown of the perfectly-matched sgRNA.
# 
# ### Also, visualize singlecell-grit

# In[1]:


import pathlib
import numpy as np
import pandas as pd
import plotnine as gg


# In[2]:


# Load perturbseq results
perturbseq_data_dir = pathlib.Path("../../1.calculate-metrics/perturb-seq/results")
gse_id = "GSE132080"
results_file = pathlib.Path(f"{perturbseq_data_dir}/{gse_id}_grit.tsv")

output_dir = "figures"

grit_df = pd.read_csv(results_file, sep="\t")

grit_df.loc[:, "gene"] = pd.Categorical(
    grit_df.gene, categories=grit_df.gene.unique()
)

print(grit_df.shape)
grit_df.head(2)


# In[3]:


# Global view
global_gg = (
    gg.ggplot(grit_df, gg.aes(x="relative_activity_day5", y="grit")) +
    gg.geom_point(size=0.8) +
    gg.theme_bw() +
    gg.xlab("Relative Activity (Day 5)") +
    gg.ylab("Grit") +
    gg.ggtitle(f"{gse_id} (Jost et al. 2020 CRISPRi)")
)

output_file = pathlib.Path(f"{output_dir}/{gse_id}_crispri_grit_relative_activity_comparison.png")
global_gg.save(output_file, dpi=500, height=3.5, width=4)

global_gg


# In[4]:


gene_gg = (
    gg.ggplot(grit_df, gg.aes(x="relative_activity_day5", y="grit")) +
    gg.geom_point(size=0.6) +
    gg.theme_bw() +
    gg.xlab("Relative Activity (Day 5)") +
    gg.ylab("Grit") +
    gg.ggtitle(f"{gse_id} (Jost et al. 2020 CRISPRi)") +
    gg.facet_wrap("~gene") +
    gg.theme(
        strip_background=gg.element_rect(colour="black", fill="#fdfff4"),
        axis_text=gg.element_text(size=6)
    )
)

output_file = pathlib.Path(f"{output_dir}/{gse_id}_crispri_grit_relative_activity_comparison.png")
gene_gg.save(output_file, dpi=500, height=5, width=6)

gene_gg


# ## Visualize single cell grit

# In[5]:


sc_results_file = pathlib.Path(f"{perturbseq_data_dir}/{gse_id}_single_cell_grit.tsv.gz")
sc_df = pd.read_csv(sc_results_file, sep="\t")

sc_df.loc[:, "gene"] = pd.Categorical(
    sc_df.gene, categories=grit_df.gene.unique()
)

sc_df.loc[:, "gene_identity"] = pd.Categorical(
    sc_df.gene_identity, categories=["neg"] + grit_df.gene.unique().tolist()
)

print(sc_df.shape)
sc_df.head()


# In[6]:


sc_df.gene_identity.value_counts()


# In[7]:


global_gg = (
    gg.ggplot(sc_df.dropna(subset=["gene"]), gg.aes(x="relative_activity_day5", y="grit")) +
    gg.geom_density_2d() +
    gg.geom_point(size=0.2, alpha=0.01) +
    gg.theme_bw() +
    gg.xlab("Relative Activity (Day 5)") +
    gg.ylab("Single Cell Grit") +
    gg.ggtitle(f"{gse_id} (Jost et al. 2020 CRISPRi)")
)

output_file = pathlib.Path(f"{output_dir}/{gse_id}_singlecell_crispri_grit_relative_activity_comparison.png")
global_gg.save(output_file, dpi=500, height=5, width=6)

global_gg


# In[8]:


gene_gg = (
    gg.ggplot(sc_df.dropna(subset=["gene"]), gg.aes(x="relative_activity_day5", y="grit")) +
    gg.geom_density_2d() +
    gg.geom_point(alpha=0.05, size=0.1) +
    gg.theme_bw() +
    gg.xlab("Relative Activity (Day 5)") +
    gg.ylab("Grit") +
    gg.ggtitle(f"{gse_id} (Jost et al. 2020 CRISPRi)") +
    gg.facet_wrap("~gene") +
    gg.theme(
        strip_background=gg.element_rect(colour="black", fill="#fdfff4"),
        axis_text=gg.element_text(size=6)
    )
)

output_file = pathlib.Path(
    f"{output_dir}/{gse_id}_singlecell_by_gene_crispri_grit_relative_activity_comparison.png"
)
gene_gg.save(output_file, dpi=500, height=5, width=6)

gene_gg


# In[9]:


umap_gg = (
    gg.ggplot(sc_df, gg.aes(x="umap_0", y="umap_1")) +
    gg.geom_point(gg.aes(color="grit"), alpha=0.2, size=1, stroke=0) +
    gg.facet_wrap("~gene_identity") +
    gg.theme_bw() +
    gg.xlab("UMAP X") +
    gg.ylab("UMAP Y") +
    gg.ggtitle(f"{gse_id} (Jost et al. 2020 CRISPRi)") +
    gg.theme(
        strip_background=gg.element_rect(colour="black", fill="#fdfff4"),
        axis_text=gg.element_text(size=6)
    )
)

output_file = pathlib.Path(f"{output_dir}/{gse_id}_singlecell_umap_grit.png")
umap_gg.save(output_file, dpi=500, height=5, width=6)

umap_gg


# In[10]:


for gene in sc_df.gene.unique():
    if gene in ["neg", "*", "nan", np.nan]:
        continue
    
    gene_embedding_df = sc_df.query("grit_gene == @gene")

    gene_embedding_df = (
        gene_embedding_df
        .assign(
            grit_facet_label=(
                gene_embedding_df.grit_gene + " " + gene_embedding_df.relative_activity_day5.round(3).astype(str)
            )
        )
    )
    gene_embedding_df.loc[gene_embedding_df.gene_identity == "neg", "grit_facet_label"] = "Negative Ctrl"

    facet_order = ["Negative Ctrl"] + [
        f"{gene_embedding_df.grit_gene.unique()[0]} "+ str(x) 
        for x in sorted(gene_embedding_df.relative_activity_day5.dropna().unique().round(3))
    ]
    
    gene_embedding_df.loc[:, "grit_facet_label"] = pd.Categorical(
        gene_embedding_df.grit_facet_label, categories=facet_order
    )
        
    gene_gg = (
        gg.ggplot(gene_embedding_df.dropna(subset=["grit_facet_label"]), gg.aes(x="umap_0", y="umap_1")) +
        gg.geom_point(gg.aes(color="grit"), size=2, stroke=0, alpha=0.5) +
        gg.facet_wrap("~grit_facet_label") +
        gg.theme_bw() +
        gg.theme(
            strip_background=gg.element_rect(colour="black", fill="#fdfff4")
        )
    )
    
    output_file = pathlib.Path(f"{output_dir}/gene_umaps/{gse_id}_{gene}_singlecell_umap_grit.png")
    gene_gg.save(output_file, dpi=500, height=5, width=6)

    print(gene_gg)

