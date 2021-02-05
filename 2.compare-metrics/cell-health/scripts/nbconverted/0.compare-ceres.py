#!/usr/bin/env python
# coding: utf-8

# ## Compare CERES readout

# In[1]:


import pathlib
import pandas as pd
import plotnine as gg

import mygene


# In[2]:


output_dir = pathlib.Path("figures")
output_dir.mkdir(exist_ok=True)
cell_health_dir = pathlib.Path("../../1.calculate-metrics/cell-health/results")


# In[3]:


# Load cell health grit scores
cell_health_grit_file = pathlib.Path(f"{cell_health_dir}/cell_health_grit.tsv")

cell_health_grit_df = pd.read_csv(cell_health_grit_file, sep="\t")
print(cell_health_grit_df.shape)
cell_health_grit_df.head()


# In[4]:


mg = mygene.MyGeneInfo()

result = mg.querymany(
    cell_health_grit_df.group.unique().tolist(),
    scopes="symbol",
    species="human",
    fields="entrezgene,symbol,ensembl.gene,",
    as_dataframe=True
)

ncbi_id_df = result.drop_duplicates(subset="_id").loc[:, ["_id"]].reset_index(drop=False)
ncbi_id_df.head(2)


# In[5]:


cell_health_grit_df = cell_health_grit_df.merge(ncbi_id_df, left_on="group", right_on="query", how="left")

print(cell_health_grit_df.shape)
cell_health_grit_df.head(2)


# In[6]:


# Load ceres data
ceres_dir = pathlib.Path("../../0.download-data/data/")
ceres_file = pathlib.Path(f"{ceres_dir}/ceres.csv")
depmap_sample_file = pathlib.Path(f"{ceres_dir}/depmap_sample_info.csv")

ceres_df = pd.read_csv(ceres_file, index_col=0)
depmap_sample_df = pd.read_csv(depmap_sample_file, index_col=0)


# In[7]:


# Clean gene name column
hgnc = [x.split(" ")[0] for x in ceres_df.columns]
ncbi_gene_id = [x.split(" ")[1].strip("()") for x in ceres_df.columns]

ceres_genes_df = (
    pd.DataFrame(
        [hgnc, ncbi_gene_id], 
        index=["HGNC", "NCBI"]
    )
    .transpose()
    .sort_values(by="HGNC")
    .reset_index(drop=True)
)

# Set column names to NCBI
ceres_df.columns = ncbi_gene_id

ceres_genes_df.head(3)


# In[8]:


# Merge the data
ceres_df = depmap_sample_df.merge(ceres_df, left_index=True, right_index=True, how="right")

print(ceres_df.shape)
ceres_df.head(3)


# In[9]:


cell_lines = cell_health_grit_df.cell_line.unique().tolist()
assert all([x in ceres_df.stripped_cell_line_name.tolist() for x in cell_lines])


# In[10]:


cols = ["stripped_cell_line_name"] + cell_health_grit_df._id.dropna().unique().tolist()
ceres_subset_df = (
    ceres_df
    .query("stripped_cell_line_name in @cell_lines")
    .loc[:, cols]
    .reset_index()
    .melt(
        id_vars=["DepMap_ID", "stripped_cell_line_name"],
        var_name="_id",
        value_name="ceres_score"
    )
)

print(ceres_subset_df.shape)
ceres_subset_df.head(3)


# In[11]:


cell_health_results_df = (
    cell_health_grit_df.merge(
        ceres_subset_df,
        left_on=["_id", "cell_line"],
        right_on=["_id", "stripped_cell_line_name"],
        how="left"
    )
)

cell_health_results_df = (
    cell_health_results_df.merge(
        (
            cell_health_results_df
            .groupby(["cell_line", "group"])["grit"]
            .mean()
            .reset_index()
            .rename({"grit": "grit_mean"}, axis="columns")
        ),
        left_on=["cell_line", "group"],
        right_on=["cell_line", "group"],
        how="left"
    )
    .sort_values(by="grit", ascending=False)
)

print(cell_health_results_df.shape)
cell_health_results_df.head(3)


# In[12]:


cell_line_colors = {
  "A549": "#861613",
  "ES2": "#1CADA8",
  "HCC44": "#2A364D"
}

grit_ceres_comparison_gg = (
    gg.ggplot(
        cell_health_results_df,
        gg.aes(x="grit", y="ceres_score")
    ) +
    gg.geom_point(gg.aes(fill="cell_line"), size=1, stroke=0.2) +
    gg.scale_fill_manual(name="Cell Line", values=cell_line_colors) +
    gg.theme_bw() +
    gg.xlab("Grit") +
    gg.ylab("CERES") +
    gg.ggtitle("Cell Health") +
    gg.facet_wrap("~barcode_control", ncol=2) +
    gg.theme(strip_background=gg.element_rect(color="black", fill="#fdfff4"))
)

output_file = pathlib.Path(f"{output_dir}/cell_health_grit_ceres_comparison.png")
grit_ceres_comparison_gg.save(output_file, dpi=500, height=3.5, width=6)

grit_ceres_comparison_gg


# In[13]:


(
    cell_health_results_df
    .query("ceres_score > -1.1")
    .query("grit_mean > 2")
    .query("barcode_control == 'cutting_control'")
    .sort_values(by="ceres_score", ascending=False)
    .reset_index(drop=True)
)


# In[14]:


# What are the perturbations with high grit and low ceres scores?
(
    cell_health_results_df
    .query("grit > 2")
    .query("ceres_score > -1")
    .query("barcode_control == 'cutting_control'")
    .sort_values(by="ceres_score")
    .reset_index(drop=True)
)


# In[15]:


control_compare_df = (
    cell_health_results_df.pivot(
        index=["perturbation", "group", "cell_line", "ceres_score"],
        columns=["barcode_control"],
        values="grit"
    ).reset_index(drop=False)
)

control_compare_df.head()


# In[16]:


grit_barcode_comparison_gg = (
    gg.ggplot(
        control_compare_df,
        gg.aes(x="cutting_control", y="perturbation_control")
    ) +
    gg.geom_point(gg.aes(fill="cell_line"), size=1, stroke=0.2) +
    gg.scale_fill_manual(name="Cell Line", values=cell_line_colors) +
    gg.theme_bw() +
    gg.xlab("Grit (cutting control)") +
    gg.ylab("Grit (perturbation control (empty))") +
    gg.ggtitle("Cell health barcode control comparison") +
    gg.geom_abline(intercept=0, slope=1, linetype="dashed", color="red") +
    gg.coord_fixed()
)

output_file = pathlib.Path(f"{output_dir}/cell_health_barcode_control_comparison.png")
grit_barcode_comparison_gg.save(output_file, dpi=500, height=3.5, width=4)

grit_barcode_comparison_gg

