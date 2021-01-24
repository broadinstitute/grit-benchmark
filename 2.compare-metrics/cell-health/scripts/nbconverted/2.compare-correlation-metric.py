#!/usr/bin/env python
# coding: utf-8

# ## Compare Pearson and Spearman correlations
# 
# In grit, the user determines which similarity metric to use to compare profiles.
# Here, we empirically compare Pearson and Spearman metrics to determine if either provides any benefit over the other.

# In[1]:


import pathlib
import pandas as pd
import plotnine as gg


# In[2]:


output_dir = pathlib.Path("figures")

# Load cell health grit scores
cell_health_dir = pathlib.Path("../../1.calculate-metrics/cell-health/results")
cell_health_grit_file = pathlib.Path(f"{cell_health_dir}/cell_health_grit.tsv")

cell_health_grit_df = pd.read_csv(cell_health_grit_file, sep="\t")
print(cell_health_grit_df.shape)
cell_health_grit_df.head()


# In[3]:


cor_metric_df = (
    cell_health_grit_df
    .pivot(
        index=["perturbation", "group", "cell_line", "barcode_control"],
        columns="cor_method",
        values="grit"
    )
    .reset_index()
)

cor_metric_df = (
    cor_metric_df.assign(differential=cor_metric_df.pearson-cor_metric_df.spearman)
)

print(cor_metric_df.shape)
cor_metric_df.head()


# In[4]:


cell_line_colors = {
  "A549": "#861613",
  "ES2": "#1CADA8",
  "HCC44": "#2A364D"
}

grit_cor_comparison_gg = (
    gg.ggplot(cor_metric_df.dropna(), gg.aes(x="pearson", y="spearman"))
    + gg.geom_point(gg.aes(fill="cell_line"), size=0.7, stroke=0.2, alpha=0.7)
    + gg.facet_wrap("~barcode_control")
    + gg.geom_abline(intercept=0, slope=1, linetype="dotted", color="red")
    + gg.scale_fill_manual(name="Cell Line", values=cell_line_colors)
    + gg.xlab("Pearson correlation")
    + gg.ylab("Spearman correlation")
    + gg.coord_fixed()
    + gg.theme_bw()
    + gg.theme(strip_background=gg.element_rect(color="black", fill="#fdfff4"))
)

output_file = pathlib.Path(f"{output_dir}/cell_health_grit_correlation_metric_comparison.png")
grit_cor_comparison_gg.save(output_file, dpi=500, height=3.2, width=5.4)

grit_cor_comparison_gg

