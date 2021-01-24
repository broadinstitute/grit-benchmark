#!/usr/bin/env python
# coding: utf-8

# # Visualize grit as compared to replicate reproducibility
# 
# Grit combines the concepts of replicate reproducibility and difference from controls.
# 
# Here, observe how grit handles the tradeoff between calculating each independently.

# In[1]:


import pathlib
import pandas as pd
import plotnine as gg


# In[2]:


output_dir = pathlib.Path("figures/replicate_reproducibility")

cell_health_dir = pathlib.Path("../../1.calculate-metrics/cell-health/results")
grit_file = pathlib.Path(f"{cell_health_dir}/cell_health_grit.tsv")
reprod_file = pathlib.Path(f"{cell_health_dir}/cell_health_replicate_reproducibility.tsv")


# In[3]:


cell_line_colors = {
  "A549": "#861613",
  "ES2": "#1CADA8",
  "HCC44": "#2A364D"
}

replicate_reproducibility_theme = gg.theme(
    strip_background=gg.element_rect(color="black", fill="#fdfff4"),
    strip_text_x=gg.element_text(size=6),
    axis_text=gg.element_text(size=7),
    axis_title=gg.element_text(size=8),
    legend_title=gg.element_text(size=6),
    legend_text=gg.element_text(size=5),
    legend_key_size=10,
    legend_key_width=10,
    legend_key_height=10,
    panel_grid=gg.element_line(size=0.35)
)


# In[4]:


# Load and process data
grit_df = pd.read_csv(grit_file, sep="\t").query("barcode_control == 'cutting_control'").query("cor_method == 'pearson'")
reprod_df = pd.read_csv(reprod_file, sep="\t")

grit_df = (
    grit_df.merge(
        reprod_df,
        on=["perturbation", "group", "cell_line"],
        how="inner"
    )
    .dropna()
)

print(grit_df.shape)
grit_df.head()


# In[5]:


replicate_gg = (
    gg.ggplot(grit_df, gg.aes(x="median_replicate_correlation", y="median_control_correlation"))
    + gg.geom_point(gg.aes(fill="grit", size="grit"), alpha=0.7)
    + gg.scale_fill_gradient(name="Grit", high="yellow", low="red")
    + gg.scale_size_continuous(guide=False, range=[0.1, 4], breaks=[0, 0.1, 2, 3, 5])
    + gg.facet_wrap("~cell_line", nrow=3)
    + gg.theme_bw()
    + gg.geom_hline(yintercept=0, linetype="dashed", color="red")
    + gg.geom_vline(xintercept=0, linetype="dashed", color="red")
    + gg.xlab("Median replicate correlation")
    + gg.ylab("Median control correlation")
    + replicate_reproducibility_theme
)

output_file = pathlib.Path(f"{output_dir}/cell_health_grit_replicate_reproduce.png")
replicate_gg.save(output_file, dpi=500, height=5, width=3)

replicate_gg

