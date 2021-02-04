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

grit_df = grit_df.assign(
    replicate_over_cor=grit_df.median_replicate_correlation / grit_df.median_control_correlation,
    replicate_over_oneminuscor=grit_df.median_replicate_correlation / (1 - grit_df.median_control_correlation)
)
grit_df.replicate_over_cor = grit_df.replicate_over_cor.abs()

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


# In[6]:


replicate_gg = (
    gg.ggplot(grit_df, gg.aes(x="grit", y="replicate_over_cor"))
    + gg.geom_point(gg.aes(fill="median_control_correlation", size="median_control_correlation"), alpha=0.7)
    + gg.scale_fill_distiller(name="Median\ncontrol\ncorrelation\n", type='div')
    + gg.scale_size_continuous(guide=False, range=[0.1, 4], breaks=[0, 0.1, 2, 3, 5])
    + gg.facet_wrap("~cell_line", nrow=3)
    + gg.theme_bw()
    + gg.geom_hline(yintercept=0, linetype="dashed", color="red")
    + gg.geom_vline(xintercept=0, linetype="dashed", color="red")
    + gg.xlab("Grit")
    + gg.ylab("Replicate correlation / control correlation")
    + replicate_reproducibility_theme
)

replicate_gg


# In[7]:


replicate_gg = (
    gg.ggplot(grit_df.query("replicate_over_cor < 100"), gg.aes(x="grit", y="replicate_over_cor"))
    + gg.geom_point(gg.aes(fill="median_control_correlation", size="median_control_correlation"), alpha=0.7)
    + gg.scale_fill_distiller(name="Median\ncontrol\ncorrelation\n", type='div')
    + gg.scale_size_continuous(guide=False, range=[0.1, 4], breaks=[0, 0.1, 2, 3, 5])
    + gg.facet_wrap("~cell_line", nrow=3)
    + gg.theme_bw()
    + gg.geom_hline(yintercept=0, linetype="dashed", color="red")
    + gg.geom_vline(xintercept=0, linetype="dashed", color="red")
    + gg.xlab("Grit")
    + gg.ylab("Replicate correlation / control correlation")
    + replicate_reproducibility_theme
)

replicate_gg


# In[8]:


replicate_gg = (
    gg.ggplot(grit_df.query("replicate_over_cor < 100"), gg.aes(x="grit", y="replicate_over_oneminuscor"))
    + gg.geom_point(gg.aes(fill="median_control_correlation", size="median_control_correlation"), alpha=0.7)
    + gg.scale_fill_distiller(name="Median\ncontrol\ncorrelation\n", type='div')
    + gg.scale_size_continuous(guide=False, range=[0.1, 4], breaks=[0, 0.1, 2, 3, 5])
    + gg.facet_wrap("~cell_line", nrow=3)
    + gg.theme_bw()
    + gg.geom_hline(yintercept=0, linetype="dashed", color="red")
    + gg.geom_vline(xintercept=0, linetype="dashed", color="red")
    + gg.xlab("Grit")
    + gg.ylab("Replicate correlation / (1 - control correlation)")
    + replicate_reproducibility_theme
)

replicate_gg


# In[9]:


replicate_gg = (
    gg.ggplot(grit_df.query("replicate_over_cor < 100"), gg.aes(x="grit", y="replicate_over_oneminuscor"))
    + gg.geom_point(gg.aes(fill="median_replicate_correlation", size="median_replicate_correlation"), alpha=0.7)
    + gg.scale_fill_distiller(name="Median\nreplicate\ncorrelation\n", type='div')
    + gg.scale_size_continuous(guide=False, range=[0.1, 4], breaks=[0, 0.1, 2, 3, 5])
    + gg.facet_wrap("~cell_line", nrow=3)
    + gg.theme_bw()
    + gg.geom_hline(yintercept=0, linetype="dashed", color="red")
    + gg.geom_vline(xintercept=0, linetype="dashed", color="red")
    + gg.xlab("Grit")
    + gg.ylab("Replicate correlation / (1 - control correlation)")
    + replicate_reproducibility_theme
)

replicate_gg

