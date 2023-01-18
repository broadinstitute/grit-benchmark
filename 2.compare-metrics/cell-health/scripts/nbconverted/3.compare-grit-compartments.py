#!/usr/bin/env python
# coding: utf-8

# # Visualize grit calculated per compartment
# 
# Previously, we calculated grit using the following subsets:
# 
# * Per compartment (Cells, Cytoplasm, Nuclei)
# * Per compartment feature group (Cells x AreaShape, Cells x Correlation, Nuclei x Texture, etc.)
# * Per channel (DNA, RNA, AGP, ER, Mito)
#     * Use all features that include any information from one of these channels
# * No feature subsetting
# 
# Here, we visualize these results.

# In[1]:


import pathlib
import pandas as pd
import plotnine as gg


# In[2]:


output_dir = pathlib.Path("figures/compartment_drop")

cell_health_dir = pathlib.Path("../../1.calculate-metrics/cell-health/results")
grit_file = pathlib.Path(f"{cell_health_dir}/cell_health_grit.tsv")
compartment_grit_file = pathlib.Path(f"{cell_health_dir}/cell_health_grit_compartments.tsv.gz")


# In[3]:


cell_line_colors = {
  "A549": "#861613",
  "ES2": "#1CADA8",
  "HCC44": "#2A364D"
}

compartment_drop_theme = gg.theme(
    strip_background=gg.element_rect(color="black", fill="#fdfff4"),
    strip_text_x=gg.element_text(size=6),
    axis_text=gg.element_text(size=7),
    axis_title=gg.element_text(size=8),
    legend_title=gg.element_text(size=6),
    legend_text=gg.element_text(size=5),
    panel_grid=gg.element_line(size=0.35)
)

grit_results = (
    pd.read_csv(grit_file, sep="\t")
    .assign(compartment="all", feature_group="all", channel="all")
    .query("barcode_control == 'cutting_control'")
    .query("grit_replicate_summary_method == 'mean'")
    .query("cor_method == 'pearson'")
)

grit_results = grit_results.assign(num_features=grit_results.shape[0])
compartment_grit_results = pd.read_csv(compartment_grit_file, sep="\t")


# ## 1) Per compartment grit

# In[4]:


per_compartment_df = (
    compartment_grit_results
    .query("compartment != 'all'")
    .query("feature_group == 'all'")
    .query("channel == 'all'")
)

print(
    pd.crosstab(
        per_compartment_df.num_features,
        per_compartment_df.compartment
    )
)

print(per_compartment_df.shape)
per_compartment_df.head()


# In[5]:


per_compartment_df = (
    per_compartment_df
    .pivot(
        index=["perturbation", "group", "cell_line", "compartment"],
        values="grit",
        columns="dropped_or_exclusive"
    )
    .reset_index()
)

per_compartment_df = (
    per_compartment_df
    .assign(channel_signal = per_compartment_df.exclusive - per_compartment_df.dropped)
    .sort_values(by="channel_signal", ascending=False)
    .reset_index(drop=True)
)

print(per_compartment_df.shape)
per_compartment_df.head(20)


# In[6]:


compartment_drop_gg = (
    gg.ggplot(per_compartment_df, gg.aes(x="dropped", y="exclusive"))
    + gg.geom_point(gg.aes(fill="cell_line"), alpha=0.4, size=1, stroke=0.1)
    + gg.facet_grid("~compartment")
    + gg.xlab("Compartment dropped (grit)")
    + gg.ylab("Compartment only (grit)")
    + gg.coord_fixed()
    + gg.geom_abline(intercept=0, slope=1, size=0.15, color="red", linetype="dashed")
    + gg.scale_fill_manual(name="Cell Line", values=cell_line_colors)
    + gg.theme_bw()
    + compartment_drop_theme
)

output_file = pathlib.Path(f"{output_dir}/cell_health_grit_compartment_drop.png")
compartment_drop_gg.save(output_file, dpi=500, height=3, width=6)

compartment_drop_gg


# ## 2) Per compartment feature group

# In[7]:


per_featuregroup_df = (
    compartment_grit_results
    .query("compartment != 'all'")
    .query("feature_group != 'all'")
    .query("channel == 'all'")
    .query("feature_group != 'Location'")
    .pivot(
        index=["perturbation", "group", "cell_line", "channel", "feature_group", "compartment"],
        values="grit",
        columns="dropped_or_exclusive")
    .reset_index()
)

per_featuregroup_df = (
    per_featuregroup_df
    .assign(channel_signal = per_featuregroup_df.exclusive - per_featuregroup_df.dropped)
    .sort_values(by="channel_signal", ascending=False)
    .reset_index(drop=True)
)

print(per_featuregroup_df.shape)
per_featuregroup_df.head(20)


# In[8]:


feature_group_drop_gg = (
    gg.ggplot(per_featuregroup_df, gg.aes(x="dropped", y="exclusive"))
    + gg.geom_point(gg.aes(fill="cell_line"), alpha=0.4, size=0.5, stroke=0.1)
    + gg.facet_grid("compartment~feature_group")
    + gg.xlab("Feature group dropped (grit)")
    + gg.ylab("Feature group only (grit)")
    + gg.coord_fixed()
    + gg.geom_abline(intercept=0, slope=1, size=0.15, color="red", linetype="dashed")
    + gg.scale_fill_manual(name="Cell Line", values=cell_line_colors)
    + gg.theme_bw()
    + compartment_drop_theme
)

output_file = pathlib.Path(f"{output_dir}/cell_health_grit_featuregroup_drop.png")
feature_group_drop_gg.save(output_file, dpi=500, height=6, width=8)

feature_group_drop_gg


# ## 3) Per channel grit

# In[9]:


per_channel_df = (
    compartment_grit_results
    .query("compartment == 'all'")
    .query("feature_group == 'all'")
    .query("channel != 'all'")
    .pivot(index=["perturbation", "group", "cell_line", "channel"], values="grit", columns="dropped_or_exclusive")
    .reset_index()
)

per_channel_df = (
    per_channel_df
    .assign(channel_signal = per_channel_df.exclusive - per_channel_df.dropped)
    .sort_values(by="channel_signal", ascending=False)
    .reset_index(drop=True)
)

print(per_channel_df.shape)
per_channel_df.head(10)


# In[10]:


channel_drop_gg = (
    gg.ggplot(per_channel_df, gg.aes(x="dropped", y="exclusive"))
    + gg.geom_point(gg.aes(fill="cell_line"), alpha=0.4, size=0.5, stroke=0.1)
    + gg.facet_grid("~channel")
    + gg.xlab("Channel dropped (grit)")
    + gg.ylab("Channel only (grit)")
    + gg.geom_abline(intercept=0, slope=1, size=0.15, color="red", linetype="dashed")
    + gg.scale_fill_manual(name="Cell Line", values=cell_line_colors)
    + gg.coord_fixed()
    + gg.theme_bw()
    + compartment_drop_theme
)

output_file = pathlib.Path(f"{output_dir}/cell_health_grit_channel_drop.png")
channel_drop_gg.save(output_file, dpi=500, height=3, width=6)

channel_drop_gg

