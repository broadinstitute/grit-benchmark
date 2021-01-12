#!/usr/bin/env python
# coding: utf-8

# In[44]:


import pathlib
import pandas as pd
import plotnine as gg


# In[45]:


output_dir = pathlib.Path("figures")

cell_health_dir = pathlib.Path("../../1.calculate-metrics/cell-health/results")
grit_file = pathlib.Path(f"{cell_health_dir}/cell_health_grit.tsv")
compartment_grit_file = pathlib.Path(f"{cell_health_dir}/cell_health_grit_compartments.tsv.gz")


# In[47]:


grit_results = (
    pd.read_csv(grit_file, sep="\t")
    .assign(compartment="all", feature_group="all")
    .query("barcode_control == 'cutting_control'")
)

grit_results = grit_results.assign(num_features=grit_results.shape[0])

compartment_grit_results = pd.read_csv(compartment_grit_file, sep="\t")


# In[49]:


grit_results = (
    pd.concat(
        [
            grit_results,
            compartment_grit_results
        ],
        axis="rows"
    )
    .dropna()
    .reset_index(drop=True)
)

print(grit_results.shape)
grit_results.head()


# In[52]:


(
    gg.ggplot(grit_results.query("group == 'ITGAV'"), gg.aes(x="compartment+feature_group", y="grit"))
    + gg.geom_boxplot()
    + gg.facet_wrap("~channel", nrow=4, scales="free_y")
    + gg.theme_bw()
    + gg.theme(axis_text_x=gg.element_text(angle=90))
    + gg.coord_flip()
)


# In[ ]:





# In[ ]:





# In[16]:


(
    gg.ggplot(grit_results, gg.aes(x="compartment", y="grit"))
    + gg.geom_boxplot()
    + gg.facet_wrap("~feature_group")
)


# In[ ]:




