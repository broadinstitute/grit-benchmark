#!/usr/bin/env python
# coding: utf-8

# # Evaluating aggregated profiles
# 
# We evaluate the well-aggregated profiles by replicate reproducibility via:
# 1. cytominer-eval package
# 2. workflow adapted from *https://github.com/jump-cellpainting/workflow_demo_analysis/blob/main/analysis_Broad/0.percent_scores.ipynb*

# ## Calculate replicate reproducibility

# In[1]:


import os
import glob
import gzip
from pathlib import Path

import random
import numpy as np
import pandas as pd
from scipy.special import softmax 

import matplotlib.pyplot as plt
import seaborn as sns
import plotnine as gg

from pycytominer.cyto_utils import infer_cp_features
from cytominer_eval import evaluate
from scripts.utils import calculate_weighted_agg
import scripts.eval_utils as utils


# In[2]:


input_folder = 'data/aggregated-profiles/'
results_folder = 'results/'

df_list = [pd.read_csv(file, sep='\t') for file in glob.glob(input_folder+'*.tsv')]
concat_df = pd.concat(df_list, axis='rows', sort=False)
print(concat_df.shape)


# In[3]:


get_ipython().run_cell_magic('time', '', '\n# save the results\n# replicate correlation at the guide-level\nf_guide_cor_list = []\nf_guide_df_list = []\n# replicate correlation at the gene-level\nf_gene_cor_list = []\nf_gene_df_list = []\n# method list\nmethod_list = []\n# cell_line list\ncell_list = []\n\nfor method in concat_df.Metadata_agg_method.unique():\n    for cell_line in concat_df.Metadata_cell_line.unique():\n        print(f"method is {method}, cell line is {cell_line}")\n        sub_df = concat_df.query("Metadata_agg_method == @method & Metadata_cell_line == @cell_line")\n\n        # guide-level replicate correlation\n        f_guide_cor_val, f_guide_cor_df = evaluate(\n            profiles = sub_df,\n            features = infer_cp_features(sub_df),\n            meta_features = infer_cp_features(sub_df, metadata=True),\n            replicate_groups=[\'Metadata_pert_name\', \'Metadata_gene_name\', \'Metadata_cell_line\', \'Metadata_agg_method\'], \n            operation = \'replicate_reproducibility\',\n            similarity_metric = \'pearson\',\n        #         replicate_reproducibility_quantile = 0.99,\n            replicate_reproducibility_return_median_cor = True\n                                  )\n        f_guide_cor_list.append(f_guide_cor_val)\n        f_guide_df_list.append(f_guide_cor_df)    \n    \n        # gene-level replicate correlation\n        f_gene_cor_val, f_gene_cor_df = evaluate(\n            profiles = sub_df,\n            features = infer_cp_features(sub_df),\n            meta_features = infer_cp_features(sub_df, metadata=True),\n            replicate_groups=[\'Metadata_gene_name\', \'Metadata_cell_line\', \'Metadata_agg_method\'],\n            operation = \'replicate_reproducibility\',\n            similarity_metric = \'pearson\',\n        #         replicate_reproducibility_quantile= 0.99,\n            replicate_reproducibility_return_median_cor = True\n                                  )\n        f_gene_cor_list.append(f_gene_cor_val)\n        f_gene_df_list.append(f_gene_cor_df)\n        \n        method_list.append(method)\n        cell_list.append(cell_line)\n\nguide_rep_cors = pd.concat(f_guide_df_list)\ngene_rep_cors = pd.concat(f_gene_df_list)')


# In[4]:


pd.DataFrame(data={'Cell Line': cell_list, 'Aggregation Method': method_list, 
                   'Guide Replicate Reproducibility': f_guide_cor_list, 'Gene Replicate Reproducibility':f_gene_cor_list}).round(4)


# #### quick visualization of well-aggregated profiles

# In[5]:


box = (
    gg.ggplot(guide_rep_cors.query("~Metadata_agg_method.str.contains('ctrls')"), gg.aes(fill="Metadata_agg_method", x="Metadata_agg_method", y="similarity_metric"))
    + gg.geom_boxplot()
    + gg.facet_wrap('Metadata_cell_line')
    + gg.labs(title = 'Guide Replicate Correlation in Profiles')
    + gg.labs(x="Aggregation Method", y="Replicate Pearson Correlation", fill="Aggregation Method")
)
display(box)

jitter = (
    gg.ggplot(guide_rep_cors.query("~Metadata_agg_method.str.contains('ctrls')"), gg.aes(fill="Metadata_agg_method", x="Metadata_agg_method", y="similarity_metric"))
    + gg.geom_jitter()
    + gg.facet_wrap('Metadata_cell_line')
    + gg.labs(title = 'Guide Replicate Correlation in Profiles')
    + gg.labs(x="Aggregation Method", y="Replicate Pearson Correlation", fill="Aggregation Method")
)
display(jitter)


# In[6]:


box = (
    gg.ggplot(gene_rep_cors.query("~Metadata_agg_method.str.contains('ctrls')"), gg.aes(fill="Metadata_agg_method", x="Metadata_agg_method", y="similarity_metric"))
    + gg.geom_boxplot()
    + gg.facet_wrap('Metadata_cell_line')
    + gg.labs(title = 'Gene Replicate Correlation in Profiles')
    + gg.labs(x="Aggregation Method", y="Replicate Pearson Correlation", fill="Aggregation Method")
)
display(box)

jitter = (
    gg.ggplot(gene_rep_cors.query("~Metadata_agg_method.str.contains('ctrls')"), gg.aes(fill="Metadata_agg_method", x="Metadata_agg_method", y="similarity_metric"))
    + gg.geom_jitter()
    + gg.facet_wrap('Metadata_cell_line')
    + gg.labs(title = 'Gene Replicate Correlation in Profiles')
    + gg.labs(x="Aggregation Method", y="Replicate Pearson Correlation", fill="Aggregation Method")
)
display(jitter)


# #### Save outputs

# In[7]:


guide_rep_cors.to_csv(Path(results_folder + "guide_replicate_correlations.tsv"), index=False, sep='\t')
gene_rep_cors.to_csv(Path(results_folder + "gene_replicate_correlations.tsv"), index=False, sep='\t')


# ## Calculate percent replicating and percent matching
# *reference: https://github.com/jump-cellpainting/workflow_demo_analysis/blob/main/analysis_Broad/0.percent_scores.ipynb*

# In[8]:


random.seed(9000)

n_samples = 1000 # Number of points to sample from the null distribution.
n_replicates = 6 # Number of replicates of each sample.

replicate_grouping_feature = 'Metadata_pert_name' 
class_grouping_feature = 'Metadata_gene_name' 


# ### Percent Replicating (previously percent strong)
# 
# A "strong signature" is determined by :
# 1. high pairwise (pearson) correlation w/ each other
# 2. median replicate correlation is "significantly different" (aka perturbation with median replicate correlation greater than 95th percentile) from null distribution (pairwise correlations of non-replicates)

# In[9]:


get_ipython().run_cell_magic('time', '', '\ncorr_replicating_df = pd.DataFrame()\nreplicate_corr_list = []\ncell_list = []\nagg_list = []\n\nfor cell_line in concat_df.Metadata_cell_line.unique(): #[\'ES2\']:\n    for agg_method in concat_df.Metadata_agg_method.unique(): #[\'median\']: \n        print(f\'cell line is \\t{cell_line} \\naggregation method is \\t{agg_method}\')\n        # Subset for experiment\n        sub_df = concat_df.query("Metadata_cell_line == @cell_line & Metadata_agg_method == @agg_method").dropna(axis=\'columns\')\n\n        print(f\'The combined experiment has {sub_df.shape[0]} wells and {sub_df.shape[1]} features\')\n\n        # Remove empty wells and negative controls\n        compound_df = sub_df.query("Metadata_gene_name != \'Chr2\'")\n\n        # Calculate list of correlations for replicates\n        replicate_corr_df = utils.corr_between_replicates(compound_df,\n                                                    replicate_grouping_feature)\n\n        # Calculate list of correlations for non-replicates \n        null_replicate_corr = list(utils.corr_between_non_replicates(compound_df,\n                                                               n_samples,\n                                                               n_replicates,\n                                                               replicate_grouping_feature))\n\n        # Calculate the percent score (% of replicates above 95th %ile of null)\n        prop_95_replicating, value_95_replicating = utils.percent_score(null_replicate_corr,\n                                                                  replicate_corr_df.correlation,\n                                                                  how=\'right\')\n\n        corr_replicating_df = corr_replicating_df.append({\'Experiment\':f\'{cell_line} {agg_method} plates\',\n                                                          \'Replicate_corr\':list(replicate_corr_df.correlation),\n                                                          \'Null_corr\':null_replicate_corr,\n                                                          \'Percent_Replicating\':\'%.3f\'%prop_95_replicating,\n                                                          \'Value_95\':value_95_replicating}, ignore_index=True)\n\n\n        cell_list.append(cell_line)\n        agg_list.append(agg_method)')


# In[10]:


display(corr_replicating_df)


# ### Percent Matching (previously percent recall)
# 
# Percent Matching is similar to Percent Replicating, but instead of pairwise correlations between replicates, we calculate the pairwise correlations between a pair of compounds that have the same MOA annotation or have been known to target the same gene. 
# 
# In this case, the null distribution is constructed from the pairwise correlations of compounds that belong to different MOA classes or target different genes.
# 

# In[11]:


get_ipython().run_cell_magic('time', '', '\ncorr_matching_df = pd.DataFrame()\ncell_list = []\nagg_list = []\n\nfor cell_line in concat_df.Metadata_cell_line.unique(): #[\'ES2\']:\n    for agg_method in concat_df.Metadata_agg_method.unique(): #[\'median\']:\n        print(f\'cell line is \\t{cell_line} \\naggregation method is \\t{agg_method}\')\n        # Subset for experiment\n        sub_df = concat_df.query("Metadata_cell_line == @cell_line & Metadata_agg_method == @agg_method").dropna(axis=\'columns\')\n\n        print(f\'The combined experiment has {sub_df.shape[0]} wells and {sub_df.shape[1]} features\')\n\n        # Remove empty wells and negative controls\n        compound_df = sub_df.query("Metadata_gene_name != \'Chr2\'")\n        \n        # Calculate list of correlations for pairs \n        matching_corr = list(utils.corr_between_perturbation_pairs(compound_df,\n                                                                   class_grouping_feature,\n                                                                   replicate_grouping_feature))\n\n\n        # Calculate list of corremations for non-pairs\n        null_matching = list(utils.corr_between_perturbation_non_pairs(compound_df,\n                                                                       n_samples,\n                                                                       class_grouping_feature,\n                                                                       replicate_grouping_feature))\n        \n        \n        # Calculate the percent score (% of replicates above 95th %ile of null)\n        prop_95_matching, value_95_matching = utils.percent_score(null_matching,\n                                                                  matching_corr,\n                                                                  how=\'right\')\n\n\n        corr_matching_df = corr_matching_df.append({\'Experiment\':f\'{cell_line} {agg_method} plates\',\n                                                    \'Matching_corr\':matching_corr,\n                                                    \'Null_Matching\':null_matching,\n                                                    \'Percent_Matching\':\'%.3f\'%prop_95_matching,\n                                                    \'Value_95\':value_95_matching}, ignore_index=True)\n\n        cell_list.append(cell_line)\n        agg_list.append(agg_method)')


# In[12]:


display(corr_matching_df)


# In[13]:


guide_sub = corr_replicating_df.drop(columns=['Null_corr', 'Replicate_corr'], axis='columns', inplace=False)
gene_sub = corr_matching_df.drop(columns=['Matching_corr', 'Null_Matching', ], axis='columns', inplace=False)
merged_df = pd.merge(guide_sub, gene_sub, how='inner', on='Experiment')
merged_df[['Cell Line', 'Aggregation Method']] = merged_df.Experiment.str.split(' ', expand=True).loc[:,0:1]
merged_df[['Cell Line', 'Aggregation Method', 
           'Percent_Replicating', 'Value_95_x', 'Percent_Matching', 'Value_95_y']].sort_values("Aggregation Method")


# #### Save outputs

# In[14]:


corr_replicating_df.to_csv(Path(results_folder + "percent_replicating.tsv"), index=False, sep='\t')
corr_matching_df.to_csv(Path(results_folder + "percent_matching.tsv"), index=False, sep='\t')

