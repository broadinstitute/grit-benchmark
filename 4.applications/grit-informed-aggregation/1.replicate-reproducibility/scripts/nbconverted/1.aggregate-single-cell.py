#!/usr/bin/env python
# coding: utf-8

# # Aggregate single cells into bulk profiles using grit information
#
# We use normalized, **single-cell profiles** from the Cell Health experiment and **single-cell grit scores** (calculated with respect to Chr2 cutting controls and using normalized, feature selected profiles).
#
# We only used the same cell-painting feature columns of Cell Health data in the **Grit Benchmarking project** to ensure consistency across analyses.
#
# Here we aggregate profiles from the Cell Health experiments using several "grit-informed" methods:
# 1. Standard median aggregation
# 2. Weighted mean, weighting by raw grit scores
# 3. Weighted mean, weighting by softmax-transformed grit scores
# 4. Weighted mean, weighting by grit scores that are clipped to 0. This assigns a minimum grit score of 0 to any cell with grit scores < 0.
# 5. Weighted mean, weighting by grit scores that are clipped to 0 then softmax-transformed.
#
# Methods 3-5 were exploratory and did not yield dramatically improved results ("improved" assessed by replicate reproducibility measures). Therefore, they are commented out in this notebook to save on runtime and compute.
#
# Note: depending on AWS instance size, some cell lines were aggregated with a python script of the same name in `scripts/`.

# In[1]:


import os
import glob
import gzip
from pathlib import Path
from datetime import datetime

import numpy as np
import pandas as pd
from scipy.special import softmax

from pycytominer import aggregate, get_na_columns
from pycytominer.cyto_utils import infer_cp_features
from cytominer_eval import evaluate
from scripts.utils import calculate_weighted_agg


# In[2]:


def merge_metadata(cell_line, level3_profile):
    # load metadata file from Cell Health  data
    commit = "07e4b40c39dd27084be36fbef4d64c5654b2960f"
    base_url = f"https://github.com/broadinstitute/cell-health/raw/{commit}"
    url = f"{base_url}/1.generate-profiles/data/metadata/platemap/DEPENDENCIES1_{cell_line}.csv"
    platemap = pd.read_csv(url, sep=",")
    platemap.columns = ["Metadata_" + str(x) for x in platemap.columns]

    # merge with the aggregated files
    meta_df = pd.merge(
        level3_profile,
        platemap,
        left_on="Metadata_Well",
        right_on="Metadata_well_position",
    )
    # # reorder columns for metadata to be in front
    meta_df = meta_df[
        sorted(meta_df, key=lambda x: x not in meta_df.filter(like="Metadata").columns)
    ]

    return meta_df


# ## Load level 2 data

# In[3]:


plate_dict = {
    "ES2": ["SQ00014613", "SQ00014614", "SQ00014615"],
    "A549": ["SQ00014610", "SQ00014611", "SQ00014612"],
    "HCC44": ["SQ00014616", "SQ00014617", "SQ00014618"],
}


# In[4]:


# take the same columns as Cell Health data in rest of Grit Benchmark project
commit = "2916770cc9cc9e75b693348b683aa398987fb9f9"
base_url = f"https://github.com/broadinstitute/grit-benchmark/raw/{commit}"
url = f"{base_url}/1.calculate-metrics/cell-health/data/cell_health_merged_feature_select.csv.gz"

df = pd.read_csv(url, sep=",")
print(df.shape)
df.head(2)

cols_to_keep = infer_cp_features(df)


# In[ ]:


get_ipython().run_cell_magic(
    "time",
    "",
    "results_folder = 'data/aggregated-profiles/'\n    \nfor cell_line in ['ES2', 'HCC44', 'A549']: \n    ####### read in single-cell grit data #######\n    start_merge = datetime.now()\n    grit_folder = '../../../1.calculate-metrics/cell-health/results/'\n    grit_files = glob.glob(grit_folder+'*single_cell_grit*.tsv.gz')\n\n    scgrit_df = []\n    for file in grit_files:\n        plate_name=file.split('/')[-1].split('_')[-2]\n        if plate_name in plate_dict[cell_line]:\n            print(f\"adding scrgrit of {plate_name} to list of {cell_line}\")\n            scgrit_plate = pd.read_csv(file, sep='\\t').assign(plate=plate_name, cell_line = cell_line)\n            print(scgrit_plate.shape)\n            scgrit_df.append(scgrit_plate)\n    scgrit_df = pd.concat(scgrit_df)\n    scgrit_df['cell_identity'] = scgrit_df.perturbation.str.split(\"_\", expand=True)[1].astype(int)\n    scgrit_df.columns = ['Metadata_'  + str(col) for col in scgrit_df.columns]\n    print(f\"total shape of of scgrit_df for {cell_line} is: {scgrit_df.shape}\")\n    \n    ####### read in single-cell cell painting profiles #######\n    profile_folder = '../../../0.download-data/data/cell_health/normalized/' \n    profile_files = glob.glob(profile_folder+'*normalized.csv.gz')\n\n    scprofiles_df = []\n    for file in profile_files:\n        plate_name=file.split('/')[-1].split('_')[0]\n        if plate_name in plate_dict[cell_line]:\n            print(f\"adding scprofiles of {plate_name} to list of {cell_line}\")\n            scprofile_plate = (pd.read_csv(file, sep=',', low_memory=False)\n                               .reset_index()\n                               .rename({'index':'Metadata_cell_identity'}, axis='columns')\n                              ).assign(cell_line = cell_line)\n            plate_cols = infer_cp_features(scprofile_plate)\n            drop_cols = [x for x in plate_cols if x not in cols_to_keep]\n            scprofile_plate.drop(columns = drop_cols, inplace=True)\n            scprofiles_df.append(scprofile_plate)\n    scprofiles_df = pd.concat(scprofiles_df, sort=False)\n    print(f'total shape of scprofiles_df for {cell_line} is: {scprofiles_df.shape}')\n\n    ####### merge scgrit scores + cell painting profiles #######\n    scprofiles_df = (pd.merge(scprofiles_df, scgrit_df, \n         left_on=['Metadata_cell_identity', 'Metadata_Plate', 'Metadata_pert_name'], \n                 right_on=['Metadata_cell_identity', 'Metadata_plate', 'Metadata_group'])\n        )\n    del scgrit_df\n    print(f\"total shape of sc_df for {cell_line} is: {scprofiles_df.shape}\")\n    # remove columns with any NA entries\n    na_cols_to_drop = get_na_columns(scprofiles_df, cutoff=0)\n    print(f\"Dropping {len(na_cols_to_drop)} columns because of missing data\")\n    scprofiles_df = scprofiles_df.drop(na_cols_to_drop, axis=\"columns\")\n    print(f\"FINAL shape of merged data {scprofiles_df.shape}\")\n\n    print(f\"TOTAL TIME constructing merged df for cell_line {cell_line} : {str(datetime.now()-start_merge)}\")\n    \n    \n    ###### standard median aggregation ######\n    start_agg = datetime.now()\n    agg_df = aggregate(\n        population_df = scprofiles_df,\n        strata = [\"Metadata_Plate\", \"Metadata_Well\"],\n        features = \"infer\",\n        operation =\"median\"\n    ).assign(Metadata_agg_method = 'median', cell_line = cell_line)\n    agg_meta_df = merge_metadata(cell_line, agg_df)\n    display(agg_meta_df.head())\n    # writing data\n    agg_meta_df.to_csv(Path(results_folder + cell_line + \"_median.tsv\"), index=False, sep='\\t')\n    \n    ###### grit-informed aggregation methods ######\n    ### raw grit as weights ###\n    agg_df = (calculate_weighted_agg(\n        population_df = scprofiles_df,\n        columns = ['Metadata_Plate', 'Metadata_Well'],\n        features = 'infer',\n        transform = 'weighted_grit', weight = 'Metadata_grit')\n                    ).assign(Metadata_agg_method = 'weighted', cell_line = cell_line)\n    agg_meta_df = merge_metadata(cell_line, agg_df)\n    display(agg_meta_df.head())\n    # writing data\n    agg_meta_df.to_csv(Path(results_folder + cell_line + \"_weighted.tsv\"), index=False, sep='\\t')\n    \n#     ### grit that is softmax-transformed as weights ###\n#     agg_df = (calculate_weighted_agg(\n#         population_df = scprofiles_df,\n#         columns = ['Metadata_Plate', 'Metadata_Well'],\n#         features = 'infer',\n#         transform = 'softmax_grit', weight = 'Metadata_grit')\n#                    ).assign(Metadata_agg_method = 'softmax', cell_line = cell_line)\n#     agg_meta_df = merge_metadata(cell_line, agg_df)\n#     # writing data\n#     agg_meta_df.to_csv(Path(results_folder + cell_line + \"_softmax.tsv\"), index=False, sep='\\t')\n    \n#     ### grit clipped to 0 (as lowest values), as weights ###\n#     agg_df = (calculate_weighted_agg(\n#         population_df = scprofiles_df,\n#         columns = ['Metadata_Plate', 'Metadata_Well'], \n#         features = 'infer',\n#         transform = 'weighted_grit', weight='Metadata_clipped_grit', lower_threshold=0)\n#                        ).assign(Metadata_agg_method = 'clipped0_weighted', cell_line = cell_line)\n#     agg_meta_df = merge_metadata(cell_line, agg_df)\n#     # writing data\n#     agg_meta_df.to_csv(Path(results_folder + cell_line + \"_clipped0_weighted.tsv\"), index=False, sep='\\t')\n    \n#     ### grit clipped to 0 (as lowest values), then softmax-transfored, as weights ###\n#     agg_df = (calculate_weighted_agg(\n#         population_df = scprofiles_df, \n#         columns = ['Metadata_Plate', 'Metadata_Well'], \n#         features = 'infer',\n#         transform = 'softmax_grit', weight='Metadata_clipped_grit', lower_threshold=0)\n#                        ).assign(Metadata_agg_method = 'clipped0_softmax')\n#     agg_meta_df = merge_metadata(cell_line, agg_df)\n#     # writing data\n#     agg_meta_df.to_csv(Path(results_folder + cell_line + \"_clipped0_softmax.tsv\"), index=False, sep='\\t')\n    \n\n    print(f\"TOTAL TIME performing aggregation for cell_line {cell_line} : {str(datetime.now()-start_agg)}\")",
)
