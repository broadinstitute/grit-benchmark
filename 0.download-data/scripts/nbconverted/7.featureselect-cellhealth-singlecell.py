#!/usr/bin/env python
# coding: utf-8

# ## Applying feature selection to single cell data

# In[1]:


import random
import pathlib
import pandas as pd

from pycytominer import feature_select
from pycytominer.cyto_utils import output


# In[2]:


feature_select_operations = [
    "variance_threshold",
    "correlation_threshold",
    "drop_na_columns",
    "blocklist",
    "drop_outliers",
]

na_cutoff = 0

compression_options = {"method": "gzip", "mtime": 1}


# In[3]:


sc_dir = pathlib.Path("data/cell_health/normalized/")

plates = [
    "SQ00014610",
    "SQ00014611",
    "SQ00014612",
    "SQ00014613",
    "SQ00014614",
    "SQ00014615",
    "SQ00014616",
    "SQ00014617",
    "SQ00014618",
]

plate_files = {
    plate: pathlib.Path(f"{sc_dir}/{plate}_normalized.csv.gz") for plate in plates
}
plate_files


# In[4]:


for plate in plate_files:
    plate_file = plate_files[plate]
    output_file = pathlib.Path(f"{sc_dir}/{plate}_normalized_featureselected.csv.gz")

    # Set console output
    print(f"Now performing feature selection for... {plate_file}")
    sc_df = pd.read_csv(plate_file, low_memory=False)
    print("Before feature selection:")
    print(sc_df.shape)

    sc_df = feature_select(
        profiles=sc_df,
        operation=feature_select_operations,
        na_cutoff=na_cutoff,
    )

    print("After feature selection:")
    print(sc_df.shape)

    # Output file to disk
    output(
        df=sc_df,
        output_filename=output_file,
        sep=",",
        float_format="%.5f",
        compression_options=compression_options,
    )

    print("Done.")
    print("\n\n")
