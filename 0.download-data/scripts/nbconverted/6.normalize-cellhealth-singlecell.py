#!/usr/bin/env python
# coding: utf-8

# ## Normalizing Cell Health single cell profiles and output normalized single cell data

# In[1]:


import random
import pathlib
import pandas as pd

from pycytominer import normalize
from pycytominer.cyto_utils import cells, output


# In[2]:


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


# In[3]:


compression_options = {"method": "gzip", "mtime": 1}


# In[4]:


# Set output directory where normalized single cell files will live
output_dir = pathlib.Path("data/cell_health/normalized/")
output_dir.mkdir(exist_ok=True, parents=True)


# In[5]:


# Load Metadata on plate info
metadata_file = "https://github.com/broadinstitute/cell-health/blob/cd91bd0daacef2b5ea25dcceb62482bb664d9de1/1.generate-profiles/data/metadata/platemap/DEPENDENCIES1_ES2.csv?raw=True"
metadata_df = pd.read_csv(metadata_file)

metadata_df.columns = [f"Metadata_{x}" for x in metadata_df.columns]

print(metadata_df.shape)
metadata_df.head()


# In[6]:


for plate in plates:
    # Set file names
    sql_file = f"sqlite:///data/cell_health/{plate}.sqlite"
    output_file = pathlib.Path(output_dir, f"{plate}_normalized.csv.gz")
    
    # Set console output
    print(f"Now processing... {output_file}")

    # Initiate single cell class
    sc = cells.SingleCells(
        file_or_conn=sql_file,
        strata=["Image_Metadata_Plate", "Image_Metadata_Well"],
    )
    
    # Merge single cells
    sc_df = sc.merge_single_cells()
    
    # Normalize data
    sc_df = normalize(
        profiles=sc_df,
        method="standardize"
    )
    
    # Merge well and plate metadata
    sc_df = (
        sc.image_df.merge(
            metadata_df,
            left_on="Image_Metadata_Well",
            right_on="Metadata_well_position",
            how="left"
        ).merge(
            sc_df,
            left_on=["TableNumber", "ImageNumber"],
            right_on=["Metadata_TableNumber", "Metadata_ImageNumber"],
            how="right"
        )
        .drop(
            [
                "TableNumber",
                "ImageNumber",
                "Metadata_WellRow",
                "Metadata_WellCol",
                "Metadata_well_position"
            ], axis="columns"
        )
        .rename(
            {
                "Image_Metadata_Plate": "Metadata_Plate",
                "Image_Metadata_Well": "Metadata_Well"
            }, axis="columns"
        )
    )

    # Print data shape
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
    print("\n")

