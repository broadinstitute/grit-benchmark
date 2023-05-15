#!/usr/bin/env python
# coding: utf-8

# ## Reprocess Cell Health profiles
#
# Use a whole-plate normalization scheme instead of normalization by controls only.
#
# We will use the control normalization in downstream analyses, but we are interested in comparing the impact of normalization strategy on grit calculations.

# In[1]:


import pathlib
import pandas as pd

from pycytominer import normalize, feature_select


# In[2]:


def normalize_profile(
    plate, output_dir, commit="cd91bd0daacef2b5ea25dcceb62482bb664d9de1"
):
    link = f"https://github.com/broadinstitute/cell-health/raw/{commit}/1.generate-profiles/data/profiles/{plate}/{plate}_augmented.csv.gz"

    annotate_df = pd.read_csv(link)

    norm_file = pathlib.Path(f"{output_dir}/{plate}_wholeplate_normalized.csv.gz")
    feat_select_file = pathlib.Path(
        f"{output_dir}/{plate}_wholeplate_normalized_feature_selected.csv.gz"
    )

    normalize(
        profiles=annotate_df,
        features="infer",
        meta_features=meta_features,
        samples="all",
        method="mad_robustize",
        output_file=norm_file,
        compression_options={"method": "gzip", "mtime": 1},
    )


# In[3]:


# Define the plates
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

# Define metadata features
meta_features = [
    "Image_Metadata_Plate",
    "Image_Metadata_Well",
    "Metadata_WellRow",
    "Metadata_WellCol",
    "Metadata_gene_name",
    "Metadata_pert_name",
    "Metadata_broad_sample",
    "Metadata_cell_line",
]

output_dir = pathlib.Path("data/cell-health/profiles")

commit = "cd91bd0daacef2b5ea25dcceb62482bb664d9de1"


# In[4]:


for plate in plates:
    normalize_profile(plate, output_dir, commit)


# ## Now form a single merged dataset to perform feature selection

# In[5]:


# Load different normalized data
plate_files = [x for x in output_dir.iterdir() if "_normalized.csv.gz" in x.name]


# In[6]:


# Concatentate all plates
x_df = (
    pd.concat([pd.read_csv(x) for x in plate_files], sort=True)
    .rename(
        {
            "Image_Metadata_Plate": "Metadata_Plate",
            "Image_Metadata_Well": "Metadata_Well",
        },
        axis="columns",
    )
    .drop(["Metadata_broad_sample"], axis="columns")
)

# Realign metadata column names
x_metadata_cols = x_df.columns[x_df.columns.str.startswith("Metadata")]
x_metadata_df = x_df.loc[:, x_metadata_cols]

x_df = x_df.drop(x_metadata_cols, axis="columns")
x_df = pd.concat([x_metadata_df, x_df], axis="columns")

print(x_df.shape)
x_df.head()


# In[7]:


# Perform feature selection
feature_select_ops = [
    "variance_threshold",
    "correlation_threshold",
    "drop_na_columns",
    "blocklist",
    "drop_outliers",
]

x_df = feature_select(profiles=x_df, operation=feature_select_ops, na_cutoff=0)

print(x_df.shape)
x_df.head(2)


# In[8]:


# Also drop Costes features
costes_cols_to_drop = [x for x in x_df.columns if "costes" in x.lower()]
print("Dropping {} costes features".format(len(costes_cols_to_drop)))
x_df = x_df.drop(costes_cols_to_drop, axis="columns")

print(x_df.shape)
x_df.head(2)


# In[9]:


# Output
profile_file = pathlib.Path(
    f"{output_dir}/cell_health_profiles_merged_wholeplate_normalized_featureselected.tsv.gz"
)
x_df.to_csv(profile_file, index=False, sep="\t")
