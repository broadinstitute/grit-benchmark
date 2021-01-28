#!/usr/bin/env python
# coding: utf-8

# ## Download CERES metrics
#
# We download CERES scores from DepMap 20Q3 public release. The data are available at https://depmap.org/portal/download/.
#
# > DepMap, Broad (2020): DepMap 20Q3 Public. figshare. Dataset doi:10.6084/m9.figshare.12931238.v1.
#
# The CERES score was developed to estimate gene dependencies from CRISPR screens accounting for copy number impact.
#
# > Robin M. Meyers, Jordan G. Bryan, James M. McFarland, Barbara A. Weir, ... David E. Root, William C. Hahn, Aviad Tsherniak. Computational correction of copy number effect improves specificity of CRISPR-Cas9 essentiality screens in cancer cells. Nature Genetics 2017 October 49:1779–1784. doi:10.1038/ng.3984

# In[1]:


import pathlib
import pandas as pd
from urllib.request import urlretrieve


# In[2]:


figshare_base_url = "https://ndownloader.figshare.com/files/"

file_info = {
    "ceres": {"file_id": "24613292", "output_file": "ceres.csv"},
    "sample_id": {"file_id": "24613394", "output_file": "depmap_sample_info.csv"},
}

output_dir = pathlib.Path("data")
output_dir.mkdir(exist_ok=True)


# In[3]:


for data in file_info:
    file_id = file_info[data]["file_id"]
    output_file = file_info[data]["output_file"]

    download_url = f"{figshare_base_url}/{file_id}"
    output_file = pathlib.Path(f"{output_dir}/{output_file}")

    urlretrieve(download_url, output_file)


# ## Preview downloaded data

# In[4]:


# Load CERES scores
output_file = file_info["ceres"]["output_file"]
df = pd.read_csv(f"{output_dir}/{output_file}")

print(df.shape)
df.head(3)


# In[5]:


# Load sample ID
output_file = file_info["sample_id"]["output_file"]
sample_df = pd.read_csv(f"{output_dir}/{output_file}")

print(sample_df.shape)
sample_df.head(3)
