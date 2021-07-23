#!/usr/bin/env python
# coding: utf-8

# ## Download and Extract Perturb-Seq Data
# 
# We dowload the perturbseq dataset from [GEO accession `GSE132080`](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE132080).
# 
# > Jost M, Santos DA, Saunders RA, Horlbeck MA et al. Titrating gene expression using libraries of systematically attenuated CRISPR guide RNAs. Nat Biotechnol 2020 Mar;38(3):355-364. PMID: 31932729
# 
# We also gunzip the data and extract to a separate folder to enable 10X data processing.

# In[1]:


import gzip
import shutil
import pathlib
import pandas as pd
from urllib.request import urlretrieve


# In[2]:


def download_file(file, base_url, output_dir):
    download_url = f"{base_url}/{file}"
    print(f"Now downloading {download_url}...")
    output_file = pathlib.Path(f"{output_dir}/{file}")

    urlretrieve(download_url, output_file)


# In[3]:


gse_id = "GSE132080"
base_url = "https://ftp.ncbi.nlm.nih.gov/geo/series/GSE132nnn/GSE132080/suppl/"
output_dir = pathlib.Path("data/perturbseq")
output_dir.mkdir(exist_ok=True)

files = {
    "barcodes": {
        "download": f"{gse_id}_10X_barcodes.tsv.gz",
        "rename": pathlib.Path(f"{gse_id}/barcodes.tsv")
    },
    "data": {
        "download": f"{gse_id}_10X_matrix.mtx.gz",
        "rename": pathlib.Path(f"{gse_id}/matrix.mtx")
    },
    "genes": {
        "download": f"{gse_id}_10X_genes.tsv.gz",
        "rename": pathlib.Path(f"{gse_id}/genes.tsv")
    },
    "other": {
        f"{gse_id}_cell_identities.csv.gz",
        f"{gse_id}_sgRNA_barcode_sequences_and_phenotypes.csv.gz",
    }
}

files


# In[4]:


for data_type in files:
    if data_type != "other":
        file = files[data_type]["download"]
        download_file(file, base_url, output_dir)
    else:
        for file in files[data_type]:
            download_file(file, base_url, output_dir)


# In[5]:


for data_type in files:
    if data_type != "other":
        file = files[data_type]["download"]
        file = pathlib.Path(f"{output_dir}/{file}")
        
        rename_and_gunzip_file = pathlib.Path(f"{output_dir}/{files[data_type]['rename']}")
        rename_and_gunzip_file.parent.mkdir(exist_ok=True)
        
        print(f"Now extracting {file} to {rename_and_gunzip_file}...")
        with gzip.open(file, "rb") as f_in:
            with open(rename_and_gunzip_file, "wb") as f_out:
                shutil.copyfileobj(f_in, f_out)

