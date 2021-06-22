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

def merge_metadata(cell_line, level3_profile):
    # load metadata file from Cell Health  data
    commit = "07e4b40c39dd27084be36fbef4d64c5654b2960f"
    base_url = f"https://github.com/broadinstitute/cell-health/raw/{commit}"
    url = f"{base_url}/1.generate-profiles/data/metadata/platemap/DEPENDENCIES1_{cell_line}.csv"
    platemap = pd.read_csv(url, sep=",")
    platemap.columns = ['Metadata_' + str(x) for x in platemap.columns]

    # merge with the aggregated files
    meta_df = pd.merge(level3_profile, platemap, left_on='Metadata_Well', right_on='Metadata_well_position')
    # # reorder columns for metadata to be in front
    meta_df = meta_df[sorted(meta_df, key = lambda x: x not in meta_df.filter(like="Metadata").columns)]

    return meta_df

### load level 2 data ###
plate_dict = {
    "ES2": ["SQ00014613","SQ00014614", "SQ00014615"], 
    "A549": ["SQ00014610", "SQ00014611", "SQ00014612"],
    "HCC44": ["SQ00014616", "SQ00014617", "SQ00014618"]
}

### take the same columns as original Cell Health paper did ###
commit = "2916770cc9cc9e75b693348b683aa398987fb9f9"
base_url = f"https://github.com/broadinstitute/grit-benchmark/raw/{commit}"
url = f"{base_url}/1.calculate-metrics/cell-health/data/cell_health_merged_feature_select.csv.gz"

df = pd.read_csv(url, sep=",")
print(df.shape)
# df.head(2)

cols_to_keep = infer_cp_features(df)

### do the aggregation ###
results_folder = '../data/aggregated-profiles/'
    
for cell_line in ['A549',]: #'HCC44', 'ES2', 
    start_process = datetime.now()
    print("starting to read files for {} at {}".format(cell_line, start_process))
          
        
    ####### read in single-cell grit data #######
    start_merge = datetime.now()
    grit_folder = '../../../../1.calculate-metrics/cell-health/results/'
    grit_files = glob.glob(grit_folder+'*single_cell_grit*.tsv.gz')

    scgrit_df = []
    for file in grit_files:
        plate_name=file.split('/')[-1].split('_')[-2]
        if plate_name in plate_dict[cell_line]:
            print(f"adding scrgrit of {plate_name} to list of {cell_line}")
            scgrit_plate = pd.read_csv(file, sep='\t').assign(plate=plate_name, cell_line = cell_line)
            print(scgrit_plate.shape)
            scgrit_df.append(scgrit_plate)
    scgrit_df = pd.concat(scgrit_df)
    scgrit_df['cell_identity'] = scgrit_df.perturbation.str.split("_", expand=True)[1].astype(int)
    scgrit_df.columns = ['Metadata_'  + str(col) for col in scgrit_df.columns]
    print(f"total shape of of scgrit_df for {cell_line} is: {scgrit_df.shape}")

    ####### read in single-cell cell painting profiles #######
    profile_folder = '../../../../0.download-data/data/cell_health/normalized/' 
    profile_files = glob.glob(profile_folder+'*normalized.csv.gz')

    scprofiles_df = []
    for file in profile_files:
        plate_name=file.split('/')[-1].split('_')[0]
        if plate_name in plate_dict[cell_line]:
            print(f"adding scprofiles of {plate_name} to list of {cell_line}")
            scprofile_plate = (pd.read_csv(file, sep=',', low_memory=False)
                               .reset_index()
                               .rename({'index':'Metadata_cell_identity'}, axis='columns')
                              ).assign(cell_line = cell_line)
            plate_cols = infer_cp_features(scprofile_plate)
            drop_cols = [x for x in plate_cols if x not in cols_to_keep]
            scprofile_plate.drop(columns = drop_cols, inplace=True)
            scprofiles_df.append(scprofile_plate)
    scprofiles_df = pd.concat(scprofiles_df, sort=False)
    print(f'total shape of scprofiles_df for {cell_line} is: {scprofiles_df.shape}')
    

    ####### merge scgrit scores + cell painting profiles #######
    scprofiles_df = (pd.merge(scprofiles_df, scgrit_df, 
         left_on=['Metadata_cell_identity', 'Metadata_Plate', 'Metadata_pert_name'], 
                 right_on=['Metadata_cell_identity', 'Metadata_plate', 'Metadata_group'])
        )
    del scgrit_df
    print(f"total shape of sc_df for {cell_line} is: {scprofiles_df.shape}")
    # remove columns with any NA entries
    na_cols_to_drop = get_na_columns(scprofiles_df, cutoff=0)
    print(f"Dropping {len(na_cols_to_drop)} columns because of missing data")
    scprofiles_df = scprofiles_df.drop(na_cols_to_drop, axis="columns")
    print(f"FINAL shape of merged data {scprofiles_df.shape}")

    print(f"TOTAL TIME constructing merged df for cell_line {cell_line} : {str(datetime.now()-start_merge)}")

    
    ###### standard median aggregation ######
    start_agg = datetime.now()
    agg_df = aggregate(
        population_df = scprofiles_df,
        strata = ["Metadata_Plate", "Metadata_Well"],
        features = "infer",
        operation ="median"
    ).assign(Metadata_agg_method = 'median', cell_line = cell_line)
    agg_meta_df = merge_metadata(cell_line, agg_df)
    # writing data
    agg_meta_df.to_csv(Path(results_folder + cell_line + "_median.tsv"), index=False, sep='\t')
    
    ###### grit-informed aggregation methods ######
    ### raw grit as weights ###
    agg_df = (calculate_weighted_agg(
        population_df = scprofiles_df,
        columns = ['Metadata_Plate', 'Metadata_Well'],
        features = 'infer',
        transform = 'weighted_grit', weight = 'Metadata_grit')
                    ).assign(Metadata_agg_method = 'weighted', cell_line = cell_line)
    agg_meta_df = merge_metadata(cell_line, agg_df)
    # writing data
    agg_meta_df.to_csv(Path(results_folder + cell_line + "_weighted.tsv"), index=False, sep='\t')
    

    print(f"TOTAL TIME performing aggregation for cell_line {cell_line} : {str(datetime.now()-start_agg)}")

    
