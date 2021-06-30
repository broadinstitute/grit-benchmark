from pathlib import Path
import numpy as np
import pandas as pd
from scipy.special import softmax

from pycytominer.cyto_utils import infer_cp_features


def calculate_weighted_agg(
    population_df,
    columns=["Metadata_Plate", "Metadata_Well"],
    features="infer",
    transform="softmax_grit",  # or 'weighted_grit'
    weight="Metadata_grit",  # or 'Metadata_clipped_grit'
    lower_threshold=0,
    shuffle=False,
):
    """
    Aggregate population dataframe as a weighted mean with grit, with different weighting schemes.
    Currently weighting schemes: 1. weighted by grit as is (includes negative values); 2. weighted by clipped grit (with an option to specify minimum threshold, set to default 0 if not specified); 3. weighted by softmax-transformed grit.

    Arguments:
    population_df: pandas DataFrame to aggregate, must contain column of grit score
    columns: columns to groupby() for aggregation; defaults to well-level aggregations
    transform: only supports ['softmax_grit', 'weighted_grit']
    features: the feature columns to aggregate; defaults to all Cell Painting features
    weight: only supports ['Metadata_clipped_grit', 'Metadata_grit'], note Metadata pre-pended for consistency to denote metadata columns; further, currently does not support using softmax transforming clipped_grit values
    lower_threshold: lower threshold of grit to clip to
    shuffle: shuffle grit column to create synthetic negative control

    Returns:
    pandas DataFrame of mean aggregated features
    """

    if features == "infer":
        features = infer_cp_features(population_df)

    # this performs the weighted mean calculation when grit is softmax() transformed
    def weighted_sum(groupby_df, weights="Metadata_grit"):
        groupby_df["Metadata_softmax_grit"] = softmax(groupby_df[weights])
        weighted_product = groupby_df.loc[:, features].multiply(
            groupby_df["Metadata_softmax_grit"], axis=0
        )
        solution = weighted_product.sum(axis=0)
        return solution

    # this performs the weighted mean calculation when grit is either its current state or a clipped version for positive weights only (havent built in a min contribution, oops)
    def weighted_mean(groupby_df, weights="Metadata_grit"):
        solution = np.average(
            groupby_df.loc[:, features], axis=0, weights=groupby_df[weights]
        )
        return solution

    # option to shuffle the original grit column
    if shuffle:
        population_df.Metadata_grit = np.random.permutation(population_df.Metadata_grit)

    # option to clip the grit
    if weight == "Metadata_clipped_grit":
        population_df["Metadata_clipped_grit"] = population_df.Metadata_grit.clip(
            lower=lower_threshold
        )

    # branch entered for softmax-transformed grit
    if transform == "softmax_grit":
        res_df = (
            population_df.groupby(columns).apply(
                lambda x: weighted_sum(x, weights=weight)
            )
        ).reset_index()

    # branch entered for mean weighted by grit
    elif transform == "weighted_grit":
        res_df = (
            population_df.groupby(columns).apply(
                lambda x: weighted_mean(x, weights=weight)
            )
        ).reset_index()

        res_df[population_df.loc[:, features].columns] = pd.DataFrame(
            np.vstack(res_df[0].values), columns=population_df.loc[:, features].columns
        )
        res_df.drop(columns=0, inplace=True)

    return res_df
