# Self-similarity Metrics Taxonomy

The similarity matrix represents a graph with vertices and edges.

Each vertex belongs to 3 nested sets

-   Level 0 "singleton" set - i.e. just itself
-   Level 1 "replicate" set - e.g. replicates of the same perturbation
-   Level 2 "group replicate" set e.g. replicates of perturbations with same MOA

We calculate self-similarity metrics – numbers that report similarity within these sets – hierarchically:

-   **Level 1-0**: similarity of elements of a Level 0 (singleton) set to elements of its Level 1 (replicates) set, except elements of its Level 0 set. In simpler terms, this is a *replicate similarity of a vertex, i.e. the* *similarity of vertex to its replicates (except itself)*. This is a **Level 0** (singleton) set metric.

-   **Level 1**: average Level 1-0 similarity across all Level 0 (singleton) sets that are nested in the Level 1 set. In simpler terms, this is the *average replicate similarity of a set of replicate vertices*. This is a **Level 1** (replicate) set metric.

-   **Level 2-1**: similarity of each element of a Level 1 set (replicates) to elements of its Level 2 set (group replicates), except elements of its Level 1 set (replicates). *In simpler terms, this is a group replicate similarity of a replicate set, i.e. similarity of elements of a replicate set to its group replicates (except to other elements of its replicate set)*. This is a **Level 1** (replicate) set metric.

-   **Level 2**: average Level 2-1 similarity across all Level 1 (replicate) sets that are nested in the Level 2 set. In simpler terms, this is the *average group replicate similarity of a set of replicate sets.* This is a a **Level 2** (group replicate) set metric.

Consider a compound perturbation experiment done in replicates in a multi-well plate. Each compound belongs to one (or more) MOAs.

-   Each **replicate well** has a Level 1-0 metric, which is the similarity of that well to its replicates.

-   Each **compound** has a Level 1 metric, which is the average Level 1-0 metric across all its replicate wells.

-   Each **compound** has a Level 2-1 metric, which is the average similarity of each of it's replicate wells to replicate wells of other compounds with the same MOA.

-   Each **MOA** has a Level 2 metric, which is the average Level 2-1 metric across all its compounds.

## Level 1-0

### Similarities

| Index | Metric       | Description                                           |
|:------|:-------------|:------------------------------------------------------|
| 1     | `sim_mean_i` | mean similarity of a vertex to its replicate vertices |

Also

| Index | Metric         | Description                                              |
|:------|:---------------|:---------------------------------------------------------|
| 2     | `sim_median_i` | median similarity of a vertex to its replicate vertices. |

### Scaling parameters

#### With respect to non-replicates

| Index | Metric                    | Description                                                  |
|:------|:--------------------------|:-------------------------------------------------------------|
| 1     | `sim_mean_stat_non_rep_i` | mean of similarity of a vertex to its non-replicate vertices |
| 2     | `sim_sd_stat_non_rep_i`   | s.d. of similarity of a vertex to its non-replicate vertices |

#### With respect to reference

| Index | Metric                | Description                                          |
|:------|:----------------------|:-----------------------------------------------------|
| 1     | `sim_mean_stat_ref_i` | mean of similarity of a vertex to reference vertices |
| 2     | `sim_sd_stat_ref_i`   | s.d. of similarity of a vertex to reference vertices |

### Scaled similarities

#### With respect to non-replicates

| Index | Metric                      | Description                                                                    |
|:------|:----------------------------|:-------------------------------------------------------------------------------|
| 1     | `sim_scaled_mean_non_rep_i` | scale `sim_mean_i` using `sim_mean_stat_non_rep_i` and `sim_sd_stat_non_rep_i` |

Also

| Index | Metric                        | Description                                                                      |
|:------|:------------------------------|:---------------------------------------------------------------------------------|
| 2     | `sim_scaled_median_non_rep_i` | scale `sim_median_i` using `sim_mean_stat_non_rep_i` and `sim_sd_stat_non_rep_i` |

#### With respect to reference

| Index | Metric                  | Description                                                            |
|:------|:------------------------|:-----------------------------------------------------------------------|
| 1     | `sim_scaled_mean_ref_i` | scale `sim_mean_i` using `sim_mean_stat_ref_i` and `sim_sd_stat_ref_i` |

Also

| Index | Metric                    | Description                                                              |
|:------|:--------------------------|:-------------------------------------------------------------------------|
| 2     | `sim_scaled_median_ref_i` | scale `sim_median_i` using `sim_mean_stat_ref_i` and `sim_sd_stat_ref_i` |

## Level 1 (summaries of Level 1-0)

### Summaries of Level 1-0 similarities

| Index | Metric              | Description                                                        |
|:------|:--------------------|:-------------------------------------------------------------------|
| 1     | `sim_mean_i_mean_i` | mean `sim_mean_i` across all replicate vertices in a replicate set |

Also

| Index | Metric                  | Description                                                            |
|:------|:------------------------|:-----------------------------------------------------------------------|
| 2     | `sim_mean_i_median_i`   | median `sim_mean_i` across all replicate vertices in a replicate set   |
| 3     | `sim_median_i_mean_i`   | mean `sim_median_i` across all replicate vertices in a replicate set   |
| 4     | `sim_median_i_median_i` | median `sim_median_i` across all replicate vertices in a replicate set |

### Summaries of Level 1-0 scaling parameters

Note: These are summaries of scaling parameters; they are not used for scaling, themselves.

#### With respect to non-replicates

| Index | Metric                           | Description                                                                     |
|:------|:---------------------------------|:--------------------------------------------------------------------------------|
| 1     | `sim_mean_stat_non_rep_i_mean_i` | mean `sim_mean_stat_non_rep_i` across all replicate vertices in a replicate set |
| 2     | `sim_sd_stat_non_rep_i_mean_i`   | mean `sim_sd_stat_non_rep_i` across all replicate vertices in a replicate set   |

Also

| Index | Metric                             | Description                                                                       |
|:------|:-----------------------------------|:----------------------------------------------------------------------------------|
| 3     | `sim_mean_stat_non_rep_i_median_i` | median `sim_mean_stat_non_rep_i` across all replicate vertices in a replicate set |
| 4     | `sim_sd_stat_non_rep_i_median_i`   | median `sim_sd_stat_non_rep_i` across all replicate vertices in a replicate set   |

#### With respect to reference

| Index | Metric                       | Description                                                                 |
|:------|:-----------------------------|:----------------------------------------------------------------------------|
| 1     | `sim_mean_stat_ref_i_mean_i` | mean `sim_mean_stat_ref_i` across all replicate vertices in a replicate set |
| 2     | `sim_sd_stat_ref_i_mean_i`   | mean `sim_sd_stat_ref_i` across all replicate vertices in a replicate set   |

Also

| Index | Metric                         | Description                                                                   |
|:------|:-------------------------------|:------------------------------------------------------------------------------|
| 3     | `sim_mean_stat_ref_i_median_i` | median `sim_mean_stat_ref_i` across all replicate vertices in a replicate set |
| 4     | `sim_sd_stat_ref_i_median_i`   | median `sim_sd_stat_ref_i` across all replicate vertices in a replicate set   |

### Summaries of Level 1-0 scaled similarities

#### With respect to non-replicates

| Index | Metric                             | Description                                                                       |
|:------|:-----------------------------------|:----------------------------------------------------------------------------------|
| 1     | `sim_scaled_mean_non_rep_i_mean_i` | mean `sim_scaled_mean_non_rep_i` across all replicate vertices in a replicate set |

Also

| Index | Metric                                 | Description                                                                         |
|:------|:---------------------------------------|:------------------------------------------------------------------------------------|
| 2     | `sim_scaled_median_non_rep_i_median_i` | median `sim_scaled_mean_non_rep_i` across all replicate vertices in a replicate set |

#### With respect to reference

| Index | Metric                         | Description                                                                   |
|:------|:-------------------------------|:------------------------------------------------------------------------------|
| 1     | `sim_scaled_mean_ref_i_mean_i` | mean `sim_scaled_mean_ref_i` across all replicate vertices in a replicate set |

Also

| Index | Metric                             | Description                                                                     |
|:------|:-----------------------------------|:--------------------------------------------------------------------------------|
| 2     | `sim_scaled_median_ref_i_median_i` | median `sim_scaled_mean_ref_i` across all replicate vertices in a replicate set |

## Level 2-1

### Similarities

| Index | Metric       | Description                                                 |
|:------|:-------------|:------------------------------------------------------------|
| 1     | `sim_mean_g` | mean similarity of a vertex to its group replicate vertices |

Also

| Index | Metric         | Description                                                   |
|:------|:---------------|:--------------------------------------------------------------|
| 2     | `sim_median_g` | median similarity of a vertex to its group replicate vertices |

### Scaling parameters

#### With respect to non-replicates

| Index | Metric                    | Description |
|:------|:--------------------------|:------------|
| 1     | `sim_mean_stat_non_rep_g` |             |
| 2     | `sim_sd_stat_non_rep_g`   |             |

#### With respect to reference

| Index | Metric                | Description |
|:------|:----------------------|:------------|
| 3     | `sim_mean_stat_ref_g` |             |
| 4     | `sim_sd_stat_ref_g`   |             |

### Scaled similarities

#### With respect to non-replicates

| Index | Metric                      | Description |
|:------|:----------------------------|:------------|
| 1     | `sim_scaled_mean_non_rep_g` |             |

Also

| Index | Metric                        | Description |
|:------|:------------------------------|:------------|
| 2     | `sim_scaled_median_non_rep_g` |             |

#### With respect to reference

| Index | Metric                  | Description |
|:------|:------------------------|:------------|
| 1     | `sim_scaled_mean_ref_g` |             |

Also

| Index | Metric                    | Description |
|:------|:--------------------------|:------------|
| 2     | `sim_scaled_median_ref_g` |             |

## Level 2 (summaries of Level 2-1)

### Summaries of Level 2-1 similarities

### Summaries of Level 2-1 scaling parameters

### Summaries of Level 2-1 scaled similarities
