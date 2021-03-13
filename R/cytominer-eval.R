sim_cols <- c("id1", "id2", "sim")

#' Get row annotations.
#' 
#' \code{get_annotation} gets row annotations.
#'
#' @param population tbl with annotations (a.k.a. metadata) and observation variables.
#' @param annotation_prefix optional character string specifying  prefix for annotation columns.
#'
#' @return row annotations of the same class as \code{population}
#' @export
#' 
#' @importFrom magrittr %>%
#'
#'
#' @examples
#' suppressMessages(suppressWarnings(library(magrittr)))
#' population <- tibble::tibble(
#'   Metadata_group = c(
#'     "control", "control", "control", "control",
#'     "experiment", "experiment", "experiment", "experiment"
#'   ),
#'   Metadata_batch = c("a", "a", "b", "b", "a", "a", "b", "b"),
#'   AreaShape_Area = c(10, 12, 15, 16, 8, 8, 7, 7)
#' )
#' cytoeval::get_annotation(population, annotation_prefix = "Metadata_")
#' @export
get_annotation <-
  function(population,
           annotation_prefix = "Metadata_") {
    population %>%
      dplyr::select(dplyr::matches(annotation_prefix)) %>%
      dplyr::mutate(id = seq(nrow(population))) %>%
      dplyr::select(id, dplyr::everything())
  }


#' Drop row annotations.
#' 
#' \code{drop_annotation} drops row annotations.
#'
#' @param population tbl with annotations (a.k.a. metadata) and observation variables.
#' @param annotation_prefix optional character string specifying prefix for annotation columns.
#'
#' @return data with all columns except row annotations of the same class as \code{population}
#' @export
#' 
#' @importFrom magrittr %>%
#'
#' @examples
#' suppressMessages(suppressWarnings(library(magrittr)))
#' population <- tibble::tibble(
#'   Metadata_group = c(
#'     "control", "control", "control", "control",
#'     "experiment", "experiment", "experiment", "experiment"
#'   ),
#'   Metadata_batch = c("a", "a", "b", "b", "a", "a", "b", "b"),
#'   AreaShape_Area = c(10, 12, 15, 16, 8, 8, 7, 7)
#' )
#' cytoeval::drop_annotation(population, annotation_prefix = "Metadata_")
#' @export 
drop_annotation <-
  function(population,
           annotation_prefix = "Metadata_") {
    population %>%
      dplyr::select(-dplyr::matches(annotation_prefix))
  }

#' Calculate melted similarity matrix.
#' 
#' \code{sim_calculate} calculates a melted similarity matrix.
#'
#' @param population tbl with annotations (a.k.a. metadata) and observation variables.
#' @param annotation_prefix optional character string specifying prefix for annotation columns.
#' @param method optional character string specifying method for \code{stats::cor} to calculate similarity.  This must be one of the strings \code{"pearson"} (default), \code{"kendall"}, \code{"spearman"}.
#'
#' @return data.frame of melted similarity matrix.
#' 
#' @importFrom magrittr %>%
#'
#' @examples
#' suppressMessages(suppressWarnings(library(magrittr)))
#' population <- tibble::tibble(
#'   Metadata_group = sample(c('a', 'b'), 4, replace = TRUE),
#'   x = rnorm(4),
#'   y = x +rnorm(4) / 100,
#'   z = y + rnorm(4) / 1000
#' )
#' cytoeval::sim_calculate(population, method = "pearson")
#' @export 
sim_calculate <-
  function(population,
           annotation_prefix = "Metadata_",
           method = "pearson") {
    # get data matrix
    data_matrix <-
      population %>%
      dplyr::select(-dplyr::matches(annotation_prefix))
    
    # get metadata
    metadata <-
      population %>%
      dplyr::select(dplyr::matches(annotation_prefix)) %>%
      tibble::rowid_to_column(var = "id")
    
    # measure similarities between treatments
    stopifnot(method %in% c("pearson", "kendall", "spearman"))
    
    sim_df <- stats::cor(t(data_matrix), method = method)
    
    colnames(sim_df) <- seq(1, ncol(sim_df))
    
    sim_df %<>%
      tibble::as_tibble() %>%
      tibble::rowid_to_column(var = "id1") %>%
      tidyr::pivot_longer(-id1, names_to = "id2", values_to = "sim") %>%
      dplyr::mutate(id2 = as.integer(id2)) %>%
      dplyr::filter(id1 != id2)
    
    # do this instead of adding a column because adding a fourth, character
    # column (<16 chars) increases the size of the data.frame by ~50%
    attr(sim_df, "method") <- method
    
    sim_df
  }

#' Annotate melted similarity matrix.
#' 
#' \code{sim_annotate} annotate a melted similarity matrix.
#'
#' @param sim_df tbl with melted similarity matrix. 
#' @param metadata tbl with metadata annotations with which to annotate the similarity matrix.
#' @param annotation_cols character vector specifying annotation columns.
#' @param index optional character string specifying whether to annotate left index, right index, or both.  This must be one of the strings \code{"both"} (default), \code{"left"}, \code{"right"}.
#'
#' @return annotated melted similarity matrix of the same class as \code{sim_df}.
#' 
#' @importFrom magrittr %>%
#'
#' @examples
#' suppressMessages(suppressWarnings(library(magrittr)))
#' population <- tibble::tibble(
#'   Metadata_group = sample(c('a', 'b'), 4, replace = TRUE),
#'   Metadata_type = sample(c('x', 'y'), 4, replace = TRUE),
#'   x = rnorm(4),
#'   y = x +rnorm(4) / 100,
#'   z = y + rnorm(4) / 1000
#' )
#' metadata <- cytoeval::get_annotation(population)
#' annotation_cols <- c("Metadata_group")
#' sim_df <- cytoeval::sim_calculate(population, method = "pearson")
#' cytoeval::sim_annotate(sim_df, metadata, annotation_cols)
#' @export
sim_annotate <-
  function(sim_df,
           metadata,
           annotation_cols,
           index = "both") {
    metadata_i <-
      metadata %>%
      select(id, any_of(annotation_cols))
    
    sim_df %<>% select(all_of(sim_cols))
    
    if (index == "left") {
      sim_df %<>%
        inner_join(metadata_i,
                   by = c("id1" = "id"),
                   suffix = c("1", "2"))
    }
    
    if (index == "right") {
      sim_df %<>%
        inner_join(metadata_i,
                   by = c("id2" = "id"),
                   suffix = c("1", "2"))
    }
    
    if (index == "both") {
      sim_df %<>%
        inner_join(metadata_i,
                   by = c("id1" = "id")) %>%
        inner_join(metadata_i,
                   by = c("id2" = "id"),
                   suffix = c("1", "2"))
    }
    
    sim_df
    
  }

#' Title
#'
#' @param sim_df
#' @param metadata
#' @param all_same_cols
#' @param annotation_cols
#' @param include_group_tag
#' @param drop_lower
#'
#' @return
#' @export
#'
#' @examples
sim_all_same <-
  function(sim_df,
           metadata,
           all_same_cols,
           annotation_cols = NULL,
           include_group_tag = FALSE,
           drop_lower = FALSE) {
    metadata_i <-
      metadata %>%
      select(id, all_of(all_same_cols)) %>%
      unite("all_same_col", all_of(all_same_cols), sep = ":")
    
    ids <-
      inner_join(metadata_i,
                 metadata_i,
                 by = "all_same_col",
                 suffix = c("1", "2"))
    
    if (include_group_tag) {
      ids %<>% select(id1, id2, group = all_same_col)
      
    } else {
      ids %<>% select(id1, id2)
    }
    
    if (drop_lower) {
      sim_df %<>% filter(id1 > id2)
    }
    
    sim_df %<>%
      inner_join(ids, by = c("id1", "id2"))
    
    if (!is.null(annotation_cols)) {
      sim_df %<>%
        sim_annotate(metadata,
                     annotation_cols,
                     index = "left")
    }
    
    sim_df
    
  }


#' Title
#'
#' @param sim_df
#' @param metadata
#' @param filter_keep
#' @param filter_side
#'
#' @return
#' @export
#'
#' @examples
sim_keep <-
  function(sim_df,
           metadata,
           filter_keep,
           filter_side) {
    stopifnot(filter_side %in% c("left", "right"))
    
    filter_ids <-
      metadata %>%
      inner_join(filter_keep, by = colnames(filter_keep)) %>%
      select(id)
    
    join_str <- c("id")
    
    # join_str will be either c("id1" = "id") or c("id2" = "id")
    names(join_str) <-
      paste0("id", ifelse(filter_side == "left", 1, 2))
    
    sim_df %<>%
      inner_join(filter_ids, by = join_str)
    
    sim_df
    
  }


#' Title
#'
#' @param sim_df
#' @param metadata
#' @param filter_drop
#' @param filter_side
#'
#' @return
#' @export
#'
#' @examples
sim_drop <-
  function(sim_df,
           metadata,
           filter_drop,
           filter_side) {
    stopifnot(filter_side %in% c("left", "right"))
    
    if (is.null(filter_drop)) {
      return(sim_df)
    }
    
    sim_df %<>%
      select(all_of(sim_cols)) %>%
      sim_annotate(metadata,
                   colnames(filter_drop),
                   index = filter_side) %>%
      anti_join(filter_drop, by = colnames(filter_drop)) %>%
      select(all_of(sim_cols))
    
    sim_df
  }



#' Title
#'
#' @param sim_df
#' @param metadata
#' @param all_same_cols
#' @param filter_keep_right
#' @param annotation_cols
#' @param drop_reference
#'
#' @return
#' @export
#'
#' @examples
sim_all_same_keep_some <-
  function(sim_df,
           metadata,
           all_same_cols,
           filter_keep_right,
           annotation_cols = NULL,
           drop_reference = TRUE)
  {
    sim_df %<>%
      sim_all_same(metadata,
                   all_same_cols) %>%
      sim_keep(metadata,
               filter_keep_right,
               "right")
    
    if (drop_reference) {
      filter_drop_left <- filter_keep_right
      
      sim_df %<>%
        sim_drop(metadata,
                 filter_drop_left,
                 "left")
    }
    
    if (!is.null(annotation_cols)) {
      sim_df %<>%
        select(all_of(sim_cols)) %>%
        sim_annotate(metadata,
                     annotation_cols,
                     index = "left")
      
    }
    
    sim_df
  }

#' Title
#'
#' @param sim_df
#' @param metadata
#' @param any_different_cols
#' @param all_same_cols
#' @param all_different_cols
#' @param filter_keep_left
#' @param filter_keep_right
#' @param annotation_cols
#'
#' @return
#' @export
#'
#' @examples
sim_some_different_drop_some <-
  function(sim_df,
           metadata,
           any_different_cols,
           all_same_cols = NULL,
           all_different_cols = NULL,
           filter_drop_left = NULL,
           filter_drop_right = NULL,
           annotation_cols = NULL) {
    stopifnot(!any(all_same_cols %in% all_different_cols))
    
    metadata_i <- metadata
    
    if (is.null(all_same_cols)) {
      # create a dummy column on which to join
      metadata_i %<>% mutate(all_same_col = 0)
      all_same_cols <- "all_same_col"
    } else {
      # create a unified column on which to join
      metadata_i %<>%
        unite("all_same_col",
              all_of(all_same_cols),
              sep = ":",
              remove = FALSE)
    }
    
    # ignore any_different_cols if superseded by all_different_cols
    if (any(all_different_cols %in% any_different_cols)) {
      any_different_cols <- NULL
    }
    
    # remove from any_different_cols its intersection with all_same_cols
    any_different_cols <- setdiff(any_different_cols, all_same_cols)
    
    # create a unified column for any_different_cols and include that new column
    # in all_different_cols
    if (!is.null(any_different_cols)) {
      metadata_i %<>%
        unite(
          "any_different_col",
          all_of(any_different_cols),
          sep = ":",
          remove = FALSE
        )
      
      all_different_cols <-
        c(all_different_cols, "any_different_col")
      
    }
    
    # create left and right metadata
    f_metadata_filter <-
      function(filter_drop) {
        if (is.null(filter_drop)) {
          metadata_i %>%
            select(id, all_same_col)
          
        } else {
          metadata_i %>%
            anti_join(filter_drop, by = colnames(filter_drop)) %>%
            select(id, all_same_col)
          
        }
      }
    
    metadata_left  <- f_metadata_filter(filter_drop_left)
    metadata_right <- f_metadata_filter(filter_drop_right)
    
    # list of rows that should be the same (weak constraint)
    ids_all_same <-
      inner_join(metadata_left,
                 metadata_right,
                 by = "all_same_col",
                 suffix = c("1", "2"))
    
    # list of rows that should be the different (strong constraint)
    ids_all_different <-
      map_df(all_different_cols,
             function(all_different_col) {
               inner_join(
                 metadata_i %>% select(id, all_of(all_different_col)),
                 metadata_i %>% select(id, all_of(all_different_col)),
                 by = all_different_col,
                 suffix = c("1", "2")
               ) %>%
                 select(id1, id2)
               
             }) %>%
      distinct()
    
    # impose strong constraint on weak constraint
    ids <-
      anti_join(ids_all_same,
                ids_all_different,
                by = c("id1", "id2"))
    
    ids %<>% select(id1, id2)
    
    # filter similarity matrix
    sim_df %<>%
      inner_join(ids, by = c("id1", "id2"))
    
    # add annotations
    if (!is.null(annotation_cols)) {
      sim_df %<>%
        sim_annotate(metadata,
                     annotation_cols,
                     index = "left")
    }
    
    sim_df
    
  }



#' Title
#'
#' @param sim_df
#' @param metadata
#' @param reference
#' @param all_same_cols_rep
#' @param all_same_cols_rep_ref
#' @param all_same_cols_ref
#' @param any_different_cols_non_rep
#' @param all_same_cols_non_rep
#' @param all_different_cols_non_rep
#' @param annotation_cols
#' @param drop_group
#'
#' @return
#' @export
#'
#' @examples
sim_munge <-
  function(sim_df,
           metadata,
           all_same_cols_rep,
           annotation_cols,
           all_same_cols_ref = NULL,
           all_same_cols_rep_ref = NULL,
           any_different_cols_non_rep = NULL,
           all_same_cols_non_rep = NULL,
           all_different_cols_non_rep = NULL,
           any_different_cols_group = NULL,
           all_same_cols_group = NULL,
           reference = NULL,
           drop_reference = FALSE,
           drop_group = NULL) {
    # ---- 0. Filter out some rows ----
    
    if (!is.null(drop_group)) {
      sim_df %<>%
        sim_drop(metadata, drop_group, "left") %>%
        sim_drop(metadata, drop_group, "right")
      
    }
    
    fetch_ref <-
      !is.null(all_same_cols_ref) &&
      !is.null(reference)
    
    fetch_rep_ref <-
      !is.null(all_same_cols_ref) &&
      !is.null(reference) &&
      !is.null(all_same_cols_rep_ref)
    
    fetch_non_rep <-
      !is.null(any_different_cols_non_rep) &&
      !is.null(all_same_cols_non_rep) &&
      !is.null(all_different_cols_non_rep)
    
    fetch_rep_group <-
      !is.null(any_different_cols_group) &&
      !is.null(all_same_cols_group)
    
    # ---- 1. Similarity to reference ----
    
    # Fetch similarities between
    # a. all rows (except, optionally those containing `reference`)
    # and
    # b. all rows containing `reference`
    # Do so only for those (a, b) pairs that
    # - have *same* values in *all* columns of `all_same_cols_ref`
    
    if (fetch_ref) {
      ref <-
        sim_df %>%
        sim_all_same_keep_some(
          metadata,
          all_same_cols_ref,
          filter_keep_right = reference,
          drop_reference = drop_reference,
          annotation_cols
        )
    }
    
    # ---- 2. Similarity to replicates (no references) ----
    
    # Fetch similarities between
    # a. all rows except `reference` rows
    # and
    # b. all rows except `reference` rows (i.e. to each other)
    #
    # Do so for only those (a, b) pairs that
    # - have *same* values in *all* columns of `all_same_cols_rep
    #
    # Keep, both, (a, b) and (b, a)
    
    rep <-
      sim_df %>%
      sim_drop(metadata, reference, "left") %>%
      sim_drop(metadata, reference, "right") %>%
      sim_all_same(metadata,
                   all_same_cols_rep,
                   annotation_cols,
                   drop_lower = FALSE)
    
    # ---- 3. Similarity to replicates (only references) ----
    
    # Fetch similarities between
    # a. all rows containing `reference`
    # and
    # b. all rows containing `reference` (i.e. to each other)
    #
    # Do so for only those (a, b) pairs that
    # - have *same* values in *all* columns of `all_same_cols_rep_ref`.
    #
    # Keep, both, (a, b) and (b, a)
    
    if (fetch_rep_ref) {
      rep_ref <-
        sim_df %>%
        sim_keep(metadata,
                 filter_keep = reference,
                 "left") %>%
        sim_keep(metadata,
                 filter_keep = reference,
                 "right") %>%
        sim_all_same(
          metadata,
          all_same_cols = all_same_cols_rep_ref,
          annotation_cols = annotation_cols,
          drop_lower = FALSE
        )
    }
    
    # ---- 4. Similarity to non-replicates ----
    
    # Fetch similarities between
    # a. all rows (except, optionally, `reference` rows)
    # and
    # b. all rows except `reference` rows
    #
    # Do so for only those (a, b) pairs that
    # - have *same* values in *all* columns of `all_same_cols_non_rep`
    # - have *different* values in *all* columns `all_different_cols_non_rep`
    # - have *different* values in *at least one* column of `any_different_cols_non_rep`
    #
    # Keep, both, (a, b) and (b, a)
    
    if (fetch_non_rep) {
      if (drop_reference) {
        reference_left <- reference
      } else {
        reference_left <- NULL
      }
      
      non_rep <-
        sim_df %>%
        sim_some_different_drop_some(
          metadata,
          any_different_cols = any_different_cols_non_rep,
          all_same_cols = all_same_cols_non_rep,
          all_different_cols = all_different_cols_non_rep,
          filter_drop_left = reference_left,
          filter_drop_right = reference,
          annotation_cols = annotation_cols
        )
    }
    
    # ---- 5. Similarity to group ----
    
    # Fetch similarities between
    # a. all rows (except, optionally, `reference` rows)
    # and
    # b. all rows (except, optionally, `reference` rows)
    #
    # Do so for only those (a, b) pairs that
    # - have *same* values in *all* columns of `all_same_cols_group`
    # - have *different* values in *at least one* column of `any_different_cols_group`
    #
    # Keep, both, (a, b) and (b, a)
    
    if (fetch_rep_group) {
      if (drop_reference) {
        reference_both <- reference
      } else {
        reference_both <- NULL
      }
      
      rep_group <-
        sim_df %>%
        sim_some_different_drop_some(
          metadata,
          any_different_cols = any_different_cols_group,
          all_same_cols = all_same_cols_group,
          filter_drop_left = reference_both,
          filter_drop_right = reference_both,
          annotation_cols = annotation_cols
        )
    }
    
    # 6. Combine
    
    combined <-
      rep %>% mutate(type = "rep")
    
    if (fetch_rep_ref) {
      combined %<>%
        bind_rows(rep_ref %>% mutate(type = "rep")) # same tag as ref
    }
    
    if (fetch_non_rep) {
      combined %<>%
        bind_rows(non_rep %>% mutate(type = "non_rep"))
    }
    
    if (fetch_ref) {
      combined %<>%
        bind_rows(ref %>% mutate(type = "ref"))
    }
    
    if (fetch_rep_group) {
      combined %<>%
        bind_rows(rep_group %>% mutate(type = "rep_group"))
    }
    
    # add attributes
    
    attr(combined, "all_different_cols_non_rep") <-
      all_different_cols_non_rep
    attr(combined, "all_same_cols_group") <- all_same_cols_group
    attr(combined, "all_same_cols_non_rep") <- all_same_cols_non_rep
    attr(combined, "all_same_cols_ref") <- all_same_cols_ref
    attr(combined, "all_same_cols_rep") <- all_same_cols_rep
    attr(combined, "all_same_cols_rep_ref") <- all_same_cols_rep_ref
    attr(combined, "annotation_cols") <- annotation_cols
    attr(combined, "any_different_cols_group") <-
      any_different_cols_group
    attr(combined, "any_different_cols_non_rep") <-
      any_different_cols_non_rep
    attr(combined, "drop_group") <- drop_group
    attr(combined, "drop_reference") <- drop_reference
    attr(combined, "reference") <- reference
    
    combined
  }


#' Title
#'
#' @param grouped_sim
#' @param sim_type
#' @param calculate_grouped
#' @param annotation_prefix
#'
#' @return
#' @export
#'
#' @examples
sim_metrics <- function(munged_sim,
                        sim_type,
                        calculate_grouped = FALSE,
                        annotation_prefix = "Metadata_") {
  if (!is.null(attr(munged_sim, "all_same_cols_rep", TRUE))) {
    rep_cols <- attr(munged_sim, "all_same_cols_rep", TRUE)
  } else {
    message("Warning: Inferring columns specifying replicates from similarity dataframe...")
    rep_cols <-
      str_subset(colnames(munged_sim), pattern = annotation_prefix)
  }
  
  helper_scale_aggregate <-
    function(summary_cols,
             sim_type_replication,
             identifier = NULL) {
      sim_stats <-
        munged_sim %>%
        filter(type == sim_type) %>%
        group_by(across(all_of(summary_cols))) %>%
        summarise(across(all_of("sim"), list(mean = mean, sd = sd)),
                  .groups = "keep") %>%
        ungroup()
      
      # scale using mean and s.d. of the `sim_type` distribution
      
      join_cols <-
        intersect(colnames(munged_sim), colnames(sim_stats))
      
      sim_norm <-
        munged_sim %>%
        filter(type == sim_type_replication) %>%
        inner_join(sim_stats, by = join_cols) %>%
        mutate(sim_scaled = (sim - sim_mean) / sim_sd)
      
      # get a summary per group
      sim_norm_agg <-
        sim_norm %>%
        group_by(across(all_of(summary_cols))) %>%
        summarise(across(all_of(c(
          "sim_scaled", "sim"
        )),
        list(mean = mean, median = median)),
        .groups = "keep") %>%
        rename_with( ~ paste(., sim_type, sep = "_"),
                     starts_with("sim_scaled")) %>%
        ungroup()
      
      sim_norm_agg %<>%
        inner_join(sim_stats %>%
                     rename_with(~ paste(., "stat", sim_type, sep = "_"),
                                 starts_with("sim")),
                   by = join_cols)
      
      if (!is.null(identifier)) {
        sim_norm_agg %<>%
          rename_with(~ paste(., identifier, sep = "_"),
                      starts_with("sim"))
      }
      
      sim_norm_agg
    }
  
  # ---- Replicates ----
  sim_norm_agg <-
    helper_scale_aggregate(c("id1", rep_cols), "rep", "i")
  
  # get a summary per set
  
  sim_norm_agg_agg <-
    sim_norm_agg %>%
    ungroup() %>%
    group_by(across(all_of(c(rep_cols)))) %>%
    summarise(across(-all_of("id1"),
                     list(mean = mean, median = median)),
              .groups = "keep")
  
  # append identified ("_i" for "individual")
  
  sim_norm_agg_agg %<>%
    rename_with( ~ paste(., "i", sep = "_"),
                 starts_with("sim"))
  
  result <-
    list(per_row = sim_norm_agg,
         per_set = sim_norm_agg_agg)
  
  # ---- Group replicates  ----
  
  if (calculate_grouped) {
    sim_norm_group_agg <-
      helper_scale_aggregate(rep_cols, "rep_group", "g")
    
    result <-
      c(result,
        list(per_set_group = sim_norm_group_agg))
  }
  
  result
}
