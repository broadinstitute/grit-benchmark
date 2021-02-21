#' Title
#'
#' @param population
#' @param annotation_prefix
#'
#' @return
#' @export
#'
#' @examples
get_annotation <-
  function(population,
           annotation_prefix = "Metadata_") {
    profiles %>%
      select(matches(annotation_prefix)) %>%
      mutate(id = seq(nrow(profiles))) %>%
      select(id, everything())
  }


#' Title
#'
#' @param population
#' @param annotation_prefix
#'
#' @return
#' @export
#'
#' @examples
drop_annotation <-
  function(population,
           annotation_prefix = "Metadata_") {
    population %>%
      select(-matches(annotation_prefix))
  }

# Alias
get_matrix <- drop_annotation

#' Title
#'
#' @param population
#' @param annotation_prefix
#' @param method
#'
#' @return
#' @export
#'
#' @examples
sim_calculate <-
  function(population,
           annotation_prefix = "Metadata_",
           method = "pearson") {
    # get data matrix
    data_matrix <-
      population %>%
      select(-matches(annotation_prefix))
    
    # get metadata
    metadata <-
      population %>%
      select(matches(annotation_prefix)) %>%
      rowid_to_column(var = "id")
    
    # measure similarities between treatments
    stopifnot(method %in% c("pearson", "kendall", "spearman"))
    
    sim_df <- cor(t(data_matrix), method = method)
    
    colnames(sim_df) <- seq(1, ncol(sim_df))
    
    sim_df %<>%
      as_tibble() %>%
      rowid_to_column(var = "id1") %>%
      gather(id2, sim, -id1) %>%
      mutate(id2 = as.integer(id2)) %>%
      filter(id1 != id2)
    
    # do this instead of adding a column because adding a fourth, character
    # column (<16 chars) increases the size of the data.frame by ~50%
    attr(sim_df, "method") <- method
    
    sim_df
  }

#' Title
#'
#' @param sim_df
#' @param metadata
#' @param grouping_cols
#' @param drop_group_tag_col
#' @param annotation_cols
#'
#' @return
#' @export
#'
#' @examples
sim_filter_groups <-
  function(sim_df,
           metadata,
           grouping_cols,
           annotation_cols = NULL,
           drop_group_tag_col = TRUE) {
    metadata_i <-
      metadata %>%
      select(id, any_of(grouping_cols))
    
    ids <-
      full_join(metadata_i,
                metadata_i,
                by = grouping_cols,
                suffix = c("1", "2")) %>%
      unite("group", any_of(grouping_cols), sep = ":") %>%
      select(id1, id2, group)
    
    sim_df %<>%
      inner_join(ids, by = c("id1", "id2"))
    
    if (!is.null(annotation_cols)) {
      sim_df %<>%
        sim_annotate(metadata,
                     annotation_cols,
                     index = "left")
      
    }
    
    if (drop_group_tag_col) {
      sim_df %<>% select(-group)
    }
    
    sim_df
    
  }

#' Title
#'
#' @param sim_df
#' @param metadata
#' @param group_left
#' @param group_right
#'
#' @return
#' @export
#'
#' @examples
sim_filter_reference <-
  function(sim_df,
           metadata,
           group_right = NULL,
           group_left = NULL) {
    if (!is.null(group_right)) {
      reference_ids <-
        metadata %>%
        inner_join(group_right, by = colnames(group_right)) %>%
        select(id2 = id)
      
      sim_df %<>%
        inner_join(reference_ids, by = c("id2"))
      
    }
    
    if (!is.null(group_left)) {
      reference_ids <-
        metadata %>%
        inner_join(group_left, by = colnames(group_left)) %>%
        select(id1 = id)
      
      sim_df %<>%
        inner_join(reference_ids, by = c("id1"))
      
    }
    
    sim_df
  }

#' Title
#'
#' @param sim_df
#' @param metadata
#' @param annotation_cols
#' @param index
#'
#' @return
#' @export
#'
#' @examples
sim_annotate <-
  function(sim_df,
           metadata,
           annotation_cols,
           index = "both") {
    metadata_i <-
      metadata %>%
      select(id, any_of(annotation_cols))
    
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
#' @param grouping_cols
#' @param reference
#' @param filter_out_reference
#' @param annotate_grouping_cols
#'
#' @return
#' @export
#'
#' @examples
sim_filter_reference_within_group <- function(sim_df,
                                              metadata,
                                              grouping_cols,
                                              reference,
                                              annotation_cols = NULL,
                                              filter_out_reference = TRUE)
{
  sim_df %<>%
    sim_filter_groups(metadata,
                      grouping_cols) %>%
    sim_filter_reference(metadata,
                         reference)
  
  if (filter_out_reference) {
    sim_df %<>%
      sim_annotate(metadata,
                   colnames(reference),
                   index = "left") %>%
      anti_join(reference, by = colnames(reference)) %>%
      select(-all_of(colnames(reference)))
  }
  
  if (!is.null(annotation_cols)) {
    sim_df %<>%
      sim_annotate(metadata,
                   annotation_cols,
                   index = "left")
    
  }
  
  sim_df
}
