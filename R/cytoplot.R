#' Plot platemap
#' 
#' \code{plot_plate} plots platemaps
#'
#' @param plate 
#' @param variable 
#' @param well_position 
#'
#' @return
#' @export
#'
#' @examples
plot_plate <- function(plate,
                       variable,
                       well_position = "well_position") {
  
  variable <- rlang::sym(variable)
  
  well_position <- rlang::sym(well_position)
  
  plate %<>%
    rowwise() %>%
    mutate(row_id =
             str_sub(!!well_position, 1, 1)[[1]]) %>%
    mutate(row_id_int =
             str_locate(paste0(LETTERS, collapse = ""), row_id)[[1]]) %>%
    mutate(col_id =
             as.integer(str_sub(!!well_position, 2, 3)))
  
  p <- 
    ggplot(plate, aes(as.factor(col_id), fct_rev(as.factor(row_id)), fill = !!variable)) +
    geom_tile(color = "white") +
    coord_equal() +
    theme_minimal() +
    scale_x_discrete(name = "", position = "top") +
    scale_y_discrete(name = "")
  
  p
}
