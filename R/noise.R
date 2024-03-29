#' This function takes a network and either adds or deletes a percentage of
#' edges randomly
#'
#' @param network Network to modify.
#' @param mode Whether to add ('add') or delete ('del') edges.
#' @param perc Percentage of edges to modify per regulon.
#' @param source Network's source column.
#' @param target Network's target genes column.
#' @param seed An integer to set the RNG state for random number generation. Use
#'    NULL for random number generation.
#' @return Modified network.
#'
#' @import dplyr
#' @import purrr
#'
#' @export
net_noise <- function(network, mode='add', perc=0.1, source='source',
                          target='target', seed=42){
  # Get unique targets of network
  targets <- unique(network[[target]])

  # For each source, change % of edges
  network %>%
    group_by(!! sym(source)) %>%
    group_split() %>%
    map(function(df){
      # Name source
      name_source <- df[[source]][1]

      # Number of genes
      n_genes <- ceiling(perc * nrow(df))

      if (mode == 'add'){
        # Get genes not in the regulon
        diff <- setdiff(targets, df[[target]])

        # Sample n edges to add
        set.seed(seed)
        sampled <- sample(diff, n_genes, replace=F)

        # Append to df
        tbl <- tibble(
          mor = rep(1, length(sampled))
          )
        tbl[[source]] <- rep(name_source, length(sampled))
        tbl[[target]] <- sampled
        tbl <- bind_rows(df, tbl)
        tbl
      } else if (mode == 'del') {
        # Sample n edges to del
        set.seed(seed)
        sampled <- sample(df[[target]], n_genes, replace=F)

        # Delete rows with sampled genes
        df %>%
          filter(! (!!sym(target)) %in% sampled)
      } else {
        stop(stringr::str_glue(
          "{mode} is an invalid mode, please select 'add' or 'del'."))
      }
    }) %>%
    bind_rows() %>%
    distinct(!! sym(source), !! sym(target), .keep_all=T)
}

