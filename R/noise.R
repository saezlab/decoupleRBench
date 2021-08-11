#' This function takes a network and either adds or deletes a percentage of
#' edges randomly
#'
#' @param network Network to modify.
#' @param mode Whether to add ('add') or delete ('del') edges.
#' @param perc Percentage of edges to modify per regulon.
#' @param tf Network's TF column.
#' @param target Network's target genes column.
#' @param seed An integer to set the RNG state for random number generation. Use
#'    NULL for random number generation.
#' @return Modified network.
#'
#' @import dplyr
#' @import purrr
#'
#' @export
net_noise <- function(network, mode='add', perc=0.1, tf='tf',
                          target='target', seed=42){
  # Get unique targets of network
  targets <- unique(network[[target]])

  # For each TF, change % of edges
  network %>%
    group_by(tf) %>%
    group_split() %>%
    map(function(df){
      # Name TF
      tf <- df[[tf]][1]

      # Number of genes
      n_genes <- ceiling(perc * nrow(df))

      if (mode == 'add'){
        # Get genes not in the regulon
        diff <- setdiff(targets, df[[target]])

        # Sample n edges to add
        set.seed(seed)
        sampled <- sample(diff, n_genes, replace=F)

        # Append to df
        tibble(tf=tf, target=sampled, mor=1, likelihood=1) %>%
          bind_rows(df, .)
      } else if (mode == 'del') {
        # Sample n edges to del
        set.seed(seed)
        sampled <- sample(df[[target]], n_genes, replace=F)

        # Delete rows with sampled genes
        df %>%
          filter(!target %in% sampled)
      } else {
        stop(stringr::str_glue(
          "{mode} is an invalid mode, please select 'add' or 'del'."))
      }
    }) %>%
    bind_rows()
}

