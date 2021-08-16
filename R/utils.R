# Helper functions for run_benchmark

#' Helper Function to to generate the bools used to check if the current
#' locations/rds objects are the same as the previous one.
#'
#' @param .design input tibble used to provide the experimental design for each
#' benchmark run
#'
#' @keywords internal
#'
#' @details This is used to limit the number of times that any of the
#' prerequsites is loaded.
#' @export
format_design <- function(.design){
  if (!"noise_crit" %in% colnames(.design)){
    .design <- .design %>%
      add_column('noise_crit' = NA)
  }

  if (!"weight_crit" %in% colnames(.design)){
    .design <- .design %>%
      add_column('weight_crit' = NA)
  }

  .design %>%
    mutate(.source_bln = .data$source_loc %>% check_preced(),
           .expr_bln = .data$bexpr_loc %>% check_preced(),
           .meta_bln = .data$bmeta_loc %>% check_preced())
}


#' Helper Function that checks if the preceding vector element is the same
#' as the current element
#'
#' @param vector_loc character vector with directory paths
#'
#' @return logical values describing whether the location of the loaded files
#' has changes
#'
#' @keywords internal
check_preced <- function(vector_loc){
  tib_loc <- tibble(current=vector_loc, behind=lag(vector_loc))

  pmap_lgl(tib_loc, function(behind, current){
    ifelse(is.na(behind) || behind!=current, FALSE, TRUE)
  })
}


#' Helper Function to filter and format the gene set resource
#'
#' @param set_source Set Source (e.g. TF regulon sets, GO:term sets, etc)
#' @inheritParams input_tibble
#' @param .minsize minimum size of each set
#' @param .silent bool whether to silence wanring messages
#'
#' @importFrom stringr str_glue
#'
#' @return returns a filtered and formatted set source
#'
#' @details Filtering can be omitted if `filter_col` is `NA`.
filter_sets <- function(set_source,
                        source_col,
                        filter_col,
                        filter_crit,
                        .minsize,
                        .silent){
  n_duprows <- sum(duplicated(set_source))
  na_bool <- is.na(filter_col)

  gs_filtered <- set_source %>%
    {
      if(na_bool){distinct(.)}
      else if(!na_bool){
        filter(., .data[[filter_col]] %in% filter_crit) %>%
          distinct_at(vars(-.data[[filter_col]]), .keep_all = FALSE)
      }
    } %>%
    group_by(.data[[source_col]]) %>%
    add_count() %>%
    filter(n >= .minsize) %>%
    ungroup()

  if (n_duprows & !.silent){
    warning(str_glue("{n_duprows} rows were duplicated in the set resource! ",
                     "{sum(duplicated(gs_filtered))} duplicated rows ",
                     "remain after filtering."))
  }
  return(gs_filtered)
}


#' `base::readRDS` helper function that enables loading files from urls
#'
#' @inheritParams base::readRDS
#' @inheritDotParams base::readRDS
#'
#' @param .url_bool bool whether the location is a url or not
#'
#' @export
readRDS_helper <- function(file, .url_bool=FALSE, ...){
  if(.url_bool){
    readRDS(url(file, "rb", ...))
  } else{
    readRDS(file, ...)
  }
}

#' Function to format benchmarking results
#'
#' @param .bench_res benchmarking results
#' @param .silent bool whether to silence warnings or not
#' @returns formatted benchmarking results
#' @importFrom rlang .data
#' @importFrom stringr str_glue_data
#'
#' @details If infinite values are present in the results, this function will
#'   notify the user.
bench_format <- function(.bench_res, .silent) {
  res_format <- .bench_res %>%
    unnest(.data$activity) %>%
    # convert filter_criteria from character to string
    rowwise() %>%
    mutate(filter_crit = paste0(unlist(.data$filter_crit), collapse = "")) %>%
    ungroup() %>%
    # get statistic name
    mutate(statistic = .data$activity %>%
             map(function(tib)
               unique(tib[["statistic"]]))) %>%
    unnest(.data$statistic) %>%
    select(.data$set_name, .data$bench_name, .data$filter_crit,
           .data$statistic, .data$activity)

  inf_sums <- res_format$activity %>%
    map(function(x) sum(is.infinite(x$score))) %>%
    setNames(
      paste(
        res_format$set_name,
        res_format$bench_name,
        res_format$statistic,
        sep = "_"
      )) %>%
    enframe() %>% unnest(.data$value)

  if (sum(inf_sums$value)) {
    res_format <- res_format %>%
      mutate(activity = .data$activity %>%
               map(function(tib)
                 tib %>%
                   mutate_at(
                     vars(.data$score), ~ replace(., is.infinite(.), 0)
                   )))
    warning(
      inf_sums %>%
        filter(.data$value > 0) %>%
        str_glue_data("{.$value} infinite values were filtered",
                      " in {.$name}. \n ")
    )

  }
  return(res_format)
}


# benchmark input --------------------------------------------------------------
#' Benchmark input tibble containing the (experimental) design for each of the
#' benchmark runs corresponding to rows
#' @name input_tibble
#'
#' @details A tibble with locations, options, and filter options for
#' the desired benchmark setting
#'
#' @param set_name user-defined name of the set resource
#' @param bench_name user-defined name of the benchmark data
#' @param stats_list List of statistics to run
#' @param opts_list Named list containing the options for each stat. method
#' @param bexpr_loc benchmark expression data location (.rds format tibble)
#' @param bmeta_loc benchmark metadata location (.rds format tibble)
#' @param source_loc set source (e.g. network resource, gene ontology sets,
#'  kinase sets, etc.) location (.rds format tibble)
#' @param source_col name of the column with the source for the set source
#' @param target_col name of the column with the targets for the set source
#' @param filter_col name of the column by which we wish to filter
#' @param filter_crit criteria by which we wish to filter the `filter_col`
NULL
