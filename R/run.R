#' Run decouple on a design matrix.
#'
#' @param .design design tibble.
#' @param .minsize Minsize to filter
#' @param .silent Whether to run without verbose.
#' @param .url_bool Wheter to download data.
#'
#' @import tibble tidyr dplyr tidyselect purrr
#'
#' @export
#'
#' @importFrom rlang .data
#' @importFrom stats reorder setNames
#' @importFrom methods new
#'
#' @return Tibble of activities.
run_methods <- function(bench_env,
                        .design,
                        .minsize = 5,
                        .silent = TRUE,
                        .url_bool = FALSE
){
  res <- .design %>%
    decoupleRBench::format_design() %>%
    dplyr::mutate(activity = pmap(.l=.,
                                  .f=function(set_name, bench_name,
                                              stats_list, opts_list,
                                              bexpr_loc, bmeta_loc, source_loc,
                                              source_col, target_col, filter_col,
                                              filter_crit, noise_crit, weight_crit,
                                              extra_name, consensus_crit, .source_bln, .expr_bln, .meta_bln){

                                    # Check prerequisites
                                    if(!.expr_bln){
                                      bench_env$gene_expression <- decoupleRBench::readRDS_helper(bexpr_loc, .url_bool) %>%
                                        as.matrix()
                                      message(stringr::str_glue("Expression loaded"))
                                    }
                                    if(!.meta_bln){
                                      bench_env$meta_data <- decoupleRBench::readRDS_helper(bmeta_loc, .url_bool)
                                      message(stringr::str_glue("Meta loaded"))
                                    }

                                    # Filter set_source/network
                                    bench_env$set_source <- decoupleRBench::readRDS_helper(source_loc, .url_bool)
                                    ss_filtered <- filter_sets(bench_env$set_source, source_col,
                                                               filter_col, filter_crit,
                                                               .minsize, .silent)
                                    message(stringr::str_glue("Network loaded"))

                                    # Add noise
                                    if (is.list(noise_crit)){
                                      message(
                                        stringr::str_glue("Modify network"))
                                      ss_filtered <- decoupleRBench::net_noise(
                                        network = ss_filtered,
                                        mode = noise_crit$mode,
                                        perc = noise_crit$perc,
                                        seed = noise_crit$seed,
                                        source = source_col,
                                        target = target_col
                                      )
                                      message(
                                        stringr::str_glue("{noise_crit$mode} {noise_crit$perc} noise"))
                                    }

                                    # Remove weight
                                    if (is.list(weight_crit)){
                                      message(
                                        stringr::str_glue("Unweight network"))
                                      ss_filtered <- decoupleRBench::net_weight(
                                        network = ss_filtered,
                                        weight_rm = weight_crit$.mor
                                      )
                                      message(
                                        stringr::str_glue("{weight_crit$.mor} set to 1 "))
                                    }
                                    consensus_crit <- consensus_crit[[1]]
                                    if (!any(is.na(consensus_crit))){
                                      message(
                                        stringr::str_glue("Setting consensus criteria {consensus_crit}"))
                                      if (consensus_crit == 'None'){
                                        consensus_score <- FALSE
                                        consensus_stats <- NULL
                                      } else{
                                        consensus_score <- TRUE
                                        consensus_stats <- consensus_crit
                                      }
                                    } else {
                                      consensus_score <- TRUE
                                      consensus_stats <- NULL
                                    }

                                    # Show Current Row/Run
                                    if(!.silent){
                                      .curr_row <- paste(set_name, bench_name,
                                                         paste0(unlist(filter_crit), collapse=""),
                                                         sep="_")
                                      message(str_glue("Currently Running: {.curr_row}"))
                                    }

                                    # Obtain Activity with decouple and format
                                    decoupleR::decouple(mat = bench_env$gene_expression, network = ss_filtered,
                                                        .source = source_col, .target = tidyselect::all_of(target_col),
                                                        statistics = stats_list, args = opts_list, consensus_score=consensus_score,
                                                        consensus_stats = consensus_stats, include_time = TRUE)
                                  }))
  return(res)
}

#' Run evaluation of results from decouple.
#' @param .design design tibble.
#' @param  res results from decouple.
#' @param .form bool whether to format or not
#' @param .perform bool whether to calculate ROC and performance summary
#' @param .downsample_pr whether to downsample precision recall curve TNs
#' @param .downsample_roc whether to downsample ROC true negatives
#' @param .downsample_times downsampling iterations
#' @param .silent Whether to run without verbose.
#'
#' @import tibble tidyr dplyr tidyselect purrr
#'
#' @export
#'
#' @importFrom rlang .data
#' @importFrom stats reorder setNames
#' @importFrom tidyr unnest
#' @importFrom methods new
#'
#' @return An S4 object of \link{BenchResult-class}
run_eval <- function(bench_env,
                     .design,
                     res,
                     .form = TRUE,
                     .perform = TRUE,
                     .downsample_pr = FALSE,
                     .downsample_roc = FALSE,
                     .downsample_times = 100,
                     .silent = F
){
  res <- res %>%
    dplyr::mutate(activity = map(activity, function(df){
      df %>%
        dplyr::rename(id=.data$condition) %>%
        ungroup() %>%
        dplyr::left_join(bench_env$meta_data, by="id")  %>%
        dplyr::group_by(.data$statistic) %>%
        dplyr::group_split(.keep=TRUE) %>%
        as.list()
    })) %>%
    {
      if(.form & !.perform) bench_format(., .silent)
      else if(.form & .perform) bench_format(., .silent) %>%
        mutate(roc = .data$activity %>%
                 map(~calc_curve(df=.x,
                                 downsampling=.downsample_roc,
                                 times=.downsample_times,
                                 curve="ROC")),
               prc = .data$activity %>%
                 map(~calc_curve(df=.x,
                                 downsampling=.downsample_pr,
                                 times=.downsample_times,
                                 curve="PR")))
      else .
    }

  if(.form & .perform){
    bench_result <- methods::new("BenchResult",
                                 bench_res=res,
                                 summary=res %>% get_bench_summary(),
                                 design=.design)
  }
  else{
    bench_result <- methods::new("BenchResult",
                                 bench_res=res,
                                 summary=list(NULL),
                                 design=.design)
  }
  return(bench_result)
}


#' Benchmark pipeline built on the statistical method wrapper \link{decouple}.
#'
#' @param .form bool whether to format or not
#' @param .perform bool whether to calculate ROC and performance summary
#' @inheritParams filter_sets
#' @param .downsample_pr whether to downsample precision recall curve TNs
#' @param .downsample_roc whether to downsample ROC true negatives
#' @param .downsample_times downsampling iterations
#' @inheritParams readRDS_helper
#'
#' @import tibble tidyr dplyr tidyselect purrr
#'
#' @seealso See \link{input_tibble} for a description of the params/columns
#'   of .design (i.e. input tibble).
#'
#' @export
#'
#' @importFrom rlang .data
#' @importFrom stats reorder setNames
#' @importFrom methods new
#'
#' @return An S4 object of \link{BenchResult-class}
run_benchmark <- function(.design,
                          .form = TRUE,
                          .perform = TRUE,
                          .minsize = 5,
                          .silent = TRUE,
                          .downsample_pr = FALSE,
                          .downsample_roc = FALSE,
                          .downsample_times = 100,
                          .url_bool = FALSE
){
  bench_env <- new.env()
  res <- run_methods(bench_env, .design,.minsize,.silent,.url_bool)
  # Merge rows with same set_name and bench_name
  res <- res %>%
    group_by(set_name, bench_name) %>%
    group_modify(., function(df, ..., .keep=TRUE){
      if(nrow(df)>1){
        first_row <- df[1,]
        extra_name <- df[2,'extra_name']
        for (i in 2:nrow(df)) {
          row <- df[i,] %>% select(activity) %>% unnest(activity) %>% mutate(statistic=paste0(statistic,'_',extra_name))
          first_row <- first_row %>%
            mutate(activity=map(activity, function(act){
              bind_rows(act, row)
            }))
        }
        return(first_row)
      } else {
        return(df)
      }
    }) %>% ungroup()
  bench_result <- run_eval(bench_env, .design,res,.form,.perform,.downsample_pr,.downsample_roc,.downsample_times)
  return(bench_result)
}
