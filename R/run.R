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
                          .minsize = 10,
                          .silent = TRUE,
                          .downsample_pr = FALSE,
                          .downsample_roc = FALSE,
                          .downsample_times = 100,
                          .url_bool = FALSE
){

  bench_env <- new.env()

  res <- .design %>%
    decoupleRBench::format_design() %>%
    dplyr::mutate(activity = pmap(.l=.,
                           .f=function(set_name, bench_name,
                                       stats_list, opts_list,
                                       bexpr_loc, bmeta_loc, source_loc,
                                       source_col, target_col, filter_col,
                                       filter_crit, noise_crit, weight_crit,
                                       .source_bln, .expr_bln, .meta_bln){

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
                                 weight_rm = c(weight_crit$.mor, weight_crit$.likelihood)
                               )
                               message(
                                 stringr::str_glue("{c(weight_crit$.mor, weight_crit$.likelihood)} set to 1 "))
                             }

                             # Show Current Row/Run
                             if(!.silent){
                               .curr_row <- paste(set_name, bench_name,
                                                  paste0(unlist(filter_crit), collapse=""),
                                                  sep="_")
                               message(str_glue("Currently Running: {.curr_row}"))
                             }

                             # Match target genes between network and mat
                             targets <- rownames(bench_env$gene_expression)
                             msk <- ss_filtered[[target_col]] %in% targets
                             ss_filtered <- ss_filtered[msk,]
                             ss_filtered <- ss_filtered %>%
                               group_by_at(source_col) %>%
                               filter(n() >= .minsize)
                             # Obtain Activity with decouple and format
                             decoupleR::decouple(mat = bench_env$gene_expression, network = ss_filtered,
                                      .source = source_col, .target = tidyselect::all_of(target_col),
                                      statistics = stats_list, args = opts_list,
                                      include_time = TRUE)  %>%
                               dplyr::rename(id=.data$condition) %>%
                               ungroup() %>%
                               dplyr::left_join(bench_env$meta_data, by="id")  %>%
                               dplyr::group_by(.data$statistic) %>%
                               dplyr::group_split(.keep=TRUE) %>%
                               as.list()
                           })) %>% {
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
