#' This function takes a network and removes weight and/or directionality of the network
#'
#' @param network Network to modify.
#' @param weight_rm Columns from network specifying weights ('mor').
#' @return Modified network.
#'
#' @import dplyr
#'
#' @export
net_weight <- function(network, weight_rm='mor'){
 network %>% mutate_at(weight_rm, function(x) replace(x, TRUE, 1))
}

