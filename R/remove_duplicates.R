#' Removes duplicate paired TCRs (by default, paired TCRs are listed twice if they
#' were called by both pairing algorithms)
#'
#' @description There will be duplicates in paired chain data when a pair is called by both
#' the T-SHELL and MAD-HYPE algorithms.
#' The input 'data' needs to be either a single data frame (paired chain) or a list
#' of data frames (paired chain)
#'
#' @family data_processing
remove_duplicates = function(data) {
  if(is.data.frame(data)) {
    out = data[!duplicated(data$alpha_beta),]
  } else {
    out = lapply(data, function(x) {
      x[!duplicated(x$alpha_beta),]
    })
  }
  return(out)
}
