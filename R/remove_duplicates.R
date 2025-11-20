#' Removes duplicate paired TCRs
#'
#' @description
#' `r lifecycle::badge('experimental')`
#' By default, paired TCRs are listed in one row in the paired chain data
#' for each pairing algorithm they were called by (either T-SHELL or MAD-HYPE). This means
#' some TCRs are present in two rows, once for each algorithm. This
#' function removes duplicate TCRs and returns a data frame where each TCR pair is
#' found in only one row.
#'
#' @param data either a single data frame (paired chain) or a list
#' of data frames
#'
#' @returns Returns a data frame (or list of data frames, depending on input) of paired TCRs
#' where each TCR is listed only once.
#'
#' @family data_processing
#' @export
#'
remove_duplicates = function(data) {
  if(.is_df(data)) {
    out = data[!duplicated(data$alpha_beta),]
  } else {
    out = lapply(data, function(x) {
      x[!duplicated(x$alpha_beta),]
    })
  }
  return(out)
}
