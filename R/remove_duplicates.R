#' Removes duplicate paired TCRs
#'
#' @description By default, paired TCRs are listed twice in the paired chain data
#' if they were called by both the T-SHELL and MAD-HYPE pairing algorithms. This
#' function removes duplicate TCRs and returns a data frame with only one of each pair.
#'
#'
#' @param data either a single data frame (paired chain) or a list
#' of data frames
#'
#'
#' @family data_processing
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
