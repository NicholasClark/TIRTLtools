#' Remove duplicate TCRs from a data frame
#'
#'
#' @description
#' `r lifecycle::badge('experimental')`
#'
#' @param df data frame with paired TCRs
#' @param verbose whether to print number of TCRs removed
#'
#' @family data_processing
#'
filter_duplicate_tcrs = function(df, verbose = TRUE) {
  df = df %>% mutate(is_duplicate = duplicated(alpha_beta))
  n_tcr_orig = nrow(df)
  n_dupes = sum(df$is_duplicate)
  pct_dupes = signif(100*n_dupes/n_tcr_orig, 2)
  msg = paste("Removed ", n_dupes %>% .add_commas(), " duplicate TCRs ", "(", pct_dupes, "%) from a total of ", n_tcr_orig %>% .add_commas(), " TCRs.", sep = "")
  if(verbose) message(msg)
  df = df %>% filter(!is_duplicate) %>% select(-is_duplicate)
  return(df)
}
