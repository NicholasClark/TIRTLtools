#' Remove TCRs with short CDR3 loops
#'
#' @description
#' `r lifecycle::badge('experimental')`
#' This function removes TCR pairs where the CDR3 amino acid segment is short.
#' By defualt, it removes TCRs where CDR3-alpha or CDR3-beta is less than 6 amino acids long.
#'
#' @param df data frame with paired TCRs
#' @param verbose whether to print number of TCRs removed
#' @param min_aa the minimum number of amino acids allowed for a CDR3 sequence
#'
#' @family data_processing
#' @export
#'
filter_short_cdr3s = function(df, verbose = TRUE, min_aa = 6) {
  df = df %>% mutate(is_short_cdr3 = nchar(cdr3a) < min_aa | nchar(cdr3b) < min_aa)
  n_tcr_orig = nrow(df)
  n_short = sum(df$is_short_cdr3)
  pct_short = signif(100*n_short/n_tcr_orig,2)
  msg = paste("Removed ", n_short %>% .add_commas(), " TCRs with short CDR3 segments ", "(", pct_short, "%) from a total of ", n_tcr_orig %>% .add_commas(), " TCRs.", sep = "")
  if(verbose) message(msg)
  df = df %>% filter(!is_short_cdr3) %>% select(-is_short_cdr3)
  return(df)
}
