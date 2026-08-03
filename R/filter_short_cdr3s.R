#' Remove TCRs with short CDR3 loops
#'
#' @description
#' `r lifecycle::badge('experimental')`
#'
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
  use_alpha = ifelse("cdr3a" %in% colnames(df), TRUE, FALSE)
  use_beta = ifelse("cdr3b" %in% colnames(df), TRUE, FALSE)
  if(use_alpha) {
    df = df %>% mutate(is_short_cdr3a = nchar(cdr3a) < min_aa)
  }
  if(use_beta) {
    df = df %>% mutate(is_short_cdr3b = nchar(cdr3b) < min_aa)
  }
  if(use_alpha && use_beta) {
    df = df %>% mutate(is_short_cdr3 = is_short_cdr3a & is_short_cdr3a)
  } else{
    if(use_alpha) df = df %>% mutate(is_short_cdr3 = is_short_cdr3a)
    if(use_beta) df = df %>% mutate(is_short_cdr3 = is_short_cdr3b)
  }

  n_tcr_orig = nrow(df)
  n_short = sum(df$is_short_cdr3)
  pct_short = signif(100*n_short/n_tcr_orig,2)
  msg = paste("Removed ", n_short %>% .add_commas(), " TCRs with short CDR3 segments ", "(", pct_short, "%) from a total of ", n_tcr_orig %>% .add_commas(), " TCRs.", sep = "")
  if(verbose) message(msg)
  df = df %>% filter(!is_short_cdr3) %>% select(-is_short_cdr3)
  return(df)
}
