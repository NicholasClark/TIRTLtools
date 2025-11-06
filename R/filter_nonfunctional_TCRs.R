#' Remove TCRs with nonfunctional CDR3 amino acid sequences
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
filter_nonfunctional_TCRs = function(df, verbose = TRUE) {
  df = .identify_non_functional_seqs_single(df)
  n_tcr_orig = nrow(df)
  n_nonfunc = sum(!df$is_functional)
  pct_nonfunc = signif(100*n_nonfunc/n_tcr_orig,2)
  msg = paste("Removed ", n_nonfunc %>% .add_commas(), " TCRs with non-functional CDR3 amino acid sequences ", "(", pct_nonfunc, "%) from a total of ", n_tcr_orig %>% .add_commas(), " TCRs.", sep = "")
  if(verbose) message(msg)
  df = df %>% filter(is_functional) %>%
    select(-alpha_has_stop_codon, -alpha_has_frameshift, -beta_has_stop_codon, -beta_has_frameshift,
           -alpha_is_functional, -beta_is_functional, -is_functional)
  return(df)
}
