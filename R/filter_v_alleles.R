#' Remove TCRs with unknown V-segments
#'
#'
#' @description
#' `r lifecycle::badge('experimental')`
#'
#' @param df data frame with paired TCRs
#' @param params a table with acceptable V-segments (default is `TIRTLtools::params`)
#' @param verbose whether to print number of TCRs removed
#'
#' @family data_processing
#'
filter_v_alleles = function(df, params=NULL, verbose = TRUE) {
  if(is.null(params)) params = TIRTLtools::params
  df = .check_v_alleles(df, params = params)
  n_tcr_orig = nrow(df)
  n_not_allowed = sum(!df$va_and_vb_allowed)
  pct_not_allowed = signif(100*n_not_allowed/n_tcr_orig,2)
  msg = paste("Removed ", n_not_allowed %>% .add_commas(), " TCRs with unknown V-segments ", "(", pct_not_allowed, "%) from a total of ", n_tcr_orig %>% .add_commas(), " TCRs.", sep = "")
  if(verbose) message(msg)
  df = df %>% filter(va_and_vb_allowed) %>% select(-va_and_vb_allowed, -va_allowed, -vb_allowed)
  return(df)
}

### check that va and vb are found in the parameters for tcrdist (but don't remove them)
.check_v_alleles = function(df, params=NULL) {
  if(is.null(params)) params = TIRTLtools::params
  df = df %>% mutate(va_allowed = va %in% params$feature, vb_allowed = vb %in% params$feature) %>%
    mutate(va_and_vb_allowed = va_allowed & vb_allowed)
  return(df)
}
