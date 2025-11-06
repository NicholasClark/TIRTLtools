#' Prepare paired TCRs for TCRdist calculation
#'
#' @description
#' `r lifecycle::badge('experimental')`
#' This function prepares a paired TCR data frame for running TCRdist,
#' adding alleles to V-segments (i.e. "*01") if they are missing and removing
#' TCRs with unrecognized V-segments or non-functional CDR3 sequences.
#'
#' @param df a data frame of paired TCRs
#' @param params (optional) a data frame with permissible amino acids and V-segments.
#' By default, the TIRTLtools::params data frame is used.
#' @param remove_duplicate_tcrs (optional) whether or not to remove duplicate TCRs (default is FALSE, do not remove duplicates).
#' @param remove_MAIT (optional) whether or not to remove MAIT cells from the data frame (default is TRUE).
#' @param verbose (optional) whether to print progress of the function (default is TRUE).
#'
#' @family data_processing
#'
prep_for_tcrdist = function(df, params=NULL, remove_duplicate_tcrs = FALSE, remove_MAIT = TRUE, verbose = TRUE) {
  if(is.null(df)) return(df)
  if(is.null(params)) params = TIRTLtools::params
  n_tcrs_orig = nrow(df)
  df = .add_alleles(df) # add "*01" as allele for va and vb if necessary
  df = filter_v_alleles(df, params = params, verbose = verbose) # remove TCRs with alleles not found in the parameter data frame
  df = filter_short_cdr3s(df, verbose = verbose) ## remove TCRs with short CDR3s
  df = filter_nonfunctional_TCRs(df, verbose = verbose)
  if(remove_MAIT) df = filter_mait(df, verbose = verbose)
  if(remove_duplicate_tcrs) df = filter_duplicate_tcrs(df, verbose = verbose)
  df = as.data.frame(df) ## convert to standard data frame
  n_tcrs_final = nrow(df)
  pct_tcrs = signif(100*n_tcrs_final/n_tcrs_orig, 2)
  msg = paste("Filtered data frame contains ", n_tcrs_final %>% .add_commas(), " TCRs ", "(", pct_tcrs, "%) of original ", n_tcrs_orig  %>% .add_commas(), " TCRs.", sep = "")
  if(verbose) message(msg)
  return(df)
}
