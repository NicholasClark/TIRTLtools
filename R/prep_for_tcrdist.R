#' Prepare paired TCRs for TCRdist calculation
#'
#' @description This function prepares a paired TCR data frame for running TCRdist,
#' adding alleles to V-segments (i.e. "*01") if they are missing and removing
#' TCRs with unrecognized V-segments or non-functional CDR3 sequences.
#'
#' @param df a data frame of paired TCRs
#' @param params (optional) a data frame with permissible amino acids and V-segments.
#' By default, the TIRTLtools::params data frame is used.
#' @param remove_MAIT (optional) whether or not to remove MAIT cells from the data frame (default is TRUE).
#' @param verbose (optional) whether to print progress of the function (default is TRUE).
#'
#' @family data_processing
#'
prep_for_tcrdist = function(df, params=NULL, remove_MAIT = TRUE, verbose = TRUE) {
  if(is.null(df)) return(df)
  if(is.null(params)) params = TIRTLtools::params
  df = .add_alleles(df) # add "*01" as allele for va and vb if necessary
  df = .filter_alleles(df, params = params) # remove alleles not found in the parameter data frame
  df = df %>% filter(nchar(cdr3a) > 5, nchar(cdr3b) > 5)
  #if(!"is_functional" %in% colnames(df))
  df = .identify_non_functional_seqs_single(df)
  df = df %>% filter(is_functional) # remove seqs w/ stop codons (*) or frameshifts (_)
  df = as.data.frame(df) ## convert to standard data frame
  if(remove_MAIT) df = filter_mait(df, verbose = verbose)
  return(df)
}
