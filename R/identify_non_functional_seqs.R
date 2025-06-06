#' @title
#' Identify non-functional CDR3 nucleotide sequences
#'
#' @description
#' \code{identify_non_functional_seqs()} identifies CDR3 nucleotide sequences in
#' the data that contain errors, either stop codons (*) or frame shifts (_) that
#' would indicate a non-functional chain.
#'
#'
#' @param df a dataframe of TIRTLseq data, either paired or single-chain pseudobulk for one sample
#'
#' @return
#' A dataframe with added columns that identify whether the CDR3 alpha and beta
#' nucleotide sequences contain any stop codons (*) or frame shifts (_) that would
#' indicate a non-functional chain.
#'
#' If column \code{is_functional} is TRUE if neither chain has a stop codon or a frame shift.
#' The columns \code{has_stop_codon} and \code{has_frameshift} are similar, but specific to
#' each kind of coding error. Other columns identify if the alpha chain or beta chain has
#' a stop codon or frameshift, and if it is functional.
#'
#' @seealso \code{\link{TCRdist}()}, \code{\link{plot_num_partners}()}
#'
#' @export
#' @examples
#' # example code
#'
#'

identify_non_functional_seqs = function(df) {
  if("cdr3a" %in% colnames(df) && "cdr3b" %in% colnames(df)) {
    ### for paired data
    df_out = df %>%
      mutate(
        alpha_has_stop_codon = grepl("\\*", cdr3a), # stop codon
        alpha_has_frameshift = grepl("_", cdr3a), # frame shift
        beta_has_stop_codon = grepl("\\*", cdr3b),
        beta_has_frameshift = grepl("_", cdr3b),
      ) %>%
      mutate(
        alpha_is_functional = !(alpha_has_stop_codon | alpha_has_frameshift),
        beta_is_functional = !(beta_has_stop_codon | beta_has_frameshift)
      ) %>%
      mutate(is_functional = alpha_is_functional & beta_is_functional)
  } else {
    ### for pseudobulk data
    df_out = df %>%
      mutate(
        has_stop_codon = grepl("\\*", aaSeqCDR3), # stop codon
        has_frameshift = grepl("_", aaSeqCDR3), # frame shift
      ) %>%
      mutate(is_functional = !(has_stop_codon & has_frameshift))
  }
  return(df_out)
}
