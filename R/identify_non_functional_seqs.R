#' Identify TCRs that contain non-functional CDR3 sequences
#'
#' @description
#' `r lifecycle::badge('deprecated')`
#'
#' This function has been deprecated because its functionality has been added to the load_tirtlseq() function
#'
#' \code{identify_non_functional_seqs()} identifies CDR3 nucleotide sequences in
#' the data that contain either stop codons (*) or frame shifts (_) that
#' would indicate a non-functional protein product.
#'
#'
#' @param data a TIRTLseqDataSet object
#'
#' @return
#' A TIRTLseqDataSet object with modified pseudobulk and paired data frames for each
#' sample. Each dataframe in the ouptut object has added columns that identify
#' whether the CDR3 alpha and beta
#' nucleotide sequences contain any stop codons (*) or frame shifts (_) that would
#' indicate a non-functional chain.
#'
#' Column \code{is_functional} is TRUE if neither chain has a stop codon or a frame shift.
#' The columns \code{has_stop_codon} and \code{has_frameshift} are similar, but specific to
#' each kind of coding error. Other columns identify if the alpha chain or beta chain has
#' a stop codon or frameshift, and if it is functional.
#'
#' @family deprecated
#'
#' @seealso \code{\link{TCRdist}()}, \code{\link{plot_num_partners}()}
#'
#' @export
#' @examples
#' # example code
#'
#'

identify_non_functional_seqs = function(data) {
  lifecycle::deprecate_warn(when = "0.2.1", what="identify_non_functional_seqs()",
                            details = "This function has been deprecated because its functionality has been added to the load_tirtlseq() function"
  )
  return(data)
  # is_list = .is.list.only(data)
  # if(is_list) {
  #   data_out = lapply(data$data, function(x) {
  #     lapply(x, .identify_non_functional_seqs_single)
  #   })
  #   data$data = data_out
  #   return(data)
  # } else {
  #   return(.identify_non_functional_seqs_single(data))
  # }
}

