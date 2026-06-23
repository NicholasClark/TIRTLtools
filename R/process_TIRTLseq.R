#' Run data processing functions on a TIRTLseqDataSet object
#'
#' @description
#' `r lifecycle::badge('deprecated')`
#' This function has been deprecated because its functionality has been added to the load_tirtlseq() function
#'
#' This function runs annotation and data cleaning functions on a
#' TIRTLseqDataSet object. Specifically, it calls the functions \code{add_single_chain_data()},
#' \code{identify_paired()}, \code{identify_non_functional_seqs()}, and (optionally) \code{clean_pairs()}.
#'
#' @param data a TIRTLseqDataSet object
#' @param clean (optional) a TRUE/FALSE value, whether or not to "clean" the paired data by removing
#' excess pairs for individual alpha and beta chains (default is FALSE).
#' @param remove_nonfunctional whether to remove non-functional TCR chains (default is FALSE)
#'
#' @return a TIRTLseqDataSet object with annotated and (optionally) cleaned data
#'
#' @family deprecated
#'
process_TIRTLseq = function(data, clean = FALSE, remove_nonfunctional = FALSE) {

  lifecycle::deprecate_warn(when = "0.2.1", what="process_TIRTLseq()",
                 details = "This function has been deprecated because its functionality has been added to the load_tirtlseq() function"
                 )
  return(data)
  # sample_names = names(data$data)
  # for(i in 1:length(data$data)) {
  #   data$data[[i]] = .process_TIRTLseq_single(data$data[[i]], clean = clean, remove_nonfunctional = remove_nonfunctional)
  # }
  # names(data$data) = sample_names
  # return(data)
}
