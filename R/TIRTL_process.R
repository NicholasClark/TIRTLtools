#' Run data processing functions on a TIRTLseqData object
#'
#' @description
#' `r lifecycle::badge('experimental')`
#' This function runs annotation and data cleaning functions on a
#' TIRTLseqData object. Specifically, it calls the functions \code{add_single_chain_data()},
#' \code{identify_paired()}, \code{identify_non_functional_seqs()}, and (optionally) \code{clean_pairs()}.
#'
#' @param data a TIRTLseqData object
#' @param clean (optional) a TRUE/FALSE value, whether or not to "clean" the paired data by removing
#' excess pairs for individual alpha and beta chains (default is FALSE).
#'
#' @family data_processing
#'
TIRTL_process = function(data, clean = FALSE) {
  if(clean) data = data %>% clean_pairs()
  data = data %>%
    add_single_chain_data() %>%
    identify_paired() %>%
    identify_non_functional_seqs()
  return(data)
}
