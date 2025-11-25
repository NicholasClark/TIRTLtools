#' Run data processing functions on a TIRTLseqDataSet object
#'
#' @description
#' `r lifecycle::badge('experimental')`
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
#' @family data_processing
#'
TIRTL_process = function(data, clean = FALSE, remove_nonfunctional = FALSE) {
  if(remove_nonfunctional) {
    data = data %>% filter_nonfunctional_TCRs()
  } else {
    data = data %>% identify_non_functional_seqs()
  }
  if(clean) data = data %>% clean_pairs()
  data = data %>%
    add_single_chain_data() %>%
    identify_paired()
  return(data)
}
