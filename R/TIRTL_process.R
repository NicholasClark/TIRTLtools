#' Run the other data processing functions in a sequence (add_single_chain_data, clean_pairs, and identify_paired)
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
