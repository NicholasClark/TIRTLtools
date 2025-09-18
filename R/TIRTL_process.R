#' Run the other data processing functions in a sequence (add_single_chain_data, clean_pairs, and identify_paired)
#'
#' @family data_processing
#'
TIRTL_process = function(data) {
  data = data %>%
    add_single_chain_data() %>%
    clean_pairs() %>%
    identify_paired()
  return(data)
}
