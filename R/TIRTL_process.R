TIRTL_process = function(data) {
  data = data %>%
    add_single_chain_data() %>%
    clean_pairs() %>%
    identify_paired()
  return(data)
}
