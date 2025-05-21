identify_paired = function(data) {
  data$data = lapply(data$data, function(x) identify_paired_single(x))
  return(data)
}


identify_paired_single = function(data) {
  paired_tmp = data$paired %>% remove_dupes_paired()
  counts_alpha = paired_tmp[, .N, by = alpha_nuc] %>% dplyr::rename(n_paired = N, targetSequences = alpha_nuc)
  counts_beta = paired_tmp[, .N, by = beta_nuc] %>% dplyr::rename(n_paired = N, targetSequences = beta_nuc)
  data$alpha = left_join(data$alpha, counts_alpha) %>%
    mutate(n_paired = ifelse(is.na(n_paired), 0, n_paired)) %>%
    mutate(is_paired = n_paired != 0)
  data$beta = left_join(data$beta, counts_beta) %>%
    mutate(n_paired = ifelse(is.na(n_paired), 0, n_paired)) %>%
    mutate(is_paired = n_paired != 0)
  return(data)
}
