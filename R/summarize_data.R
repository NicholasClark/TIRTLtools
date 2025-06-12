## input is an object loaded by load_tirtlseq()
summarize_data = function(data) {
  df = lapply(1:length(data$data), function(i) {
    x = data$data[[i]]
    name = names(data$data)[i]
    df_tmp = .summarize_data_single(x) %>% mutate(sample_id = name)
  }) %>% bind_rows()
  return(df)
}

## input is one sample -- list with 3 data frames: alpha, beta, and paired
.summarize_data_single = function(data) {
  n_alpha_clones = dim(data$alpha)[1]
  n_beta_clones = dim(data$beta)[1]
  df_paired_no_dupes = remove_dupes_paired(data$paired)
  n_pairs = dim(df_paired_no_dupes)[1]
  out = tibble(n_alpha = n_alpha_clones, n_beta = n_beta_clones, n_pairs = n_pairs)
  return(out)
}
