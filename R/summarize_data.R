#' Create a summary table with number of reads and unique alpha/beta chains observed for each sample
#'
#' @param data a TIRTLseqData object
#'
#' @family repertoire_analysis
#'
summarize_data = function(data) {
  df = lapply(1:length(data$data), function(i) {
    x = data$data[[i]]
    name = names(data$data)[i]
    df_tmp = .summarize_data_single(x) %>% mutate(sample_id = name) %>% select(sample_id, everything())
  }) %>% bind_rows()
  return(df)
}

## input is one sample -- list with 3 data frames: alpha, beta, and paired
.summarize_data_single = function(data) {
  n_alpha_clones = dim(data$alpha)[1]
  n_beta_clones = dim(data$beta)[1]
  df_paired_no_dupes = remove_duplicates(data$paired)
  n_pairs = dim(df_paired_no_dupes)[1]
  n_reads_beta = sum(data$beta$readCount)
  n_reads_alpha = sum(data$alpha$readCount)
  out = tibble(n_alpha_clones = n_alpha_clones, n_beta_clones = n_beta_clones, n_clone_pairs = n_pairs, n_reads_alpha = n_reads_alpha, n_reads_beta = n_reads_beta)
  return(out)
}
