plot_n_reads = function(data, chain = c("both","beta", "alpha"), samples = NULL, color_scheme = NULL) {
  chain = chain[1]
  df = summarize_data(data)
  if(!is.null(samples)) df = df %>% filter(sample_id %in% samples)
  if(chain == "both") {
    df1 = df %>% select(n_reads_alpha, sample_id) %>% dplyr::rename(n_reads = n_reads_alpha) %>% mutate(chain = "alpha")
    df2 = df %>% select(n_reads_beta, sample_id) %>% dplyr::rename(n_reads = n_reads_beta) %>% mutate(chain = "beta")
    df = bind_rows(df1, df2)
    n_reads_col = "n_reads"
  } else {
    n_reads_col = paste("n_reads_", chain, sep = "")
  }
    brks = .get_log_labels_pos(df[[n_reads_col]])
    gg = ggplot(df) + geom_col(aes(x=sample_id, y=!!sym(n_reads_col), fill = !!sym(n_reads_col))) +
      scale_y_log10(breaks = brks$brks, labels = brks$labels) +
      theme_classic() +
      scale_fill_gradientn(colors = .tirtl_colors_gradient(palette=color_scheme)) +
      .rotate_x_labels()
  if(chain == "both") gg = gg + facet_wrap(~chain)
  return(gg)
}
