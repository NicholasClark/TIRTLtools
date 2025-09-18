#' Bar plot of the fraction of single chains that were paired within different
#' read fraction ranges for each sample.
#'
#' @family plotting
plot_paired_by_read_fraction_range = function(df, chain = c("beta","alpha"), cutoffs = 10^(-6:-1), freq=FALSE, samples = NULL, color_scheme = NULL) {
  chain = chain[1]
  if(.is.list.only(df)) df = get_paired_by_read_fraction_range(df, chain, cutoffs)
  #gg_df = get_paired_by_read_fraction_range(data)
  #gg_df = left_join(gg_df, data$meta)
  if(!is.null(samples)) df = df %>% filter(sample_id %in% samples)
  n=length(unique(df$range))
  if(!freq) {
    gg = ggplot(df) + geom_col(aes(x=sample_id, y=fraction_paired, fill = range), position = "dodge") + theme_classic() +
      scale_fill_manual(values = .tirtl_colors_gradient(palette=color_scheme, n=n)) +
      .rotate_x_labels()
  } else {
    gg = ggplot(df) + geom_col(aes(x=sample_id, y=n_paired_sum, fill = range), position = "dodge") + theme_classic() +
      scale_fill_manual(values = .tirtl_colors_gradient(palette=color_scheme, n=n)) +
      .rotate_x_labels() +
      scale_y_log10()
  }
  return(gg)
}
