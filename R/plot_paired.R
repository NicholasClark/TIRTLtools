#' Stacked bar plot of the number alpha/beta chains paired by each pairing algorithm
#'
#' @family plotting
#'
plot_paired = function(df, chain = c("paired", "alpha", "beta"), samples = NULL, color_scheme = NULL) {
  if(.is.list.only(df)) df = get_pair_stats(df, verbose = FALSE)
  chain1 = chain[1]
  df$log10Freq = log10(df$Freq)
  df_sub = df %>% filter(chain == chain1)
  if(!is.null(samples)) df_sub = df_sub %>% filter(sample_id %in% samples)
  #y_col = sym("Freq")
  #if(log_scale) y_col = sym("log10Freq")
  gg = ggplot(df_sub) + geom_col(aes(x=sample_id, y=Freq, fill = category)) + theme_classic() +
    scale_fill_manual(values = .tirtl_colors_gradient(palette=color_scheme, n=length(unique(df_sub$category)) )) +
    theme(axis.text.x = element_text(angle = 90, vjust = 0.5))
  return(gg)
}
