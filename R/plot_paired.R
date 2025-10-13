#' Stacked bar plot of the number alpha/beta chains paired by each pairing algorithm
#'
#' @param data a TIRTLseqData object or a data frame created using \code{\link{get_pair_stats}()}
#' @param chain the TCR chain to plot (default is "paired")
#' @param samples (optional) the samples to include in the plot (default is all)
#' @param color_scheme (optional) the color scheme for the plot
#' @family plotting
#'
plot_paired = function(
    data,
    chain = c("paired", "alpha", "beta"),
    samples = NULL,
    color_scheme = NULL) {
  if(.is.list.only(data)) data = get_pair_stats(data, verbose = FALSE)
  chain1 = chain[1]
  data$log10Freq = log10(data$Freq)
  data_sub = data %>% filter(chain == chain1)
  if(!is.null(samples)) data_sub = data_sub %>% filter(sample_id %in% samples)
  #y_col = sym("Freq")
  #if(log_scale) y_col = sym("log10Freq")
  gg = ggplot(data_sub) + geom_col(aes(x=sample_id, y=Freq, fill = category)) + theme_classic() +
    scale_fill_manual(values = .tirtl_colors_gradient(palette=color_scheme, n=length(unique(data_sub$category)) )) +
    theme(axis.text.x = element_text(angle = 90, vjust = 0.5))
  return(gg)
}
