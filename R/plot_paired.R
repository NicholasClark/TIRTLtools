#' Stacked bar plot of the number of TCRs paired by each algorithm
#'
#' @description
#' `r lifecycle::badge('experimental')`
#' This function plots the number of TCRs paired by each algorithm. It can also plot the number
#' of single-chains (alpha or beta) paired by the algorithms.
#'
#' @param data a TIRTLseqDataSet object or a data frame created using \code{\link{get_pair_stats}()}
#' @param chain the TCR chain to plot (default is "paired")
#' @param samples (optional) the samples to include in the plot (default is all)
#' @param by_method (optional) whether to bars by method (T-SHELL vs. MAD-HYPE)
#' @param color_scheme (optional) the color scheme for the plot
#' @family qc
#'
#' @returns a ggplot object with a stacked bar chart of the number of pairs called for each
#' sample by each algorithm.
#'
#' @examples
#' load_example_data(dataset = "SJTRC_longitudinal")
#'
#' plot_paired(SJTRC_longitudinal)
#'
#'
plot_paired = function(
    data,
    chain = c("paired", "alpha", "beta"),
    samples = NULL,
    by_method = TRUE,
    color_scheme = NULL) {
  if(.is.list.only(data)) data = get_pair_stats(data, by_method = by_method, verbose = FALSE)
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
