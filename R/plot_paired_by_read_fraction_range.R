#' Bar plot of the fraction of paired single chains by frequency
#'
#' @description
#' `r lifecycle::badge('experimental')`
#' This function returns a bar plot of the fraction of single-chains (alpha or beta, default is beta)
#' that were paired within different frequency ranges
#' (default is ```[10^-1, 10^-2], [10^-2, 10^-3], ... , [10^-5,10^-6], [10^-6, 0]```) for each sample.
#'
#'
#' @param data a TIRTLseqData object or a data frame created using \code{\link{get_paired_by_read_fraction_range}()}
#' @param chain the TCR chain to plot (default is "beta")
#' @param cutoffs a vector of cutoffs for the read fraction ranges
#' @param freq if TRUE, plot the number of pairs, if FALSE plot the fraction paired
#' (default is FALSE, plot the fraction paired)
#' @param samples (optional) the samples to include in the plot
#' @param color_scheme (optional) the color scheme to use in the plot
#'
#' @returns a ggplot object with the bar plot of fraction of paired chains by frequency
#'
#'
#' @family qc
#' @examples
#' folder = system.file("extdata/SJTRC_TIRTL_seq_longitudinal", package = "TIRTLtools")
#' ts_data = load_tirtlseq(folder, meta_columns = c("marker", "timepoint", "version"), sep = "_", verbose = FALSE)
#'
#' plot_paired_by_read_fraction_range(ts_data, chain = "beta")
#'
plot_paired_by_read_fraction_range = function(
    data,
    chain = c("beta","alpha"),
    cutoffs = 10^(-6:-1),
    freq=FALSE,
    samples = NULL,
    color_scheme = NULL
    ) {
  chain = chain[1]
  if(.is.list.only(data)) data = get_paired_by_read_fraction_range(data, chain, cutoffs)
  #gg_df = get_paired_by_read_fraction_range(data)
  #gg_df = left_join(gg_df, data$meta)
  df = data
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
