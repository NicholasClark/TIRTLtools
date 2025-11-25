#' Line plot of clonotype rank vs. read fraction for each sample
#'
#' @description
#' `r lifecycle::badge('experimental')`
#' \code{plot_ranks()} returns a line plot of TCR rank vs. read fraction for a given sample.
#'
#' @param data a TIRTLseqDataSet object
#' @param chain the chain to plot, "alpha" or "beta" (default is "beta")
#' @param column the variable to plot (default is "readFraction")
#' @param color_scheme (optional) the color scheme for the plot
#'
#' @returns a ggplot object with the rank plot
#'
#' @family qc
#'
#' @examples
#' folder = system.file("extdata/SJTRC_TIRTL_seq_longitudinal", package = "TIRTLtools")
#' ts_data = load_tirtlseq(folder, meta_columns = c("marker", "timepoint", "version"), sep = "_", verbose = FALSE)
#'
#' plot_ranks(ts_data, chain = "beta")
#'
#'

plot_ranks = function(
    data,
    chain=c("beta", "alpha"),
    column = "readFraction",
    color_scheme = NULL
) {
  chain = chain[1]
  if(chain=="paired") stop("This function currently only works with alpha or beta chains")
  data=data$data
  df_all = lapply(1:length(data), function(i) {
    x = data[[i]]
    df_tmp = x[[chain]]
    df_tmp$rank = rank(-df_tmp[[column]])
    df_tmp$Sample = names(data)[i]
    df_tmp$targetSequences = as.character(df_tmp$targetSequences)
    return(df_tmp)
  }) %>% bind_rows()
  xx = .get_log_labels_pos(df_all$rank)
  yy = .get_log_labels_neg(df_all[[column]])

  gg = ggplot(df_all) + geom_line(aes(x=rank, y=!!sym(column), color = Sample), linewidth = 1.5) +
    xlab("rank") +
    #ylab("readFraction") +
    theme_classic() +
    theme(plot.title = element_text(hjust = 0.5)) +
    theme(text = element_text(size = 16)) +
    scale_y_log10(breaks=yy$brks,labels=yy$labels) +
    scale_x_log10(breaks=xx$brks,labels=xx$labels)
  gg = gg + scale_color_manual(values = .tirtl_colors_distinct(palette=color_scheme))
  return(gg)
}
