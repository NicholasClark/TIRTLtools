#' @title
#' Line plot of TCR rank vs. frequency
#'
#' @description
#' \code{rank_plot()} returns a line plot of TCR rank vs. read fraction for a given sample.
#'
#' @details
#' Fill in later.
#'
#' @param param1 the first parameter
#' @param param2 the second parameter
#'
#' @return
#' A line plot (ggplot object) of TCR rank vs. read fraction
#'
#'
#' @export
#' @examples
#' # example code
#'
#'

rank_plot = function(
    data,
    chain=c("alpha", "beta"),
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
    return(df_tmp)
  }) %>% bind_rows()
  xx = .get_log_labels_pos(df_all$rank)
  yy = .get_log_labels_neg(df_all[[column]])

  gg = ggplot(df_all) + geom_line(aes(x=rank, y=!!sym(column), color = Sample), size = 1.5) +
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

.rank_plot_simple = function(df) {
  df$rank = rank(-df$readFraction)
  xx = .get_log_labels_pos(df$rank)
  yy = .get_log_labels_neg(df$readFraction)

  ggplot(df) + geom_line(aes(x=rank, y=readFraction)) +
    xlab("rank") + ylab("readFraction") +
    theme_classic() +
    theme(plot.title = element_text(hjust = 0.5)) +
    theme(text = element_text(size = 16)) +
    scale_y_log10(breaks=yy$brks,labels=yy$labels) +
    scale_x_log10(breaks=xx$brks,labels=xx$labels)

}
