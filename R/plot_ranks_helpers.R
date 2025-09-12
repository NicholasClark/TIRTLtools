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
