
#' @param data a list created by the `diversity` function with diversity metrics for each sample
#' @param metric the diversity metric to use (e.g. shannon, simpson, etc.)

plot_diversity = function(data, metric, group_col = NULL, label_col = "Sample", flip = FALSE, facet = FALSE) {
  ## assume data is a list of diversity metrics for each sample

  vals = get_div_metric(data, metric)
  meta = data$meta
  gg_df = meta %>% mutate(value = vals)
  if(label_col == "Sample") {
    labels = meta[[1]]
  } else {
    labels = meta[[label_col]]
  }
  if(length(unique(labels)) != dim(meta)[1]) labels = paste(1:dim(meta)[1], labels)
  gg_df$Sample = factor(labels, levels = labels)

  if(is.null(group_col)) {
    gg = ggplot(gg_df) + geom_col(aes(y=value, x= Sample)) +
      ylab(metric) +
      theme_classic()
  } else {
    gg_df$Group = gg_df[[group_col]]
    if(facet) { ### separate panels for each group
      gg = ggplot(gg_df, aes(y=value, x= Sample)) + geom_col(aes(fill = Group)) +
        facet_wrap(~Group, scales = "free_y") +
        ylab(metric) +
        theme_classic()
    } else { ### one bar with mean of each group
      gg_summ = gg_df %>% group_by(Group) %>%
        summarize(n= n(), mean = mean(value), sd = sd(value)) %>%
        mutate( se=sd/sqrt(n))  %>%
        mutate( ic=se * qt((1-0.05)/2 + .5, n-1))
      gg = ggplot(gg_summ, aes(y=mean, x= Group)) + geom_col(aes(fill = Group)) +
        geom_errorbar( aes(ymin=mean-sd, ymax=mean+sd), width=0.4, colour="grey20", alpha=1, size=0.75) +
        ylab(metric) +
        theme_classic()
    }
  }
  if(flip) return(gg + coord_flip())
  return( gg + scale_x_discrete(guide = guide_axis(angle = 90)) )
}

get_div_metric = function(data, metric) {
  ll= data$result
  if(length(ll[[1]][[metric]]) == 1) {
    return( sapply(ll, function(x) x[[metric]]) )
  } else {
    stop()
  }

}
