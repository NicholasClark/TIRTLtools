
#' @param data a list created by the `diversity` function with diversity metrics for each sample
#' @param metric the diversity metric to use (e.g. shannon, simpson, etc.)
#' @param q (optional) for 'renyi' and 'hill' metrics, the order q of the diversity index
#' @param percent (optional) for 'dXX' metric, the percentage 'XX' between 0 and 100

plot_diversity = function(
    data, metric=.get_all_div_metrics(), q=2, percent=90, group_col = NULL,
    label_col = "Sample", flip = FALSE, facet = FALSE, log_scale = FALSE,
    return_data = FALSE
    ) {
  metric = metric[1]
  ## assume data is a list of diversity metrics for each sample

  vals = get_div_metric(data, metric, q=q, percent=percent)
  meta = data$meta
  gg_df = meta %>% mutate(value = vals)
  if(label_col == "Sample") {
    labels = meta[[1]]
  } else {
    labels = meta[[label_col]]
  }
  if(length(unique(labels)) != dim(meta)[1]) labels = paste(1:dim(meta)[1], labels)
  gg_df$Sample = factor(labels, levels = labels)
  y_label = .get_ylabel(metric=metric, q=q, percent=percent)
  y_label = paste(y_label, "|", data$call_args$type_column)
  plot_title = case_when(
    metric == "d50" ~ "d50 - The minimum number of types (clones)\nneeded to comprise 50 percent of the data",
    .default = ""
  )
  if(is.null(group_col)) {
    gg = ggplot(gg_df) +
      geom_col(aes(y=value, x= Sample)) +
      ylab(y_label) +
      ggtitle(plot_title) +
      theme_classic()
  } else {
    if(length(group_col) == 1) {
      gg_df$Group = gg_df[[group_col]]
    } else {
      gg_df$Group = apply(gg_df[, group_col], 1, paste, collapse = " | ")
    }
    if(facet) { ### separate panels for each group
      gg = ggplot(gg_df, aes(y=value, x= Sample)) +
        geom_col(aes(fill = Group)) +
        ylab(y_label) +
        theme_classic()
      if(length(group_col) == 1) {
        gg = gg + facet_wrap(~Group, scales = ifelse(flip, "free_y", "free_x"))
      } else {
        formula1 = paste(group_col, collapse = "+")
        facet_formula <- as.formula(paste("~", formula1, collapse = ""))
        gg = gg + ggh4x::facet_nested_wrap(facet_formula, scales = ifelse(flip, "free_y", "free_x"))
      }
    } else { ### one bar with mean of each group
      gg_summ = gg_df %>% group_by(Group) %>%
        summarize(n= n(), mean = mean(value), sd = sd(value)) %>%
        mutate( se=sd/sqrt(n))  %>%
        mutate( ic=se * qt((1-0.05)/2 + .5, n-1))
      gg = ggplot(gg_summ, aes(y=mean, x= Group)) + geom_col(aes(fill = Group)) +
        geom_errorbar( aes(ymin=mean-sd, ymax=mean+sd), width=0.4, colour="grey20", alpha=1, size=0.75) +
        ylab(y_label) +
        theme_classic()
    }
  }
  if(flip) {
    gg = gg +
      scale_y_continuous(trans = ifelse(log_scale, "log10", "identity")) +
      coord_flip()
  } else {
    gg = gg +
      scale_y_continuous(trans = ifelse(log_scale, "log10", "identity")) +
      scale_x_discrete(guide = guide_axis(angle = 90))
  }
  #gg = gg + scale_x_discrete(guide = guide_axis(angle = 90))
  #if(log_scale) gg = gg + scale_y_log10(labels = scales::label_log())
  #if(flip) gg = gg + coord_flip()
  gg = gg + ggtitle(plot_title) + theme(plot.title = element_text(hjust = 0.5))
  res_list = list(plot = gg)
  if(return_data) res_list$data = gg_df
  return(res_list)
}

get_div_metric = function(data, metric, q=NULL, percent = NULL) {
  ll= data$result
  if(!is.list(ll[[1]]$diversity[[metric]])) {
    vec = sapply(ll, function(x) x$diversity[[metric]])
  } else {
    if(metric %in% c("d50", "dXX")) {
      if(metric == "d50") {
        percent = 50
      } else if(metric == "dXX") {
        percent = .check_percent(percent)
      }
      vec = sapply(ll, function(x) {
        tbl = x$diversity[[metric]][[1]]
        tbl$n_clones[match(percent, tbl[["percent_required"]])]
      })
    } else if(metric %in% c("renyi", "hill") ) {
      vec = sapply(ll, function(x) {
        tbl = x$diversity[[metric]][[1]]
        tbl$value[match(q, tbl[["q"]])]
      })
    }
  }
  return(vec)
}

.check_percent = function(xx) {
  if(length(xx) > 1) stop("Error: 'percent' must be a vector of length one.")
  if(is.null(xx) || is.na(xx)) stop("Error: 'percent' is NULL. It must be a number between 0 and 100.")
  if(xx > 100 || xx <= 0) stop("Error: 'percent' must be a number between 0 and 100")
  return(xx)
}

.get_ylabel = function(metric, q, percent) {
  if(metric == "dXX") {
    return(paste("d", percent, sep = ""))
  } else if(metric %in% c("hill", "renyi")) {
    return(paste(metric,", q=", q, sep = ""))
  } else {
    return(metric)
  }
}
