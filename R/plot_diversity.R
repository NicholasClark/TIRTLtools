#' @title
#' Plotting of clonal diversity metrics
#'
#' @description
#' The \code{plot_diversity()} plots the requested clonal diversity metric
#'
#' @details
#' This function can plot a variety of clonal diversity metrics for a dataset (richness,
#' Simpson diversity index, Shannon-Wiener index, etc.). See \code{\link{get_all_div_metrics}()} for
#' all available options. By default it return a barplot with one bar for each sample in
#' the dataset. If a grouping column (of the metadata) is supplied, then samples will be grouped and
#' bar heights will reflect the average diversity metric across the group, with error bars
#' showing the standard deviation.
#'
#' @param div a list created by the \code{\link{diversity}()} function with diversity metrics for each sample
#' @param metric the diversity metric to use (e.g. shannon, simpson, etc.)
#' @param q (optional) for 'renyi' and 'hill' metrics, the order q of the diversity index
#' @param percent (optional) for 'dXX' metric, the percentage 'XX' between 0 and 100
#' @param group_col (optional) if supplied, a column of the metadata that will be used
#' to group samples
#' @param label_col (optional) labels for the samples
#' @param flip (optional) if TRUE, flip the x and y-axes (default is FALSE)
#' @param facet (optional) if TRUE, plot with separate facets, or sub-plots for each group (default is FALSE)
#' @param log_scale (optional) if TRUE, use log-scale for the y-axis (default is FALSE)
#' @param return_data (optional) if TRUE, return the data used to make the plot (default is FALSE)
#'
#' @return
#' A list with two objects:
#'
#' \code{$plot} - a ggplot object with the plot of the requested diversity metric
#'
#' \code{$data} - if return_data is TRUE, the data frame used to make the plot
#'
#' @seealso \code{\link{diversity}()}, \code{\link{get_all_div_metrics}()}
#'
#' @export
#' @examples
#' # example code
#' # data = load_tirtlseq("your_directory/")
#' # div = diversity(data)
#' # plot_diversity(div, metric = "richness")
#'

plot_diversity = function(
    div, metric=get_all_div_metrics(), q=2, percent=90, group_col = NULL,
    label_col = "Sample", flip = FALSE, facet = FALSE, log_scale = FALSE,
    return_data = FALSE
    ) {
  metric = metric[1]
  #print(metric)
  #chain = div$call_args$chain
  #print(chain)
  ## assume div is a list of diversity metrics for each sample

  vals = get_div_metric(div, metric, q=q, percent=percent)
  meta = div$meta
  gg_df = meta %>% mutate(value = vals)
  if(label_col == "Sample") {
    labels = meta[[1]]
  } else {
    labels = meta[[label_col]]
  }
  if(length(unique(labels)) != dim(meta)[1]) labels = paste(1:dim(meta)[1], labels)
  gg_df$Sample = factor(labels, levels = labels)
  y_label = .get_ylabel(metric=metric, q=q, percent=percent)
  y_label = paste(y_label, "|", div$call_args$type_column)
  plot_title = case_when(
    metric == "simpson" ~ "The Simpson diversity index - the probability that two clones chosen at random represent the same type" %>% split_string_multiline(),
    metric == "gini" ~ "The Gini index/coefficient" %>% split_string_multiline(),
    metric == "gini.simpson" ~ "The Gini-Simpson index - the probability that two clones chosen at random represent different types" %>% split_string_multiline(),
    metric == "inv.simpson" ~ "The Inverse Simpson index - the effective number of types corresponding to the weighted arithmetic mean of proportional abundances" %>% split_string_multiline(),
    metric == "shannon" ~ "The Shannon-Wiener diversity index - the effective number of types corresponding to the weighted geometric mean of proportional abundances" %>% split_string_multiline(),
    metric == "berger.parker" ~ "The Berger-Parker index - the proportion of the most abundant type in the dataset" %>% split_string_multiline(),
    metric == "richness" ~ "The species richness - the total number of unique types observed in the data" %>% split_string_multiline(),
    metric == "d50" ~ "d50 - The minimum number of types (clones) needed to comprise 50 percent of the data" %>% split_string_multiline(),
    metric == "dXX" ~ "dXX - The minimum number of types (clones) needed to comprise XX percent of the data" %>% split_string_multiline(),
    metric == "renyi" ~ "The Renyi entropy - a generalization of Shannon entropy (Shannon-Wiener index)" %>% split_string_multiline(),
    metric == "hill" ~ "Hill numbers - the 'true diversity' or effective number of types" %>% split_string_multiline(),
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

get_div_metric = function(div, metric, q=NULL, percent = NULL) {
  ll= div$result
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
