#' @title
#' Plotting of read fractions for most frequent clonotypes
#'
#' @description
#' `plot_clonotype_indices()` creates a stacked bar chart containing the fraction of
#' reads for the top 10 most frequent clonotypes, the 11th to 100th most frequent
#' clonotypes, the 101st to 1000th most frequent, and so on, for each sample in the dataset.
#'
#' @details
#' Fill in later.
#'
#' @param data a dataset created by `load_tirtlseq()` and possible `filter_dataset()`
#' @param chain < fill in later >
#' @param type_column < fill in later >
#' @param proportion_column < fill in later >
#' @param cutoffs the indices used for the end of each group in the bar chart.
#' The default is 10^(1:5), i.e. the 10th most-frequent clone, the 100th most-frequent clone,
#' the 1,000th, the 10,000th, and the 100,000th.
#' @param group_col (optional) if supplied, a column of the metadata that will be used
#' to group samples
#' @param label_col (optional) labels for the samples
#' @param flip (optional) if TRUE, flip the x and y-axes (default is FALSE)
#' @param return_data (optional) if TRUE, return the data used to make the plot (default is FALSE)
#'
#' @return
#' A list with two objects:
#'
#' $plot - a stacked bar chart of relative frequencies of most-frequent clonotypes
#'
#' $data - if return_data is TRUE, the data frame used to create the plot
#'
#' @seealso [diversity()], [plot_diversity()]
#'
#' @export
#' @examples
#' # example code
#'
#'

plot_clonotype_indices = function(
  data, chain = c("beta", "alpha"),
  type_column = "auto", proportion_column="auto",
  cutoffs = 10^(1:5), group_col = NULL, label_col = "Sample", flip = FALSE,
  return_data = FALSE) {

  is_annotated = data$is_annotated
  meta = data$meta
  data2 = data$data
  chain = chain[1]
  if(chain == "paired") stop("Must use pseudobulk data (alpha or beta chain) for visualizing clonotype proportions.")
  if(!chain %in% c("alpha","beta")) stop("'chain' must be 'alpha' or 'beta'")

  data3 = lapply(data2, function(x) x[[chain]]) %>% setNames(names(data2))
  labels = get_labels_from_col(meta, label_col)

  is_paired = is.paired(data3[[1]])
  type_column = get_type_column(type_column, is_paired)
  proportion_column = get_proportion_column(proportion_column, is_paired, is_annotated)

  if(is.null(group_col)) {
    meta$Group = meta[[1]]
  } else if(length(group_col) == 1) {
    meta$Group = meta[[group_col]]
  } else {
    meta$Group = apply(meta[, group_col], 1, paste, collapse = " | ")
  }
  if(chain == "paired") data3 = remove_dupes_paired(data3)
  #props_list = calculate_proportions_list(data3, type_column=type_column, proportion_column=proportion_column, return_list = TRUE)

  gg_df = lapply(1:length(data3), function(i) {
    tmp_df = summarize_clonotype_indices_single(data3[[i]],
                                                proportion_column = proportion_column, cutoffs = cutoffs)
    tmp_df$Sample = labels[i]
    if(!is.null(group_col)) {
      tmp_df$Group = meta$Group[i]
    }
    return(tmp_df)
  }) %>% bind_rows()
  gg_df$Sample = factor(gg_df$Sample, levels = labels)
  y_label = "Proportion"
  plot_title = "Clonal proportion by indices"

  if(is.null(group_col)) {
    gg = ggplot(gg_df, aes(y=prop, x= Sample)) +
      geom_col(aes(fill = clonotype_indices),
               position = position_stack(reverse = TRUE), color="black", size = 0.25) +
      ylab(y_label) +
      ggtitle(plot_title)
  } else {
    gg = ggplot(gg_df, aes(y=prop, x= Sample)) +
      geom_col(aes(fill = clonotype_indices),
               position = position_stack(reverse = TRUE), color="black", size = 0.25)
      ylab(y_label)
    if(length(group_col) == 1) {
      gg = gg + facet_wrap(~Group, scales = ifelse(flip, "free_y", "free_x"))
    } else {
      formula1 = paste(group_col, collapse = "+")
      facet_formula <- as.formula(paste("~", formula1, collapse = ""))
      gg = gg + ggh4x::facet_nested_wrap(facet_formula, scales = ifelse(flip, "free_y", "free_x"))
    }
  }
  if(flip) {
    gg = gg +
      #scale_y_continuous(trans = ifelse(log_scale, "log10", "identity")) +
      coord_flip()
  } else {
    gg = gg +
      #scale_y_continuous(trans = ifelse(log_scale, "log10", "identity")) +
      scale_x_discrete(guide = guide_axis(angle = 90))
  }
  #gg = gg + scale_x_discrete(guide = guide_axis(angle = 90))
  #if(log_scale) gg = gg + scale_y_log10(labels = scales::label_log())
  #if(flip) gg = gg + coord_flip()
  gg = gg + ggtitle(plot_title) + theme(plot.title = element_text(hjust = 0.5))
  gg = gg + theme_classic()
  res_list = list(plot = gg)
  if(return_data) res_list$data = gg_df
  return(res_list)
}


summarize_clonotype_indices_single = function(data, proportion_column="auto", cutoffs = 10^(1:4), normalize = FALSE) {
  starts = c(1, cutoffs+1)
  ends = c(cutoffs, Inf)
  grps = paste("[",starts, ":", ends, "]", sep = "")
  data = data %>% arrange(desc(!!sym(proportion_column)))
  grp_props = sapply(1:length(grps), function(i) {
    start = starts[i]
    end = ends[i]
    if(start > dim(data)[1]) return(0)
    if(end > dim(data)[1]) end = dim(data)[1]
    sum(data[[proportion_column]][start:end])
  })
  tmp_df = tibble(clonotype_indices = factor(grps, levels = grps), prop = grp_props)
  if(normalize) tmp_df$prop = tmp_df$prop/sum(tmp_df$prop)
  return(tmp_df)
}
