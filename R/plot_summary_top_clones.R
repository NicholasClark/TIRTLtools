

plot_clonotype_indices = function(
    #data, meta, type_column = "auto", proportion_column="auto",
  data, chain = c("paired", "alpha", "beta"),
  type_column = "auto", proportion_column="auto",
  cutoffs = 10^(1:5), group_col = NULL, label_col = "Sample", flip = FALSE,
  facet = FALSE, return_data = FALSE) {

  meta = data$meta
  data = data$data
  chain = chain[1]

  data = lapply(data, function(x) x[[chain]]) %>% setNames(names(data))

  labels = get_labels_from_col(meta, label_col)

  if(length(group_col) == 1) {
    meta$Group = meta[[group_col]]
  } else {
    meta$Group = apply(meta[, group_col], 1, paste, collapse = " | ")
  }

  props_list = calculate_proportions_list(data, type_column=type_column, proportion_column=proportion_column, return_list = TRUE)
  gg_df = lapply(1:length(props_list), function(i) {
    tmp_df = summarize_clonotype_indices_single(props_list[[i]], type_column = type_column,
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


summarize_clonotype_indices_single = function(data, type_column = "auto", proportion_column="auto", cutoffs = 10^(1:4)) {
  starts = c(1, cutoffs+1)
  ends = c(cutoffs, Inf)
  grps = paste("[",starts, ":", ends, "]", sep = "")
  data = data %>% arrange(desc(prop))
  grp_props = sapply(1:length(grps), function(i) {
    start = starts[i]
    end = ends[i]
    if(start > dim(data)[1]) return(0)
    if(end > dim(data)[1]) end = dim(data)[1]
    sum(data$prop[start:end])
  })
  tmp_df = tibble(clonotype_indices = factor(grps, levels = grps), prop = grp_props)
  tmp_df$prop = tmp_df$prop/sum(tmp_df$prop)
  return(tmp_df)
}
