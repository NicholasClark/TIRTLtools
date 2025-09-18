#' Line plot of clone read fraction across multiple samples
#'
#' @family plotting
#'
plot_clone_size_across_samples = function(
    data,
    clones,
    chain = c("beta", "alpha"),
    pseudo = 1e-6,
    group_vec = NULL,
    sum_readFraction = TRUE,
    samples=NULL,
    return_data = FALSE,
    label_zero = FALSE
    ) {
  chain = chain[1]
  group_is_null = is.null(group_vec)
  if(group_is_null) group_vec = clones
  check1 = length(clones) == length(group_vec)
  if(!check1) stop("'group_vec' needs to be the same length as 'clones'")
  if(is.null(samples)) samples = names(data$data)
  df = lapply(1:length(samples), function(i) {
    sample = samples[i]
    df_tmp = data$data[[sample]][[chain]]
    out_tmp = tibble(targetSequences = clones,
                     group = group_vec,
                     readFraction = df_tmp$readFraction[match(clones, df_tmp$targetSequences)]) %>%
      mutate(readFraction = .na_to0(readFraction), source = sample)
  }) %>% bind_rows()
  if(sum_readFraction && (!group_is_null)) {
    gg_df = df %>%
      group_by(group, source) %>%
      summarize(sum_readFraction = sum(readFraction),
                n = n()
                ) %>%
      ungroup() %>%
      mutate(
        log10_sum_readFraction = log10(sum_readFraction),
        sd_log10_sum_readFraction = sd(log10_sum_readFraction),
             ) %>%
      mutate(se_readFraction = sd_log10_sum_readFraction/sqrt(n))
    grp = sym("group")
    y_col = sym("sum_readFraction")
  } else {
    gg_df = df
    grp = sym("targetSequences")
    y_col = sym("readFraction")
  }
  log_labs_y = .get_log_labels_neg(gg_df[[as.character(y_col)]], pseudo, label_zero = label_zero, max_val = 1)
  #print(log_labs_y)
  gg_df = gg_df %>%
    mutate(source = factor(source, levels = samples)) %>%
    mutate(sample_id = source) %>%
    left_join(., data$meta)
  gg = ggplot(gg_df, aes(x=source, y=!!y_col+pseudo)) +
    geom_line(aes(group = !!grp, color = group)) +
    geom_point(aes(group = !!grp, color = group)) +
    #scale_y_log10() +
    theme_classic() +
    theme(axis.text.x = element_text(angle = 90, vjust = 0.5)) +
    scale_y_log10(breaks = log_labs_y$brks, labels = log_labs_y$labels, limits = c(min(log_labs_y$brks), max(log_labs_y$brks)) )
  res = gg
  if(return_data) {
    res = list(plot = gg)
    res$plot_data = gg_df
    res$all_data = df
  }
  return(res)
}
