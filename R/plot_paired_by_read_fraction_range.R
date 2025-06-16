### bar plot of number/percentage paired by readFraction
### i.e for clones between 0.1 and 0.01 readFraction, how many paired, and paired by each algorithm
### same for between 0.01 and 0.001, etc.

### data is a list of tirtlseq samples, created by load_tirtlseq()
plot_paired_by_read_fraction_range = function(df, freq=FALSE, samples = NULL, color_scheme = NULL) {
  #gg_df = get_paired_by_read_fraction_range(data)
  #gg_df = left_join(gg_df, data$meta)
  if(!is.null(samples)) df = df %>% filter(sample_id %in% samples)
  n=length(unique(df$range))
  if(!freq) {
    gg = ggplot(df) + geom_col(aes(x=sample_id, y=fraction_paired, fill = range), position = "dodge") + theme_classic() +
      scale_fill_manual(values = .tirtl_colors_gradient(palette=color_scheme, n=n)) +
      rotate_x_labels()
  } else {
    gg = ggplot(df) + geom_col(aes(x=sample_id, y=n_paired_sum, fill = range), position = "dodge") + theme_classic() +
      scale_fill_manual(values = .tirtl_colors_gradient(palette=color_scheme, n=n)) +
      rotate_x_labels() +
      scale_y_log10()
  }
  return(gg)
}

### data is a list of tirtlseq samples, created by load_tirtlseq()
get_paired_by_read_fraction_range = function(data, chain = c("beta","alpha"), cutoffs = 10^(-6:-1)) {
  chain = chain[1]
  df = lapply(1:length(data$data), function(i) {
    x = data$data[[i]]
    df_tmp = .get_paired_by_read_fraction_range_single(x, chain = chain, cutoffs = cutoffs) %>%
      mutate(sample_id = names(data$data)[i])
    return(df_tmp)
    }) %>% bind_rows()
  df = left_join(df, data$meta)
  return(df)
}

### data is for one sample, a list of three data frames: alpha, beta, and paired
.get_paired_by_read_fraction_range_single = function(data, chain = c("beta","alpha"), cutoffs = 10^(-6:-1)) {
  chain = chain[1]
  if(!"is_paired" %in% colnames(data$beta)) {
    data = .annotate_paired_single(data)
  }
  df = data[[chain]]

  starts = c(0, cutoffs)
  ends = c(cutoffs, 1)
  grps = paste("[",starts, ":", ends, "]", sep = "")

  cutoffs_full = c(cutoffs, 1)
  df$range = cut(
    df$readFraction,
    breaks = c(0, cutoffs_full),      # must start with 0
    include.lowest = TRUE,
    right = FALSE#,               # so intervals are [a, b)
    #labels = paste0("[", c(0, head(cutoffs, -1)), ", ", cutoffs, ")")
  )
  df_summ = df %>% group_by(range) %>% summarize(
    n_paired_sum = sum(is_paired),
    n_paired_sum_tshell = sum(is_paired_tshell),
    n_paired_sum_madhype = sum(is_paired_madhype),
    n_total = length(n_paired)
  ) %>% mutate(
    fraction_paired = n_paired_sum/n_total,
    fraction_paired_tshell = n_paired_sum_tshell/n_total,
    fraction_paired_madhype = n_paired_sum_madhype/n_total
  )
  return(df_summ)
}
