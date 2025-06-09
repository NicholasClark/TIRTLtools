
plot_clone_size_across_samples = function(
    data,
    clones,
    chain = c("beta", "alpha"),
    pseudo = 1e-6,
    group_vec = NULL,
    sum_readFraction = TRUE,
    samples=NULL
    ) {
  chain = chain[1]
  group_is_null = is.null(group_vec)
  if(group_is_null) group_vec = clones
  check1 = length(clones) == length(group_vec)
  if(!check1) stop("'group_vec' needs to be the same length as 'clones'")
  if(is.null(samples)) samples = names(data$data)
  gg_df = lapply(1:length(samples), function(i) {
    sample = samples[i]
    df_tmp = data$data[[sample]][[chain]]
    out_tmp = tibble(targetSequences = clones,
                     group = group_vec,
                     readFraction = df_tmp$readFraction[match(clones, df_tmp$targetSequences)]) %>%
      mutate(readFraction = na_to0(readFraction) + pseudo, source = sample)
  }) %>% bind_rows()
  grp = sym("targetSequences")
  if(sum_readFraction && !group_is_null) {
    gg_df = gg_df %>% group_by(group, source) %>%
    summarize(readFraction = sum(readFraction),
              sd_readFraction = sd(readFraction),
              n = n()
              ) %>%
    ungroup() %>%
    mutate(se_readFraction = sd_readFraction/sqrt(n))
    grp = sym("group")
  }
  gg = ggplot(gg_df, aes(x=source, y=readFraction)) +
    geom_path(aes(group = !!grp, color = group)) +
    geom_point(aes(group = !!grp, color = group)) +
    scale_y_log10() + theme_classic() +
    theme(axis.text.x = element_text(angle = 90, vjust = 0.5))

  return(gg)
}
