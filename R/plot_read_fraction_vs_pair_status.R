
plot_read_fraction_vs_pair_status = function(
    data,
    chain = c("both","beta", "alpha"),
    n_max = 100,
    show_num_partners = F,
    color_scheme = NULL
  ) {
  chain = chain[1]

  if(chain == "alpha") df = make_df_pair_vs_rank(data$alpha, n_max, return_melted = FALSE) %>% mutate(chain = "alpha")
  if(chain == "beta") df = make_df_pair_vs_rank(data$beta, n_max, return_melted = FALSE) %>% mutate(chain = "beta")
  if(chain == "both") {
    df_alpha = make_df_pair_vs_rank(data$alpha, n_max, return_melted = FALSE) %>% mutate(chain = "alpha")
    df_beta = make_df_pair_vs_rank(data$beta, n_max, return_melted = FALSE) %>% mutate(chain = "beta")
    df = bind_rows(df_alpha, df_beta)
  }

  df$label_unpaired = ifelse(df$is_paired, NA, "unpaired")
  yy = .get_log_labels_neg(df$readFraction)

  color_col = ifelse(show_num_partners, "n_paired_factor", "paired_status")
  df$n_paired_factor = factor(df$n_paired, levels = sort(unique(df$n_paired)))
  gg = ggplot(df, aes(x=rank, y=readFraction, color = !!sym(color_col), group = chain)) +
    #geom_step(size = 1) +
    geom_point() +
    geom_point(aes(shape = label_unpaired), size = 5) +
    #xlab("Clonal rank by readFraction") +
    #ylab(ylabel) +
    #geom_point(aes(shape = method, color = method)) +
    scale_y_log10(breaks=yy$brks,labels=yy$labels) +
    facet_wrap(~chain) +
    theme_classic() +
    scale_shape_manual(values = 3) +
    #scale_color_manual(values = .tirtl_colors_distinct(palette=color_scheme))
    scale_color_manual(values = .tirtl_colors_gradient(palette=color_scheme, n=length(unique(df[[color_col]]))) %>% rev() )
  return(gg)


}
