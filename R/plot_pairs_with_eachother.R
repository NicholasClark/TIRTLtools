#' Point plot of read fraction vs. rank for the N most frequent alpha (left, mirrored)
#' and beta (right) chains with lines between alpha and beta chains indicating a
#' pair and a cross indicating an unpaired single-chain
#'
#' @family plotting
#'
plot_pairs_with_eachother = function(
    data,
    sample = 1,
    n_max = 100,
    show_num_partners = FALSE,
    color_scheme = NULL
  ) {

  if(!"is_paired" %in% colnames(data$data[[1]])) data = identify_paired(data, verbose = FALSE)
  if(is.numeric(sample)) sample = names(data$data)[[sample]]
  data = data$data[[sample]]

  if("beta_readFraction" %in% colnames(data$paired)) {
    df_pair = data$paired
  } else {
    df_pair = .add_single_chain_data_simple(data)
  }

  df_alpha = .make_df_pair_vs_rank(data$alpha, n_max, return_melted = FALSE) %>%
    mutate(chain = "alpha") %>%
    mutate(rank = -rank)
  df_beta = .make_df_pair_vs_rank(data$beta, n_max, return_melted = FALSE) %>%
    mutate(chain = "beta")
  df = bind_rows(df_alpha, df_beta)

  df$label_unpaired = ifelse(df$is_paired, NA, "unpaired")
  yy = .get_log_labels_neg(df$readFraction)

  df_pair$rank_beta = df_beta$rank[match(df_pair$beta_nuc, df_beta$targetSequences)]
  df_pair$rank_alpha = df_alpha$rank[match(df_pair$alpha_nuc, df_alpha$targetSequences)]
  df_pair = df_pair %>% filter(!is.na(rank_beta), !is.na(rank_alpha))

  df$n_paired_factor = factor(df$n_paired, levels = sort( union(0:2, unique(df$n_paired)) ) )
  color_col = ifelse(show_num_partners, "n_paired_factor", "paired_status")
  df$paired_status = factor(df$paired_status, levels = c("T-SHELL only", "MAD-HYPE only", "both", "neither") )
  gg = ggplot() +
    geom_point(data = df, aes(x=rank, y=readFraction, color = !!sym(color_col), group = chain), size = 1) +
    geom_point(data = df, aes(x=rank, y=readFraction, color = !!sym(color_col), group = chain, shape = label_unpaired), size = 5) +
    geom_segment(data = df_pair,
                 aes(x = rank_alpha, xend = rank_beta,
                     y=alpha_readFraction, yend = beta_readFraction),
                 size = 0.125) +
    #xlab("Clonal rank by readFraction") +
    #ylab(ylabel) +
    #geom_point(aes(shape = method, color = method)) +
    geom_vline(xintercept = 0, linetype = "dashed") +
    scale_y_log10(breaks=yy$brks,labels=yy$labels) +
    theme_classic() +
    ggtitle(sample) +
    scale_shape_manual(values = 3) +
    #scale_color_manual(values = .tirtl_colors_distinct(palette=color_scheme))
    scale_color_manual(values = .tirtl_colors_gradient(palette=color_scheme, n=length(levels(df[[color_col]]))) %>% set_names(levels(df[[color_col]])) )
  return(gg)
}
