#' Line plot of the number of unpaired single-chains for the N most frequent single-chains
#'
#' @family plotting
plot_paired_vs_rank = function(
    data,
    sample = 1,
    y_axis = c("n_not_paired", "n_paired"),
    chain = c("both", "beta", "alpha"),
    n_max = 100,
    color_scheme = NULL
    ) {
  chain = chain[1]
  y_axis = y_axis[1]
  if(!"is_paired" %in% colnames(data$data[[1]]$beta)) data = identify_paired(data, verbose = FALSE)
  if(is.numeric(sample)) sample = names(data$data)[[sample]]
  data = data$data[[sample]]

  if(chain == "alpha") df_melt = .make_df_pair_vs_rank(data$alpha, n_max, value = y_axis) %>% mutate(chain = "alpha")
  if(chain == "beta") df_melt = .make_df_pair_vs_rank(data$beta, n_max, value = y_axis) %>% mutate(chain = "beta")
  if(chain == "both") {
    df_melt_alpha = .make_df_pair_vs_rank(data$alpha, n_max, value = y_axis) %>% mutate(chain = "alpha")
    df_melt_beta = .make_df_pair_vs_rank(data$beta, n_max, value = y_axis) %>% mutate(chain = "beta")
    df_melt = bind_rows(df_melt_alpha, df_melt_beta)
  }

  ylabel = ifelse(y_axis == "n_not_paired", "Number of clones not paired", "Number of clones paired")
  gg = ggplot(df_melt, aes(x=rank, y=!!sym(y_axis), color = method, shape = method)) +
    geom_step() +
    xlab("Clonal rank by readFraction") +
    ylab(ylabel) +
    #geom_point(aes(shape = method, color = method)) +
    facet_wrap(~chain) +
    theme_classic() +
    ggtitle(sample) +
    scale_color_manual(values = .tirtl_colors_distinct(palette=color_scheme))
  return(gg)
}


.make_df_pair_vs_rank = function(df, n_max = 1000, return_melted = TRUE, value = c("n_not_paired", "n_paired")) {
  value = value[1]
  if(value == "n_not_paired") { ### make dataframe of number of chains NOT PAIRED
    df = df %>%
      arrange(desc(readFraction)) %>%
      mutate(
        rank = 1:length(readFraction),
        n_not_paired_both = cumsum(!is_paired),
        n_not_paired_madhype = cumsum(!is_paired_madhype),
        n_not_paired_tshell = cumsum(!is_paired_tshell)
      )
    df_sub = df[1:n_max,]
    cols_keep = c("targetSequences", "readFraction", "rank", "n_not_paired_both", "n_not_paired_madhype", "n_not_paired_tshell")
    df_melt = df_sub[,cols_keep, with = FALSE] %>%
      reshape2::melt(id.vars = c("targetSequences", "readFraction", "rank"),
                     value.name = "n_not_paired", variable.name = "method")
  } else { ### make dataframe of number of chains PAIRED
    df = df %>%
      arrange(desc(readFraction)) %>%
      mutate(
        rank = 1:length(readFraction),
        n_paired_both = cumsum(is_paired),
        n_paired_madhype = cumsum(is_paired_madhype),
        n_paired_tshell = cumsum(is_paired_tshell),
      )
    df_sub = df[1:n_max,]
    cols_keep = c("targetSequences", "readFraction", "rank", "n_paired_both", "n_paired_madhype", "n_paired_tshell")
    df_melt = df_sub[,cols_keep, with = FALSE] %>%
      reshape2::melt(id.vars = c("targetSequences", "readFraction", "rank"),
                     value.name = "n_paired", variable.name = "method")
  }
  if(return_melted) res = df_melt
  if(!return_melted) res = df_sub
  return(res)
}
