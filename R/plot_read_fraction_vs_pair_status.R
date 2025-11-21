#' A point plot of read fraction vs. rank for the most frequent alpha/beta chains
#'
#' @description
#' `r lifecycle::badge('experimental')`
#' This function creates two point plots (alpha left, beta right) of the `n_max` most frequent
#' single-chains found in a sample.
#'
#' This function creates a plot similar to \code{\link{plot_pairs_with_eachother}()},
#' but does not mirror the alpha chains plot and add connections between paired alpha and beta chains.
#'
#' @param data a TIRTLseqData object
#' @param sample the sample to plot (either by number or sample id)
#' @param chain the chain to plot, "alpha", "beta", or "both" (default is "both")
#' @param n_max the number of most frequent single-chains to plot
#' @param show_num_partners whether to show the number of partners for each single-chain (default is FALSE)
#' @param color_scheme (optional) the color scheme for the plot
#'
#' @returns A ggplot object with the point plot
#'
#'
#' @family qc
#'
#' @examples
#' folder = system.file("extdata/SJTRC_TIRTL_seq_longitudinal", package = "TIRTLtools")
#' ts_data = load_tirtlseq(folder, meta_columns = c("marker", "timepoint", "version"), sep = "_", verbose = FALSE)
#'
#' plot_read_fraction_vs_pair_status(ts_data, sample = 1, n_max = 100)
#'
#' plot_read_fraction_vs_pair_status(ts_data, sample = 1, n_max = 100,
#' show_num_partners = TRUE)
#'
plot_read_fraction_vs_pair_status = function(
    data,
    sample,
    chain = c("both","beta", "alpha"),
    n_max = 100,
    show_num_partners = F,
    color_scheme = NULL
  ) {
  chain = chain[1]

  if(!"is_paired" %in% colnames(data$data[[1]]$beta)) data = identify_paired(data, verbose = FALSE)
  if(is.numeric(sample)) sample = names(data$data)[[sample]]
  data = data$data[[sample]]

  if(chain == "alpha") df = .make_df_pair_vs_rank(data$alpha, n_max, return_melted = FALSE) %>% mutate(chain = "alpha")
  if(chain == "beta") df = .make_df_pair_vs_rank(data$beta, n_max, return_melted = FALSE) %>% mutate(chain = "beta")
  if(chain == "both") {
    df_alpha = .make_df_pair_vs_rank(data$alpha, n_max, return_melted = FALSE) %>% mutate(chain = "alpha")
    df_beta = .make_df_pair_vs_rank(data$beta, n_max, return_melted = FALSE) %>% mutate(chain = "beta")
    df = bind_rows(df_alpha, df_beta)
  }

  df$label_unpaired = ifelse(df$is_paired, NA, "unpaired")
  yy = .get_log_labels_neg(df$readFraction)

  color_col = ifelse(show_num_partners, "n_paired_factor", "paired_status")
  df$n_paired_factor = factor(df$n_paired, levels = sort(unique(df$n_paired)))
  df = as_tibble(df)
  gg = ggplot() +
    #geom_step(size = 1) +
    geom_point(data = df, aes(x=rank, y=readFraction, color = !!sym(color_col), group = chain)) +
    geom_point(data = df %>% filter(!is.na(label_unpaired)), aes(x=rank, y=readFraction, color = !!sym(color_col), group = chain, shape = label_unpaired), size = 5) +
    #xlab("Clonal rank by readFraction") +
    #ylab(ylabel) +
    #geom_point(aes(shape = method, color = method)) +
    scale_y_log10(breaks=yy$brks,labels=yy$labels) +
    facet_wrap(~chain) +
    theme_classic() +
    ggtitle(sample) +
    scale_shape_manual(values = 3) +
    #scale_color_manual(values = .tirtl_colors_distinct(palette=color_scheme))
    scale_color_manual(values = .tirtl_colors_gradient(palette=color_scheme, n=length(unique(df[[color_col]]))) %>% rev() )
  return(gg)


}
