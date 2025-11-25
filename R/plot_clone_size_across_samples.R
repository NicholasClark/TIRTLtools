#' Line plot of clone read fraction across multiple samples
#'
#' @description
#' `r lifecycle::badge('experimental')`
#' This function creates a line plot of clone read fraction for the specified clones
#' across multiple samples. The function can color the lines by groups of clones.
#' If `sum_readFraction` is TRUE, the read fraction for each group will be summed
#' and presented in one line.
#'
#'
#' @param data a TIRTLseqDataSet object
#' @param clones a list of nucleotide sequences of TCR clones
#' @param chain the TCR chain used (default is "beta")
#' @param pseudo the value of the pseudocount added to all read fractions (default is 10^-6)
#' @param group_vec (optional) a vector of "groups" for the clones
#' @param sum_readFraction if TRUE, plot the sum of read fractions of clones in each group.
#' If FALSE, plot each clone with a separate line.
#' @param samples (optional) which samples to use in the plot (default is all samples)
#' @param return_data whether to return the data used for plotting (default is FALSE)
#' @param label_zero whether to label zero on the y-axis (default is FALSE)
#' @param show_legend whether to show the legend (default is TRUE)
#' @param log_scale (optional) if TRUE, use log-scale for the y-axis (default is FALSE)
#' @param x_var a column of metadata for grouping samples in the plot. The default is NULL,
#' which considers each sample its own group.
#'
#' @returns a ggplot object with a line plot of clone read fractions across samples.
#'
#' @family longitudinal
#' @examples
#' folder = system.file("extdata/SJTRC_TIRTL_seq_longitudinal",
#'   package = "TIRTLtools")
#' sjtrc = load_tirtlseq(folder,
#'   meta_columns = c("marker", "timepoint", "version"), sep = "_",
#'   verbose = FALSE)
#'
#' top_clones1 = sjtrc$data$cd8_tp1_v2$beta %>%
#'   dplyr::arrange(desc(readFraction)) %>%
#'   head(5) %>%
#'   magrittr::extract2("targetSequences") %>% as.character()
#' top_clones2 = sjtrc$data$cd8_tp2_v2$beta %>%
#'   dplyr::arrange(desc(readFraction)) %>%
#'   head(5) %>%
#'   magrittr::extract2("targetSequences") %>% as.character()
#'
#' plot_clone_size_across_samples(sjtrc,
#'   clones = c(top_clones1, top_clones2), chain = "beta")
#'
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
    label_zero = FALSE,
    show_legend = TRUE,
    log_scale = TRUE,
    x_var = NULL
    ) {
  chain = chain[1]
  group_is_null = is.null(group_vec)
  if(group_is_null) group_vec = clones
  check1 = length(clones) == length(group_vec)
  if(!check1) stop("'group_vec' needs to be the same length as 'clones'")
  if(is.null(samples)) samples = names(data$data)

  if("SimpleList" %in% class(data)) {

  } else {
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
      left_join(., data$meta, by = "sample_id")
    xx = sym("source")
    if(!is.null(x_var)) xx = sym(x_var)
    gg_df = as_tibble(gg_df)
    gg = ggplot(gg_df, aes(x=!!xx, y=!!y_col+pseudo)) +
      geom_line(aes(group = !!grp, color = group)) +
      geom_point(aes(group = !!grp, color = group)) +
      theme_classic() +
      theme(axis.text.x = element_text(angle = 90, vjust = 0.5))
    if(log_scale) {
      gg = gg + scale_y_log10(breaks = log_labs_y$brks, labels = log_labs_y$labels, limits = c(min(log_labs_y$brks), max(log_labs_y$brks)) )
    }
    if(!show_legend) gg = gg + theme(legend.position = "none")
    res = gg
    if(return_data) {
      res = list(plot = gg)
      res$plot_data = gg_df
      res$all_data = df
    }
    return(res)
  }
}
