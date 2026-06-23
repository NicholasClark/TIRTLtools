#' Count the number of pairs called by each algorithm
#'
#' @description
#' `r lifecycle::badge('experimental')`
#' This function tabulates the number of TCR pairs called for each sample by the MAD-HYPE
#' algorithm, the T-SHELL algorithm, or both. It is used in the \code{\link{plot_paired}()} function.
#'
#' @details
#' The function can count the number of alpha chains paired, or beta chains, or the number of pairs overall.
#' By default, it counts the number of pairs overall.
#'
#'
#' @param data a TIRTLseqDataSet object
#' @param verbose (optional) whether to print progress of the function (default is TRUE).
#' @param by_method (optional) whether to get stats for each pairing method
#'
#' @family qc
#'
#' @returns A data frame with the number of pairs called by each algorithm for each sample.
#'
#' @examples
#' load_example_data(dataset = "SJTRC_longitudinal")
#'
#' get_pair_stats(SJTRC_longitudinal)
#'
get_pair_stats = function(data, verbose = TRUE, by_method = TRUE) {
  pair_stats = lapply(1:length(data$data), function(i) {
    x = data$data[[i]]
    sample_name = names(data$data)[i]
    if(verbose) {
      msg = paste("Calculating pairing stats for sample", i) %>% .add_newline()
      cat(msg)
    }
    df = .get_pair_stats_single(x, by_method = by_method) %>% mutate(sample_id = sample_name)
  }) %>% bind_rows()
  pair_stats = pair_stats %>% left_join(data$meta, by = "sample_id")
  pair_stats$sample_id = factor(pair_stats$sample_id, levels = data$meta$sample_id)
  return(pair_stats)
}

## input is data for one sample/experiment
.get_pair_stats_single = function(data, by_method) {
  if(!"is_paired" %in% colnames(data$beta)) {
    data = .annotate_paired_single(data, by_method = by_method)
  }
  paired_df = data$paired_alt$paired_status %>% table() %>%
    as.data.frame.table() %>% as_tibble() %>%
    set_colnames(c("category", "Freq")) %>%
    mutate(chain = "paired")
  alpha_df = data$alpha$paired_status %>% table() %>%
    as.data.frame.table() %>% as_tibble() %>%
    set_colnames(c("category", "Freq")) %>%
    mutate(chain = "alpha")
  beta_df = data$beta$paired_status %>% table() %>%
    as.data.frame.table() %>% as_tibble() %>%
    set_colnames(c("category", "Freq"))  %>%
    mutate(chain = "beta")
  all_df = bind_rows(list(paired_df, alpha_df, beta_df))
  return(all_df)
}
