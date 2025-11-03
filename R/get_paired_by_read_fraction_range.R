#' Calculate the number and fraction of single chains that were paired within
#' different read fraction ranges
#'
#' @description
#' `r lifecycle::badge('experimental')`
#'
#'
#' @param data a TIRTLseqData object
#' @param chain the TCR chain to use for read fraction (default is "beta")
#' @param cutoffs a vector of cutoffs for the read fraction ranges
#'
#' @family repertoire_analysis
get_paired_by_read_fraction_range = function(data, chain = c("beta","alpha"), cutoffs = 10^(-6:-1)) {
  chain = chain[1]
  df = lapply(1:length(data$data), function(i) {
    x = data$data[[i]]
    df_tmp = .get_paired_by_read_fraction_range_single(x, chain = chain, cutoffs = cutoffs) %>%
      mutate(sample_id = names(data$data)[i])
    return(df_tmp)
  }) %>% bind_rows()
  df = left_join(df, data$meta, by = "sample_id")
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
