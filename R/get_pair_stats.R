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
#'
#' @family qc
#'
#' @returns A data frame with the number of pairs called by each algorithm for each sample.
#'
#' @examples
#' folder = system.file("extdata/SJTRC_TIRTL_seq_longitudinal", package = "TIRTLtools")
#' ts_data = load_tirtlseq(folder, meta_columns = c("marker", "timepoint", "version"), sep = "_", verbose = FALSE)
#'
#' get_pair_stats(ts_data)
#'
get_pair_stats = function(data, verbose = TRUE) {
  pair_stats = lapply(1:length(data$data), function(i) {
    x = data$data[[i]]
    sample_name = names(data$data)[i]
    if(verbose) {
      msg = paste("Calculating pairing stats for sample", i) %>% .add_newline()
      cat(msg)
    }
    df = .get_pair_stats_single(x) %>% mutate(sample_id = sample_name)
  }) %>% bind_rows()
  pair_stats = pair_stats %>% left_join(data$meta, by = "sample_id")
  return(pair_stats)
}
