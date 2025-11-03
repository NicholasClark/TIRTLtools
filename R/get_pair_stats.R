#' Count the number of single-chains that were paired by each algorithm
#'
#' @description
#' `r lifecycle::badge('experimental')`
#'
#' @param data a TIRTLseqData object
#' @param verbose (optional) whether to print progress of the function (default is TRUE).
#'
#' @family repertoire_analysis
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
