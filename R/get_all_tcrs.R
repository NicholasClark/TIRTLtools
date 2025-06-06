#' @title
#' Aggregate all of the T-cell receptors from a dataset
#'
#' @description
#' The \code{get_all_tcrs()} function aggregates the TCRs from all samples of a dataset
#' and puts them into one table.
#'
#' @details
#' A pair of TCRs is included twice in the TIRTLseq data if it is recognized by both
#' the T-SHELL and MAD-HYPE algorithms. If remove_duplicates is TRUE (default) the
#' function will only return one of these pairs of TCRs.
#'
#' @param data a TIRTLseq dataset created by \code{\link{load_tirtlseq}()}
#' @param chain the TCR chain, "alpha", "beta", or "paired" (default is paired)
#' @param remove_duplicates only return one TCR pair for TCRs recognized by both
#' the T-SHELL and MAD-HYPE algorithms (default is TRUE).
#'
#' @return
#' A dataframe including all of the TCRs in a dataset.
#'
#' @seealso \code{\link{load_tirtlseq}()}
#'
#' @export
#' @examples
#' # example code
#'
#'

get_all_tcrs = function(data, chain = c("paired", "alpha", "beta"), remove_duplicates = TRUE) {
  chain = chain[1]
  df_all = lapply(1:length(data$data), function(i) {
    sample_df = data$data[[i]][[chain]]
    sample_df = bind_cols(sample_df, data$meta[i,])
    return(sample_df)
  }) %>% bind_rows()
  if(chain == "paired" && remove_duplicates) df_all = remove_dupes_paired(df_all)
  return(df_all)
}
