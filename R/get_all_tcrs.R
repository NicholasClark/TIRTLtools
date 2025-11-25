#' Returns all of the paired TCRs from all samples in a dataset
#'
#' @description
#' `r lifecycle::badge('experimental')`
#' The \code{get_all_tcrs()} function aggregates the TCRs from all samples of a dataset
#' and puts them into one table.
#'
#' @details
#' A pair of TCRs is included twice in the TIRTLseq data if it is recognized by both
#' the T-SHELL and MAD-HYPE algorithms. If remove_duplicates is TRUE (default) the
#' function will only return one of these pairs of TCRs.
#'
#' @param data a TIRTLseqDataSet object created by \code{\link{load_tirtlseq}()}
#' @param chain the TCR chain, "alpha", "beta", or "paired" (default is paired)
#' @param remove_duplicates only return one TCR for TCRs paired by both
#' the T-SHELL and MAD-HYPE algorithms (default is TRUE).
#'
#' @returns
#' A dataframe including all of the TCRs in a dataset.
#'
#' @family qc
#' @seealso \code{\link{load_tirtlseq}()}
#'
#' @export
#' @examples
#' folder = system.file("extdata/SJTRC_TIRTL_seq_longitudinal", package = "TIRTLtools")
#' ts_data = load_tirtlseq(folder, meta_columns = c("marker", "timepoint", "version"), sep = "_", verbose = FALSE)
#'
#' get_all_tcrs(ts_data, chain = "paired")
#'
#'

get_all_tcrs = function(data, chain = c("paired", "alpha", "beta"), remove_duplicates = TRUE) {
  chain = chain[1]
  df_all = lapply(1:length(data$data), function(i) {
    sample_df = data$data[[i]][[chain]]
    if(chain == "paired" && remove_duplicates) sample_df = remove_duplicates(sample_df)
    sample_df = bind_cols(sample_df, data$meta[i,])
    return(sample_df)
  }) %>% bind_rows()
  return(df_all)
}
