#' Add single-chain read counts/fractions to the paired TCR data
#'
#' @description
#' The \code{add_single_chain_data()} function adds read counts and proportions
#' from the single-chain pseudobulk data to the paired data frame for each sample
#' of a dataset.
#'
#' @param data a TIRTLseq dataset created by \code{\link{load_tirtlseq}()}
#' @param verbose (optional) whether to print progress of the function (default is TRUE).
#'
#' @return
#' A TIRTLseq dataset object where the paired data frames for each sample have
#' added columns for read counts and proportions from the single-chain pseudobulk data.
#'
#' @family data_processing
#'
#' @export
#' @examples
#' # example code
#'
#'
add_single_chain_data = function(data, verbose = TRUE) {
  if(!is.null(data$is_annotated)) {
    if(data$is_annotated) {
      message("Single-chain data already added")
      return(data)
    }
  }
  data_tmp = lapply(1:length(data$data), function(i) {
    x = data$data[[i]]
    if(verbose) {
      msg = paste("Adding single-chain data to paired dataframe for sample", i) %>% .add_newline()
      cat(msg)
    }
    df_out = .add_single_chain_data_simple(x)
    return(list(alpha = x$alpha, beta=x$beta, paired = df_out))
  }) %>% set_names(names(data$data))
  data$data = data_tmp
  data$is_annotated = TRUE
  return(data)
}



.add_single_chain_data_simple = function(data) {
  dt_pair = data$paired
  dtA = data$alpha
  dtB = data$beta
  colnames(dtA) = paste("alpha_", colnames(dtA), sep = "")
  colnames(dtB) = paste("beta_", colnames(dtB), sep = "")
  colnames(dtA)[1] = "alpha_nuc"
  colnames(dtB)[1] = "beta_nuc"
  if(!is.data.table(dtA)) dtA = as.data.table(dtA)
  if(!is.data.table(dtB)) dtB = as.data.table(dtB)
  if(!is.data.table(dt_pair)) dt_pair = as.data.table(dt_pair)
  dt_join1 = data.table::merge.data.table(x=dt_pair, y=dtA, by = "alpha_nuc", all.x=TRUE, all.y=FALSE)
  dt_join2 = data.table::merge.data.table(x=dt_join1, y=dtB, by = "beta_nuc", all.x=TRUE, all.y=FALSE)

  return(dt_join2)
}
