#' Add single-chain read counts/fractions to the paired TCR data
#'
#'
#' @description
#' `r lifecycle::badge('deprecated')`
#' This function has been deprecated because its functionality has been added to the load_tirtlseq() function
#'
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
#' @family deprecated
#'
#' @export
#'
#'
add_single_chain_data = function(data, verbose = TRUE) {
  lifecycle::deprecate_warn(when = "0.2.1", what="add_single_chain_data()",
                            details = "This function has been deprecated because its functionality has been added to the load_tirtlseq() function"
  )
  return(data)
  # version = "new"
  # if(!is.null(data$is_annotated)) {
  #   if(data$is_annotated) {
  #     message("Single-chain data already added")
  #     return(data)
  #   }
  # }
  # if(version == "old") {
  #   data_tmp = lapply(1:length(data$data), function(i) {
  #     x = data$data[[i]]
  #     if(verbose) {
  #       msg = paste("Adding single-chain data to paired dataframe for sample", i) %>% .add_newline()
  #       cat(msg)
  #     }
  #     df_out = .add_single_chain_data_simple(x)
  #     return(list(alpha = x$alpha, beta=x$beta, paired = df_out))
  #   }) %>% set_names(names(data$data))
  #   data$data = data_tmp
  # } else {
  #   for(i in 1:length(data$data)) {
  #     if(verbose) {
  #       msg = paste("Adding single-chain data to paired dataframe for sample", i) %>% .add_newline()
  #       cat(msg)
  #     }
  #     data$data[[i]]$paired = .add_single_chain_data_simple(data$data[[i]])
  #   }
  # }
  # data$is_annotated = TRUE
  # return(data)
}
