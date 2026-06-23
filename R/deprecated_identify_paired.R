#' Identify which single chains were paired
#'
#' @description
#' `r lifecycle::badge('deprecated')`
#'
#' This function has been deprecated because its functionality has been added to the load_tirtlseq() function
#'
#' For each sample in the dataset, \code{identify_paired()} annotates the alpha
#' and beta pseudobulk data with the number of distinct pairs each chain is a part of
#' in the paired data as well as a TRUE/FALSE column indicating whether the chain
#' is paired with any partners.
#'
#' @param data a TIRTLseq dataset created by \code{\link{load_tirtlseq}()}
#' @param verbose (optional) whether to print progress of the function (default is TRUE).
#' @param by_method (optional) whether to get stats for each pairing method
#'
#' @return
#' A dataset similar to that created by \code{\link{load_tirtlseq}()}, but
#' with added columns \code{is_paired} and \code{n_paired} in the alpha and beta
#' pseudobulk data frames.
#'
#' \code{is_paired} is TRUE if the chain is found in the
#' paired data. \code{n_paired} is the number of distinct chains that the particular
#' chain is paired with.
#'
#' @family deprecated
#' @seealso \code{\link{load_tirtlseq}()}
#'
#' @export
#' @examples
#' # example code
#' # paired = load_tirtlseq("path_to/your_directory", sep = "_", meta_columns = c("cell_type", "timepoint"))
#' # paired = identify_paired(paired)
#'
identify_paired = function(data, verbose = TRUE, by_method=TRUE) {
  lifecycle::deprecate_warn(when = "0.2.1", what="identify_paired()",
                            details = "This function has been deprecated because its functionality has been added to the load_tirtlseq() function"
  )
  return(data)
  # data$data = lapply(1:length(data$data), function(i) {
  #   if(verbose) {
  #     msg = paste("Annotating data with pairing status by MAD-HYPE and T-SHELL algorithms for sample", i) %>% .prepend_newline()
  #     cat(msg)
  #   }
  #   x=data$data[[i]]
  #   .annotate_paired_single(x, by_method=by_method)
  #   }) %>% set_names(names(data$data))
  # return(data)
}
